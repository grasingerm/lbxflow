# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

include("constitutive.jl");
include("forcing.jl");
include("equilibrium.jl");
include(joinpath("..", "lattice.jl"));
include("mrt_matrices.jl");
include(joinpath("..", "multiscale.jl"));
include(joinpath("..", "numerics.jl"));

typealias EntropyCache Dict{Tuple{Int, Int}, Real};

#! Scale based on normalized median of neighboorhood non-equilibrium entropy
#!
#! \param   sim            Simulation object
#! \param   i              Index of node in i-direction
#! \param   j              Index of node in j-direction
#! \param   feq_f          Equilibrium particle distribution function
#! \param   noneq_entropy  Non-equilibrium entropy density
#! \param   entropy_cache  Cache of non-equilibrium entropy densities
#! \return                 Scale for entropic filtering
function scale_root_median(sim::Sim, i::Int, j::Int, feq_f::Function,
                           noneq_entropy::Real, entropy_cache::EntropyCache)

  const ni, nj      = size(sim.msm.rho);
  noneq_entropies   = Array{Float64}(sim.lat.n-1);
  nvalid            = 0; 

  for k=1:sim.lat.n-1
    const i_nbr     = i + sim.lat.c[1,k];
    const j_nbr     = j + sim.lat.c[2,k];
    const key       = (i_nbr, j_nbr);

    if haskey(entropy_cache, key)
      nvalid                 += 1;
      noneq_entropies[nvalid] = entropy_cache[key];
    else
      if i_nbr < 1 || ni < i_nbr || j_nbr < 1 || nj < j_nbr; continue; end

      nvalid                 += 1;
      feq                     = Array{Float64}(sim.lat.n);

      for p=1:sim.lat.n
        feq[p]                    = feq_f(sim.lat, sim.msm.rho[i_nbr, j_nbr],
                                          sim.msm.u[:, i_nbr, j_nbr], p);
      end
        const f                   = sim.lat.f[:, i_nbr, j_nbr];
        const noneq_entropy_k     = entropy_quadratic(f, feq, f - feq);
        noneq_entropies[nvalid]   = noneq_entropy_k;
        entropy_cache[key]        = noneq_entropy_k;
    end
  end

  return sqrt( median(noneq_entropies[1:nvalid]) / noneq_entropy ); 
end

#! Scale based on normalized median of neighboorhood non-equilibrium entropy
#!
#! \param   sim              Simulation object
#! \param   i                Index of node in i-direction
#! \param   j                Index of node in j-direction
#! \param   noneq_entropies  Non-equilibrium entropy density
#! \return                   Scale for entropic filtering
function scale_root_median(sim::Sim, i::Int, j::Int, 
                           noneq_densities::Matrix{Float64})

  const ni, nj      = size(sim.msm.rho);
  noneq_entropies   = Array{Float64}(sim.lat.n-1);
  nvalid            = 0; 

  for k=1:sim.lat.n-1 # last vector is a "rest" particle
    const i_nbr     = i + sim.lat.c[1,k];
    const j_nbr     = j + sim.lat.c[2,k];

    if i_nbr < 1 || ni < i_nbr || j_nbr < 1 || nj < j_nbr; continue; end

    nvalid                   += 1;
    noneq_entropies[nvalid]   = noneq_densities[i_nbr, j_nbr];
  end

  return sqrt( median(noneq_entropies[1:nvalid]) / noneq_densities[i, j] ); 
end

# Filtering constants
const __FILTER    = scale_root_median;
const __DS        = 1e-6;
const __STDS      = 2.7;


#TODO should collision functions have their own type ...
#     (with equilibrium function as data member)?

#! Collision function with entropic filtering wrapped around it
#!
#! \param   inner_col_f   Inner collision function
#! \param   scale         Scale function
#! \param   feq_f         Equilibrium particle distribution function
#! \param   ds_threshold  Nonequilibrium density threshold for filtering 
function init_col_filter_cache_fixds(inner_col_f!::LBXFunction, 
                                     scale::LBXFunction=__FILTER;
                                     feq_f::LBXFunction=feq_incomp,
                                     ds_threshold::Real=__DS)
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    const nbounds   = size(bounds, 2);
    noneq_densities = Array{Float64}(size(sim.msm.rho));
    cache           = EntropyCache();
    nfiltered       = 0;
    ncollided       = 0;
    inner_col_f!(sim, bounds);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      ncollided += (i_max - i_min) * (j_max - j_min);
      for j = j_min:j_max, i = i_min:i_max
        const key = (i, j);
        if haskey(cache, key)
          ds      = cache[key];
        else

          rhoij       =   sim.msm.rho[i,j];
          uij         =   sim.msm.u[:,i,j];
          const f     =   sim.lat.f[:,i,j];
          feq         =   Array(Float64, sim.lat.n); 
          fneq        =   Array(Float64, sim.lat.n);

          for k = 1:sim.lat.n 
            feq[k]      =   feq_f(sim.lat, rhoij, uij, k);
            fneq[k]     =   sim.lat.f[k,i,j] - feq[k];
          end

          ds          =   entropy_quadratic(f, feq, fneq);
          cache[key]  =   ds; # cache this expensive operation
        end

        if ds > ds_threshold
          const delta   =   scale(sim, i, j, feq_f, ds, cache);
          nfiltered    +=   1;

          for k = 1:sim.lat.n
            sim.lat.f[k, i, j] = feq[k] + delta * fneq[k]; 
          end 
        end

      end
    end

    const percent_filtered  = nfiltered / ncollided;
    if percent_filtered > 0.01
      warn("More than 1% of the nodes had their non-equilibrium volume " *
           "collapsed due to entropic filtering. This can produce nonphysical "*
           "results.");
      info("Percent collapsed: $percent_filtered");
    end

  end
end

#! Collision function with entropic filtering wrapped around it
#!
#! \param   inner_col_f   Inner collision function
#! \param   scale         Scale function
#! \param   feq_f         Equilibrium particle distribution function
#! \param   stds          Standard deviations away from mean in which to filter
function init_col_filter_std(inner_col_f!::LBXFunction, 
                             scale::LBXFunction=__FILTER;
                             feq_f::LBXFunction=feq_incomp,
                             stds::Real=__STDS)
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    const ni, nj    = size(sim.msm.rho);
    const nbounds   = size(bounds, 2);
    noneq_densities = Array{Float64}(size(sim.msm.rho));
    nfiltered       = 0;
    ncollided       = 0;
    inner_col_f!(sim, bounds);

    # Calculate nonequilibrium entropy densities
    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      ncollided += (i_max - i_min) * (j_max - j_min);
      for j = j_min:j_max, i = i_min:i_max

        rhoij                 =   sim.msm.rho[i,j];
        uij                   =   sim.msm.u[:,i,j];
        const f               =   sim.lat.f[:,i,j];
        feq                   =   Array(Float64, sim.lat.n); 
        fneq                  =   Array(Float64, sim.lat.n);

        #TODO consider caching these values... but only 1%ish need recalced?
        for k = 1:sim.lat.n 
          feq[k]                  =   feq_f(sim.lat, rhoij, uij, k);
          fneq[k]                 =   sim.lat.f[k,i,j] - feq[k];
        end

        noneq_densities[i, j] =   entropy_quadratic(f, feq, fneq);
      end

    end

    const mean_neq_entropy  = mean(noneq_densities);
    const std_neq_entropy   = std(noneq_densities);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        if noneq_densities[i, j] > mean_neq_entropy + stds * std_neq_entropy
          const delta              =  scale(sim, i, j, noneq_densities);
          nfiltered               +=  1;

          rhoij                    =   sim.msm.rho[i,j];
          uij                      =   sim.msm.u[:,i,j];

          for k = 1:sim.lat.n 
            const feq                 =   feq_f(sim.lat, rhoij, uij, k);
            const fneq                =   sim.lat.f[k, i, j] - feq;
            sim.lat.f[k, i, j]        =   feq + delta * fneq; 
          end

        end
      end
    end

    const percent_filtered  = nfiltered / ncollided;
    if percent_filtered > 0.01
      warn("More than 1% of the nodes had their non-equilibrium volume " *
           "collapsed due to entropic filtering. This can produce nonphysical "*
           "results.");
      info("Percent collapsed: $percent_filtered");
    end

  end
end
