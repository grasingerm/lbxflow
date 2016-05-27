# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

typealias EntropyCache Dict{Tuple{Int, Int}, Real};

const __DEFAULT_DELTA     =   0.0;
const __SENTINAL          =   -maxintfloat(Float64);
__NSFLTR(x)               =   x != __SENTINAL; 

#TODO clean these scale functions up with kernal functions

#! Scale based on normalized median of neighboorhood non-equilibrium entropy
#!
#! \param   sim            Simulation object
#! \param   i              Index of node in i-direction
#! \param   j              Index of node in j-direction
#! \param   feq_f          Equilibrium particle distribution function
#! \param   noneq_entropy  Non-equilibrium entropy density
#! \param   entropy_cache  Cache of non-equilibrium entropy densities
#! \return                 Scale for entropic filtering
function scale_root_median(sim::AbstractSim, i::Int, j::Int, feq_f::LBXFunction,
                           noneq_entropy::Real, entropy_cache::EntropyCache)

  const ni, nj      = size(sim.msm.rho);
  nbr_densities     = Array{Float64}(sim.lat.n-1);
  nvalid            = 0; 

  for k=1:sim.lat.n-1
    const i_nbr     = i + sim.lat.c[1,k];
    const j_nbr     = j + sim.lat.c[2,k];
    const key       = (i_nbr, j_nbr);

    if haskey(entropy_cache, key) || nbr_densities[nvalid+1] >= 0
      nvalid                 += 1;
      nbr_densities[nvalid]   = entropy_cache[key];
    else
      if i_nbr < 1 || ni < i_nbr || j_nbr < 1 || nj < j_nbr; continue; end

      feq                       = Array{Float64}(sim.lat.n);

      for p=1:sim.lat.n
        feq[p]                    = feq_f(sim.lat, sim.msm.rho[i_nbr, j_nbr],
                                          sim.msm.u[:, i_nbr, j_nbr], p);
      end

      const f                   = sim.lat.f[:, i_nbr, j_nbr];
      const nbr_density_k       = entropy_quadratic(f, feq, f - feq);

      nvalid                   += 1;
      nbr_densities[nvalid]     = nbr_density_k;
      entropy_cache[key]        = nbr_density_k;
    end
  end

  return if nvalid > 0
            const   med_noneq_density   =   median(noneq_entropies[1:nvalid]);
            (sign(med_noneq_density) * 
             sqrt( abs(med_noneq_density) / noneq_densities[i, j] ));
         else
            0;
         end
end

#! Scale based on normalized median of neighboorhood non-equilibrium entropy
#!
#! \param   sim            Simulation object
#! \param   i              Index of node in i-direction
#! \param   j              Index of node in j-direction
#! \param   feq_f          Equilibrium particle distribution function
#! \param   noneq_entropy  Non-equilibrium entropy density
#! \param   entropy_cache  Cache of non-equilibrium entropy densities
#! \return                 Scale for entropic filtering
function scale_root_median(sim::FreeSurfSim, i::Int, j::Int, feq_f::LBXFunction,
                           noneq_entropy::Real, entropy_cache::EntropyCache)

  const ni, nj      = size(sim.msm.rho);
  nbr_densities     = Array{Float64}(sim.lat.n-1);
  nvalid            = 0; 

  for k=1:sim.lat.n-1
    const i_nbr     = i + sim.lat.c[1,k];
    const j_nbr     = j + sim.lat.c[2,k];
    const key       = (i_nbr, j_nbr);

    if haskey(entropy_cache, key) || nbr_densities[nvalid+1] >= 0
      nvalid                 += 1;
      nbr_densities[nvalid] = entropy_cache[key];
    else
      if !(i_nbr < 1 || ni < i_nbr || j_nbr < 1 || nj < j_nbr ||
           sim.tracker.state[i_nbr, j_nbr] == GAS)

        feq                       = Array{Float64}(sim.lat.n);

        for p=1:sim.lat.n
          feq[p]                    = feq_f(sim.lat, sim.msm.rho[i_nbr, j_nbr],
                                            sim.msm.u[:, i_nbr, j_nbr], p);
        end

        const f                   = sim.lat.f[:, i_nbr, j_nbr];
        const nbr_density_k       = entropy_quadratic(f, feq, f - feq);

        nvalid                   += 1;
        nbr_densities[nvalid]     = nbr_density_k;
        entropy_cache[key]        = nbr_density_k;
      end
    end
  end

  return if nvalid > 0
            const   med_noneq_density   =   median(noneq_entropies[1:nvalid]);
            (sign(med_noneq_density) * 
             sqrt( abs(med_noneq_density) / noneq_densities[i, j] ));
         else
            0;
         end
end

#! Scale based on normalized median of neighboorhood non-equilibrium entropy
#!
#! \param   sim              Simulation object
#! \param   i                Index of node in i-direction
#! \param   j                Index of node in j-direction
#! \param   noneq_densities  Non-equilibrium entropy density
#! \return                   Scale for entropic filtering
function scale_root_median(sim::AbstractSim, i::Int, j::Int, 
                           noneq_densities::Matrix{Float64})

  const ni, nj      = size(sim.msm.rho);
  nbr_densities     = Array{Float64}(sim.lat.n-1);
  nvalid            = 0; 

  for k=1:sim.lat.n-1 # last vector is a "rest" particle
    const i_nbr     = i + sim.lat.c[1,k];
    const j_nbr     = j + sim.lat.c[2,k];

    if (i_nbr < 1 || ni < i_nbr || j_nbr < 1 || nj < j_nbr ||
        noneq_densities[i_nbr, j_nbr] == __SENTINAL)
      continue;
    end

    nvalid                   += 1;
    nbr_densities[nvalid]     = noneq_densities[i_nbr, j_nbr];
  end

  @assert(noneq_densities[i, j] != __SENTINAL);
  return if nvalid > 0
            const   med_noneq_density   =   median(nbr_densities[1:nvalid]);
            (sign(med_noneq_density) * 
             sqrt( abs(med_noneq_density) / noneq_densities[i, j] ));
         else
            0;
         end

end

# Filtering constants
const __SCALE             =   scale_root_median;
const __DS                =   1e-4;
const __STDS              =   2.7;
const __FLTR_THRSH_WARN   =   0.01;

#TODO consider moving this kernal function to another file and including it
#! Kernal function for calculating uij, rhoij, f, feq, fneq
function __uij_rhoij_f_feq_fneq_kernal(lat::Lattice, msm::MultiscaleMap,
                                       feq_f::LBXFunction, i::Int, j::Int)

  rhoij       =   msm.rho[i,j];
  uij         =   msm.u[:,i,j];
  const f     =   lat.f[:,i,j];
  feq         =   Array(Float64, lat.n); 
  fneq        =   Array(Float64, lat.n);

  for k = 1:lat.n 
    @inbounds feq[k]      =   feq_f(lat, rhoij, uij, k);
    @inbounds fneq[k]     =   lat.f[k,i,j] - feq[k];
  end

  return rhoij, uij, f, feq, fneq;

end

#! Filtered collision function with fixed DS threshold
type FltrFixedDSCol <: ColFunction
  feq_f::LBXFunction;
  inner_col_f!::ColFunction;
  scale::LBXFunction;
  ds_threshold::Real;
  fltr_thrsh_warn::Real;

  function FltrFixedDSCol(inner_col_f!::ColFunction; scale::LBXFunction=__SCALE,
                          ds_threshold::Real=__DS, 
                          fltr_thrsh_warn::Real=__FLTR_THRSH_WARN)
    new(inner_col_f.feq_f, inner_col_f!, ds_threshold, fltr_thrsh_warn);
  end
end

#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function call(col_f::FltrFixedDSCol, sim::AbstractSim, bounds::Matrix{Int64})
  const nbounds   =   size(bounds, 2);
  const feq_f     =   col_f.feq_f; # alias
  noneq_densities =   Array{Float64}(size(sim.msm.rho));
  cache           =   EntropyCache();
  nfiltered       =   0;
  ncollided       =   0;
  col_f.inner_col_f!(sim, bounds);

  for r = 1:nbounds
    @inbounds i_min, i_max, j_min, j_max = bounds[:, r];
    ncollided += (i_max - i_min + 1) * (j_max - j_min + 1);
    for j = j_min:j_max, i = i_min:i_max

      const key = (i, j);
      is_cached = false;

      if haskey(cache, key)
        ds        = cache[key];
        is_cached = true;
      else
        rhoij, uij, f, feq, fneq = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                 sim.msm,
                                                                 feq_f, i, j);

        ds          =   entropy_quadratic(f, feq, fneq);
        cache[key]  =   ds; # cache this expensive operation
      end

      if ds > ds_threshold
        const delta   =   col_f.scale(sim, i, j, feq_f, ds, cache);
        nfiltered    +=   1;

        if is_cached
          rhoij, uij, f, feq, fneq = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                   sim.msm,
                                                                   feq_f, 
                                                                   i, j);
        end

        @simd for k = 1:sim.lat.n
          @inbounds sim.lat.f[k, i, j] = feq[k] + delta * fneq[k]; 
        end 
      end

    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their " *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $percent_filtered");
  end
end

#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function call(col_f::FltrFixedDSCol, sim::FreeSurfSim, bounds::Matrix{Int64})
  const nbounds   =   size(bounds, 2);
  const feq_f     =   col_f.feq_f; # alias
  noneq_densities =   Array{Float64}(size(sim.msm.rho));
  cache           =   EntropyCache();
  nfiltered       =   0;
  ncollided       =   0;
  col_f.inner_col_f!(sim, bounds);

  for r = 1:nbounds
    @inbounds i_min, i_max, j_min, j_max = bounds[:, r];
    for j = j_min:j_max, i = i_min:i_max
      @inbounds if sim.tracker.state[i, j] != GAS
        ncollided += 1;
        const key = (i, j);
        is_cached = false;

        if haskey(cache, key)
          ds        = cache[key];
          is_cached = true;
        else
          rhoij, uij, f, feq, fneq = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                   sim.msm,
                                                                   feq_f, i, j);

          ds          =   entropy_quadratic(f, feq, fneq);
          cache[key]  =   ds; # cache this expensive operation
        end

        if ds > ds_threshold
          const delta   =   col_f.scale(sim, i, j, feq_f, ds, cache);
          nfiltered    +=   1;

          if is_cached
            rhoij, uij, f, feq, fneq = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                     sim.msm,
                                                                     feq_f, 
                                                                     i, j);
          end

          @simd for k = 1:sim.lat.n
            @inbounds sim.lat.f[k, i, j] = feq[k] + delta * fneq[k]; 
          end 
        end
      end

    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their " *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $percent_filtered");
  end
end

#! Calling filtered collision function
#!
#! \param   sim           Simulation
#! \param   active_cells  Active flags for domain
function call(col_f::FltrFixedDSCol, sim::AbstractSim, 
              active_cells::Matrix{Bool})
  const ni, nj    =   size(sim.msm.rho);
  const feq_f     =   col_f.feq_f; # alias
  noneq_densities =   Array{Float64}(size(sim.msm.rho));
  cache           =   EntropyCache();
  nfiltered       =   0;
  ncollided       =   0;
  col_f.inner_col_f!(sim, bounds);

  for j = 1:nj, i = 1:ni

    @inbounds if active_cells[i, j]
      ncollided += 1;

      const key = (i, j);
      is_cached = false;

      if haskey(cache, key)
        ds        = cache[key];
        is_cached = true;
      else
        rhoij, uij, f, feq, fneq = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                 sim.msm,
                                                                 feq_f, i, j);

        ds          =   entropy_quadratic(f, feq, fneq);
        cache[key]  =   ds; # cache this expensive operation
      end

      if ds > ds_threshold
        const delta   =   col_f.scale(sim, i, j, feq_f, ds, cache);
        nfiltered    +=   1;

        if is_cached
          rhoij, uij, f, feq, fneq = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                   sim.msm,
                                                                   feq_f, 
                                                                   i, j);
        end

        @simd for k = 1:sim.lat.n
          @inbounds sim.lat.f[k, i, j] = feq[k] + delta * fneq[k]; 
        end 
      end

    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their " *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $percent_filtered");
  end
end

#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   active_cells  Active flags for domain
function call(col_f::FltrFixedDSCol, sim::FreeSurfSim,
              active_cells::Matrix{Bool})
  const ni, nj    =   size(sim.msm.rho);
  const feq_f     =   col_f.feq_f; # alias
  noneq_densities =   Array{Float64}(size(sim.msm.rho));
  cache           =   EntropyCache();
  nfiltered       =   0;
  ncollided       =   0;
  col_f.inner_col_f!(sim, bounds);

  for j = 1:nj, i = 1:ni
    @inbounds if active_cells[i, j] && sim.tracker.state[i, j] != GAS
      ncollided += 1;
      const key = (i, j);
      is_cached = false;

      if haskey(cache, key)
        ds        = cache[key];
        is_cached = true;
      else
        rhoij, uij, f, feq, fneq = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                 sim.msm,
                                                                 feq_f, i, j);

        ds          =   entropy_quadratic(f, feq, fneq);
        cache[key]  =   ds; # cache this expensive operation
      end

      if ds > ds_threshold
        const delta   =   col_f.scale(sim, i, j, feq_f, ds, cache);
        nfiltered    +=   1;

        if is_cached
          rhoij, uij, f, feq, fneq = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                   sim.msm,
                                                                   feq_f, 
                                                                   i, j);
        end

        @simd for k = 1:sim.lat.n
          @inbounds sim.lat.f[k, i, j] = feq[k] + delta * fneq[k]; 
        end 
      end
    end

  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their " *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $percent_filtered");
  end
end

#! Filtered collision function by standard deviation of noneq entropy
type FltrStdCol <: ColFunction
  feq_f::LBXFunction;
  inner_col_f!::ColFunction;
  scale::LBXFunction;
  stds::Real;
  fltr_thrsh_warn::Real;

  function FltrStdCol(inner_col_f!::ColFunction; scale::LBXFunction=__SCALE,
                      stds::Real=__DS, fltr_thrsh_warn::Real=__FLTR_THRSH_WARN)
    new(inner_col_f!.feq_f, inner_col_f!, scale, stds, fltr_thrsh_warn);
  end
end

#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function call(col_f::FltrStdCol, sim::AbstractSim, bounds::Matrix{Int64})
  const ni, nj    = size(sim.msm.rho);
  const nbounds   = size(bounds, 2);
  const feq_f     = col_f.inner_col_f!.feq_f;
  noneq_densities = fill(__SENTINAL, size(sim.msm.rho));
  nfiltered       = 0;
  ncollided       = 0;
  col_f.inner_col_f!(sim, bounds);

  # Calculate nonequilibrium entropy densities
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    ncollided += (i_max - i_min + 1) * (j_max - j_min + 1);
    for j = j_min:j_max, i = i_min:i_max
      rhoij, uij, f, feq, fneq  = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                sim.msm,
                                                                feq_f, i, j);

      noneq_densities[i, j]     =   entropy_quadratic(f, feq, fneq);
    end
  end

  fltrd_noneq_densities   = filter(__NSFLTR, noneq_densities);
  const mean_neq_entropy  = mean(fltrd_noneq_densities);
  const std_neq_entropy   = std(fltrd_noneq_densities);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max
      if (noneq_densities[i, j] != __SENTINAL && 
          noneq_densities[i, j] > mean_neq_entropy + col_f.stds * std_neq_entropy)
        const delta              =  col_f.scale(sim, i, j, noneq_densities);
        nfiltered               +=  1;

        rhoij                    =   sim.msm.rho[i,j];
        uij                      =   sub(sim.msm.u, :, i, j);

        for k = 1:sim.lat.n 
          const feq                 =   feq_f(sim.lat, rhoij, uij, k);
          const fneq                =   sim.lat.f[k, i, j] - feq;
          @inbounds sim.lat.f[k, i, j]  =   feq + delta * fneq; 
        end

      end
    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their "       *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $percent_filtered");
  end

end


#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function call(col_f::FltrStdCol, sim::FreeSurfSim, bounds::Matrix{Int64})
  const ni, nj    = size(sim.msm.rho);
  const nbounds   = size(bounds, 2);
  const feq_f     = col_f.inner_col_f!.feq_f;
  noneq_densities = fill(__SENTINAL, size(sim.msm.rho));
  nfiltered       = 0;
  ncollided       = 0;
  col_f.inner_col_f!(sim, bounds);

  # Calculate nonequilibrium entropy densities
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max
      @inbounds if sim.tracker.state[i, j] != GAS
        ncollided += 1;
        rhoij, uij, f, feq, fneq  = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                  sim.msm,
                                                                  feq_f, i, j);

        noneq_densities[i, j]     =   entropy_quadratic(f, feq, fneq);
      end
    end
  end

  fltrd_noneq_densities   = filter(__NSFLTR, noneq_densities);
  const mean_neq_entropy  = mean(fltrd_noneq_densities);
  const std_neq_entropy   = std(fltrd_noneq_densities);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max
      if (noneq_densities[i, j] != __SENTINAL && 
          noneq_densities[i, j] > mean_neq_entropy + col_f.stds * std_neq_entropy)
        const delta              =  col_f.scale(sim, i, j, noneq_densities);
        nfiltered               +=  1;

        rhoij                    =   sim.msm.rho[i,j];
        uij                      =   sub(sim.msm.u, :, i, j);

        for k = 1:sim.lat.n 
          const feq                 =   feq_f(sim.lat, rhoij, uij, k);
          const fneq                =   sim.lat.f[k, i, j] - feq;
          @inbounds sim.lat.f[k, i, j]  =   feq + delta * fneq; 
        end

      end
    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their "       *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $percent_filtered");
  end

end

#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   active_cells  Active flags for domain
function call(col_f::FltrStdCol, sim::AbstractSim, active_cells::Matrix{Bool})
  const ni, nj    = size(sim.msm.rho);
  const feq_f     = col_f.inner_col_f!.feq_f;
  noneq_densities = fill(__SENTINAL, size(sim.msm.rho));
  nfiltered       = 0;
  ncollided       = 0;
  col_f.inner_col_f!(sim, bounds);

  # Calculate nonequilibrium entropy densities
  for j=1:nj, i=1:ni
    @inbounds if active_cells[i, j]
      ncollided += 1;
      rhoij, uij, f, feq, fneq  = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                sim.msm,
                                                                feq_f, i, j);

      noneq_densities[i, j]     =   entropy_quadratic(f, feq, fneq);
    end
  end

  fltrd_noneq_densities   = filter(__NSFLTR, noneq_densities);
  const mean_neq_entropy  = mean(fltrd_noneq_densities);
  const std_neq_entropy   = std(fltrd_noneq_densities);

  for j=1:nj, i=1:ni
    if (noneq_densities[i, j] != __SENTINAL && 
        noneq_densities[i, j] > mean_neq_entropy + col_f.stds * std_neq_entropy)
      const delta              =  col_f.scale(sim, i, j, noneq_densities);
      nfiltered               +=  1;

      rhoij                    =   sim.msm.rho[i,j];
      uij                      =   sub(sim.msm.u, :, i, j);

      for k = 1:sim.lat.n 
        const feq                 =   feq_f(sim.lat, rhoij, uij, k);
        const fneq                =   sim.lat.f[k, i, j] - feq;
        @inbounds sim.lat.f[k, i, j]  =   feq + delta * fneq; 
      end

    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their "       *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $percent_filtered");
  end

end

#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   active_cells  Active flags for domain
function call(col_f::FltrStdCol, sim::FreeSurfSim, active_cells::Matrix{Bool})
  const ni, nj    = size(sim.msm.rho);
  const feq_f     = col_f.inner_col_f!.feq_f;
  noneq_densities = fill(__SENTINAL, size(sim.msm.rho));
  nfiltered       = 0;
  ncollided       = 0;
  col_f.inner_col_f!(sim, bounds);

  # Calculate nonequilibrium entropy densities
  for j=1:nj, i=1:ni
    @inbounds if active_cells[i, j] && sim.tracker.state[i, j] != GAS
      ncollided += 1;
      rhoij, uij, f, feq, fneq  = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                sim.msm,
                                                                feq_f, i, j);

      noneq_densities[i, j]     =   entropy_quadratic(f, feq, fneq);
    end
  end

  fltrd_noneq_densities   = filter(__NSFLTR, noneq_densities);
  const mean_neq_entropy  = mean(fltrd_noneq_densities);
  const std_neq_entropy   = std(fltrd_noneq_densities);

  for j=1:nj, i=1:ni
    if (noneq_densities[i, j] != __SENTINAL && 
        noneq_densities[i, j] > mean_neq_entropy + col_f.stds * std_neq_entropy)
      const delta              =  col_f.scale(sim, i, j, noneq_densities);
      nfiltered               +=  1;

      rhoij                    =   sim.msm.rho[i,j];
      uij                      =   sub(sim.msm.u, :, i, j);

      for k = 1:sim.lat.n 
        const feq                 =   feq_f(sim.lat, rhoij, uij, k);
        const fneq                =   sim.lat.f[k, i, j] - feq;
        @inbounds sim.lat.f[k, i, j]  =   feq + delta * fneq; 
      end

    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their "       *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $percent_filtered");
  end

end
