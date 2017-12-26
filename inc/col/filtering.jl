# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

typealias EntropyCache Dict{Tuple{Int, Int}, Real};

const __DEFAULT_DELTA     =   0.0;
const __SENTINAL          =   -maxintfloat(Float64);
__NSFLTR(x)               =   x != __SENTINAL; 

#! Contract toward equilibrium
function contract_eq!(lat::Lattice, i::Int, j::Int, delta::Real, 
                      feq::AbstractArray{Float64, 1})
  for k=lat.n
    @inbounds lat.f[k, i, j] = feq[k] + delta * (lat.f[k, i, j] - feq[k]); 
  end
end

#! Dissipate energy flux
function contract_qx!(lat::Lattice, i::Int, j::Int, delta::Real,
                      feq::AbstractArray{Float64, 1}; weight_a::Real=0.5)
  f = sub(lat.f, :, i, j);
  qxij = qx_neq(f, feq, f-feq);
  weight_b = 1.0 - weight_a;

  lat.f[3, i, j] += sign(qxij) * delta;
  lat.f[6, i, j] -= sign(qxij) * delta * weight_a;
  lat.f[7, i, j] -= sign(qxij) * delta * weight_b;

  lat.f[1, i, j] -= sign(qxij) * delta;
  lat.f[5, i, j] += sign(qxij) * delta * weight_b;
  lat.f[8, i, j] += sign(qxij) * delta * weight_a;
end

#TODO clean these scale functions up with kernal functions;
#     especially with one to find neighborhood of non-equilibrium entropy densities
#TODO consider adding ensemble filters

#! Scale based on normalized median of neighboorhood non-equilibrium entropy
#!
#! \param   sim            Simulation object
#! \param   i              Index of node in i-direction
#! \param   j              Index of node in j-direction
#! \param   feq_f          Equilibrium particle distribution function
#! \param   noneq_entropy  Non-equilibrium entropy density
#! \param   entropy_cache  Cache of non-equilibrium entropy densities
#! \return                 Scale for entropic filtering
function scale_root_median(sim::AbstractSim, i::Int, j::Int, metric::LBXFunction,
                           feq_f::LBXFunction, noneq_entropy::Real, 
                           entropy_cache::EntropyCache)

  const ni, nj      = size(sim.msm.rho);
  nbr_densities     = Array{Float64}(sim.lat.n-1);
  nvalid            = 0; 

  for k=1:sim.lat.n-1
    const i_nbr     = i + sim.lat.c[1,k];
    const j_nbr     = j + sim.lat.c[2,k];
    const key       = (i_nbr, j_nbr);

    if haskey(entropy_cache, key)
      nvalid                 += 1;
      nbr_densities[nvalid]   = entropy_cache[key];
    else
      if i_nbr < 1 || ni < i_nbr || j_nbr < 1 || nj < j_nbr; continue; end

      feq                       = Array{Float64}(sim.lat.n);

      for p=1:sim.lat.n
        feq[p]                    = feq_f(sim.lat, sim.msm, 
                                          view(sim.msm.u, :, i_nbr, j_nbr), 
                                          i_nbr, j_nbr, p);
      end

      const f                   = view(sim.lat.f, :, i_nbr, j_nbr);
      const nbr_density_k       = metric(f, feq, f - feq);

      nvalid                   += 1;
      nbr_densities[nvalid]     = nbr_density_k;
      entropy_cache[key]        = nbr_density_k;
    end
  end

  if !haskey(entropy_cache, (i, j))
    feq                       = Array{Float64}(sim.lat.n);

    for p=1:sim.lat.n
      feq[p]                    = feq_f(sim.lat, sim.msm,
                                        view(sim.msm.u, :, i, j), i, j, p);
    end

    const f                   = view(sim.lat.f, :, i, j);
    entropy_cache[(i, j)]     = metric(f, feq, f - feq);
  end

  return if nvalid > 0 && entropy_cache[(i, j)] != 0.0
            const   med_nbr_density   =   median(nbr_densities[1:nvalid]);
            (sign(med_nbr_density) * 
             sqrt( abs(med_nbr_density) / entropy_cache[(i, j)] ));
         else
            0;
         end
end

#! Scale based on normalized median of neighboorhood non-equilibrium entropy
#!
#! \param   sim            Simulation object
#! \param   i              Index of node in i-direction
#! \param   j              Index of node in j-direction
#! \param   metric         Measure of nonequilibrium
#! \param   feq_f          Equilibrium particle distribution function
#! \param   noneq_entropy  Non-equilibrium entropy density
#! \param   entropy_cache  Cache of non-equilibrium entropy densities
#! \return                 Scale for entropic filtering
function scale_root_median(sim::FreeSurfSim, i::Int, j::Int, metric::LBXFunction,
                           feq_f::LBXFunction, noneq_entropy::Real, 
                           entropy_cache::EntropyCache)

  const ni, nj      = size(sim.msm.rho);
  nbr_densities     = Array{Float64}(sim.lat.n-1);
  nvalid            = 0; 

  for k=1:sim.lat.n-1
    const i_nbr     = i + sim.lat.c[1,k];
    const j_nbr     = j + sim.lat.c[2,k];
    const key       = (i_nbr, j_nbr);

    if haskey(entropy_cache, key)
      nvalid                 += 1;
      nbr_densities[nvalid] = entropy_cache[key];
    else
      if !(i_nbr < 1 || ni < i_nbr || j_nbr < 1 || nj < j_nbr ||
           sim.tracker.state[i_nbr, j_nbr] == GAS)

        feq                       = Array{Float64}(sim.lat.n);

        for p=1:sim.lat.n
          feq[p]                    = feq_f(sim.lat, sim.msm,
                                            view(sim.msm.u, :, i_nbr, j_nbr), 
                                            i_nbr, j_nbr, p);
        end

        const f                   = view(sim.lat.f, :, i_nbr, j_nbr);
        const nbr_density_k       = metric(f, feq, f - feq);

        nvalid                   += 1;
        nbr_densities[nvalid]     = nbr_density_k;
        entropy_cache[key]        = nbr_density_k;
      end
    end
  end

  if !haskey(entropy_cache, (i, j))
    feq                       = Array{Float64}(sim.lat.n);

    for p=1:sim.lat.n
      feq[p]                    = feq_f(sim.lat, sim.msm,
                                        view(sim.msm.u, :, i, j), 
                                        i, j, p);
    end

    const f                   = view(sim.lat.f, :, i, j);
    entropy_cache[(i, j)]     = metric(f, feq, f - feq);
  end

  return if nvalid > 0 && entropy_cache[(i, j)] != 0.0
            const   med_nbr_density   =   median(nbr_densities[1:nvalid]);
            (sign(med_nbr_density) * 
             sqrt( abs(med_nbr_density) / entropy_cache[(i, j)] ));
         else
            0;
         end
end

#! Scale based on normalized median of neighboorhood non-equilibrium entropy
#!
#! \param   sim              Simulation object
#! \param   i                Index of node in i-direction
#! \param   j                Index of node in j-direction
#! \param   metric           Measure of nonequilibrium
#! \param   noneq_densities  Non-equilibrium entropy density
#! \return                   Scale for entropic filtering
function scale_root_median(sim::AbstractSim, i::Int, j::Int,
                           metric::LBXFunction,
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
  return if nvalid > 0 && noneq_densities[i, j] != 0.0
            const   med_nbr_density   =   median(nbr_densities[1:nvalid]);
            (sign(med_nbr_density) * 
             sqrt( abs(med_nbr_density) / noneq_densities[i, j] ));
         else
            0;
         end

end

#! An Ehrenfests' regularization step returns f to equilibrium
scale_ehrenfests_step(args...) = 0.0;

# Filtering constants
const __METRIC            =   entropy_quadratic;
const __SCALE             =   scale_root_median;
const __DISS              =   contract_eq!;
const __DS                =   1e-4;
const __STDS              =   2.7;
const __FLTR_THRSH_WARN   =   0.01;

#TODO consider moving this kernal function to another file and including it
#! Kernal function for calculating uij, rhoij, f, feq, fneq
function __uij_rhoij_f_feq_fneq_kernal(lat::Lattice, msm::MultiscaleMap,
                                       feq_f::LBXFunction, i::Int, j::Int)

  @inbounds rhoij       =   msm.rho[i, j];
  @inbounds uij         =   view(msm.u, : ,i ,j);
  @inbounds const f     =   view(lat.f, :, i, j);
  feq         =   Array(Float64, lat.n); 
  fneq        =   Array(Float64, lat.n);

  for k = 1:lat.n 
    @inbounds feq[k]      =   feq_f(lat, msm, uij, i, j, k);
    @inbounds fneq[k]     =   lat.f[k,i,j] - feq[k];
  end

  return rhoij, uij, f, feq, fneq;

end

#! Filtered collision function with fixed DS threshold
type FltrFixedDSCol <: FltrColFunction
  feq_f::LBXFunction;
  inner_col_f!::ColFunction;
  metric::LBXFunction;
  scale::LBXFunction;
  diss!::LBXFunction;
  ds_threshold::Real;
  fltr_thrsh_warn::Real;

  function FltrFixedDSCol(inner_col_f!::ColFunction; metric::LBXFunction=__METRIC,
                          scale::LBXFunction=__SCALE, diss::LBXFunction=__DISS,
                          ds_threshold::Real=__DS, 
                          fltr_thrsh_warn::Real=__FLTR_THRSH_WARN)
    new(inner_col_f!.feq_f, inner_col_f!, metric, scale, diss, ds_threshold, 
        fltr_thrsh_warn);
  end

  function FltrFixedDSCol(inner_col_f!::ColFunction, feq_f::LBXFunction; 
                          metric::LBXFunction=__METRIC,
                          scale::LBXFunction=__SCALE, diss::LBXFunction=__DISS,
                          ds_threshold::Real=__DS, 
                          fltr_thrsh_warn::Real=__FLTR_THRSH_WARN)
    new(feq_f, inner_col_f!, metric, scale, diss, ds_threshold, 
        fltr_thrsh_warn);
  end
end

#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function (col_f::FltrFixedDSCol)(sim::AbstractSim, bounds::Matrix{Int64})
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

        ds          =   col_f.metric(f, feq, fneq);
        cache[key]  =   ds; # cache this expensive operation
      end

      if ds > col_f.ds_threshold
        const delta   =   col_f.scale(sim, i, j, col_f.metric, feq_f, ds, 
                                      cache);
        nfiltered    +=   1;

        if is_cached
          rhoij, uij, f, feq, fneq = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                   sim.msm,
                                                                   feq_f, 
                                                                   i, j);
        end

        col_f.diss!(sim.lat, i, j, delta, feq); 
      end

    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their " *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $(percent_filtered * 100)");
  end
end

#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function (col_f::FltrFixedDSCol)(sim::FreeSurfSim, bounds::Matrix{Int64})
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

          ds          =   col_f.metric(f, feq, fneq);
          cache[key]  =   ds; # cache this expensive operation
        end

        if ds > col_f.ds_threshold
          const delta   =   col_f.scale(sim, i, j, col_f.metric, feq_f, ds, 
                                        cache);
          nfiltered    +=   1;

          if is_cached
            rhoij, uij, f, feq, fneq = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                     sim.msm,
                                                                     feq_f, 
                                                                     i, j);
          end

          col_f.diss!(sim.lat, i, j, delta, feq); 
        end
      end

    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their " *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $(percent_filtered * 100)");
  end
end

#! Calling filtered collision function
#!
#! \param   sim           Simulation
#! \param   active_cells  Active flags for domain
function (col_f::FltrFixedDSCol)(sim::AbstractSim, 
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

        ds          =   col_f.metric(f, feq, fneq);
        cache[key]  =   ds; # cache this expensive operation
      end

      if ds > col_f.ds_threshold
        const delta   =   col_f.scale(sim, i, j, col_f.metric, feq_f, ds, 
                                      cache);
        nfiltered    +=   1;

        if is_cached
          rhoij, uij, f, feq, fneq = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                   sim.msm,
                                                                   feq_f, 
                                                                   i, j);
        end

        col_f.diss!(sim.lat, i, j, delta, feq); 
      end

    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their " *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $(percent_filtered * 100)");
  end
end

#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   active_cells  Active flags for domain
function (col_f::FltrFixedDSCol)(sim::FreeSurfSim,
              active_cells::Matrix{Bool})
  const ni, nj    =   size(sim.msm.rho);
  const feq_f     =   col_f.feq_f; # alias
  noneq_densities =   Array{Float64}(size(sim.msm.rho));
  cache           =   EntropyCache();
  nfiltered       =   0;
  ncollided       =   0;
  col_f.inner_col_f!(sim, active_cells);

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

        ds          =   col_f.metric(f, feq, fneq);
        cache[key]  =   ds; # cache this expensive operation
      end

      if ds > col_f.ds_threshold
        const delta   =   col_f.scale(sim, i, j, col_f.metric, feq_f, ds, 
                                      cache);
        nfiltered    +=   1;

        if is_cached
          rhoij, uij, f, feq, fneq = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                   sim.msm,
                                                                   feq_f, 
                                                                   i, j);
        end

        col_f.diss!(sim.lat, i, j, delta, feq); 
      end
    end

  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their " *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $(percent_filtered * 100)");
  end
end

#! Filtered collision function by standard deviation of noneq entropy
type FltrStdCol <: FltrColFunction
  feq_f::LBXFunction;
  inner_col_f!::ColFunction;
  metric::LBXFunction;
  scale::LBXFunction;
  diss!::LBXFunction;
  stds::Real;
  ds_threshold::Real;
  fltr_thrsh_warn::Real;

  function FltrStdCol(inner_col_f!::ColFunction; metric::LBXFunction=__METRIC,
                      scale::LBXFunction=__SCALE, diss::LBXFunction=__DISS,
                      stds::Real=__STDS, 
                      ds_threshold::Real=0.0, 
                      fltr_thrsh_warn::Real=__FLTR_THRSH_WARN)
    new(inner_col_f!.feq_f, inner_col_f!, metric, scale, diss, stds, ds_threshold,
        fltr_thrsh_warn);
  end

  function FltrStdCol(inner_col_f!::ColFunction, feq_f::LBXFunction; 
                      metric::LBXFunction=__METRIC,
                      scale::LBXFunction=__SCALE, diss::LBXFunction=__DISS,
                      stds::Real=__STDS, 
                      ds_threshold::Real=0.0, 
                      fltr_thrsh_warn::Real=__FLTR_THRSH_WARN)
    new(feq_f, inner_col_f!, metric, scale, diss, stds, ds_threshold,
        fltr_thrsh_warn);
  end
end

#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function (col_f::FltrStdCol)(sim::AbstractSim, bounds::Matrix{Int64})
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

      noneq_densities[i, j]     =   col_f.metric(f, feq, fneq);
    end
  end

  fltrd_noneq_densities   = filter(__NSFLTR, noneq_densities);
  const mean_neq_entropy  = mean(fltrd_noneq_densities);
  const std_neq_entropy   = std(fltrd_noneq_densities);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max
      if (noneq_densities[i, j] != __SENTINAL && 
          noneq_densities[i, j] > mean_neq_entropy + col_f.stds * std_neq_entropy
          && noneq_densities[i, j] > col_f.ds_threshold)
        const delta              =  col_f.scale(sim, i, j, col_f.metric, 
                                                noneq_densities);
        nfiltered               +=  1;

        @inbounds rhoij          =   sim.msm.rho[i, j];
        @inbounds uij            =   view(sim.msm.u, :, i, j);
        feq                      =   map(k -> feq_f(sim.lat, sim.msm, 
                                                    uij, i, j, k), 
                                         1:sim.lat.n);
        
        col_f.diss!(sim.lat, i, j, delta, feq); 
      end
    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their "       *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $(percent_filtered * 100)");
  end

end


#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function (col_f::FltrStdCol)(sim::FreeSurfSim, bounds::Matrix{Int64})
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

        noneq_densities[i, j]     =   col_f.metric(f, feq, fneq);
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
          noneq_densities[i, j] > mean_neq_entropy + col_f.stds * std_neq_entropy
          && noneq_densities[i, j] > col_f.ds_threshold)
        const delta              =  col_f.scale(sim, i, j, col_f.metric, 
                                                noneq_densities);
        nfiltered               +=  1;

        @inbounds rhoij          =   sim.msm.rho[i, j];
        @inbounds uij            =   view(sim.msm.u, :, i, j);
        feq                      =   map(k -> feq_f(sim.lat, msm, uij, i, j, k), 
                                         1:sim.lat.n);
        
        col_f.diss!(sim.lat, i, j, delta, feq); 

      end
    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their "       *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $(percent_filtered * 100)");
  end

end

#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   active_cells  Active flags for domain
function (col_f::FltrStdCol)(sim::AbstractSim, active_cells::Matrix{Bool})
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

      noneq_densities[i, j]     =   col_f.metric(f, feq, fneq);
    end
  end

  fltrd_noneq_densities   = filter(__NSFLTR, noneq_densities);
  const mean_neq_entropy  = mean(fltrd_noneq_densities);
  const std_neq_entropy   = std(fltrd_noneq_densities);

  for j=1:nj, i=1:ni
    if (noneq_densities[i, j] != __SENTINAL && 
        noneq_densities[i, j] > mean_neq_entropy + col_f.stds * std_neq_entropy
        && noneq_densities[i, j] > col_f.ds_threshold)
      const delta              =  col_f.scale(sim, i, j, col_f.metric, 
                                              noneq_densities);
      nfiltered               +=  1;

      @inbounds rhoij          =   sim.msm.rho[i, j];
      @inbounds uij            =   view(sim.msm.u, :, i, j);
      feq                      =   map(k -> feq_f(sim.lat, msm, uij, i, j, k), 
                                       1:sim.lat.n);
      
      col_f.diss!(sim.lat, i, j, delta, feq); 

    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their "       *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $(percent_filtered * 100)");
  end

end

#! Calling filtered collision function
#!
#! \param   sim     Simulation
#! \param   active_cells  Active flags for domain
function (col_f::FltrStdCol)(sim::FreeSurfSim, active_cells::Matrix{Bool})
  const ni, nj    = size(sim.msm.rho);
  const feq_f     = col_f.inner_col_f!.feq_f;
  noneq_densities = fill(__SENTINAL, size(sim.msm.rho));
  nfiltered       = 0;
  ncollided       = 0;
  col_f.inner_col_f!(sim, active_cells);

  # Calculate nonequilibrium entropy densities
  for j=1:nj, i=1:ni
    @inbounds if active_cells[i, j] && sim.tracker.state[i, j] != GAS
      ncollided += 1;
      rhoij, uij, f, feq, fneq  = __uij_rhoij_f_feq_fneq_kernal(sim.lat,
                                                                sim.msm,
                                                                feq_f, i, j);

      noneq_densities[i, j]     =   col_f.metric(f, feq, fneq);
    end
  end

  fltrd_noneq_densities   = filter(__NSFLTR, noneq_densities);
  const mean_neq_entropy  = mean(fltrd_noneq_densities);
  const std_neq_entropy   = std(fltrd_noneq_densities);

  for j=1:nj, i=1:ni
    if (noneq_densities[i, j] != __SENTINAL && 
        noneq_densities[i, j] > mean_neq_entropy + col_f.stds * std_neq_entropy
        && noneq_densities[i, j] > col_f.ds_threshold)
      const delta              =  col_f.scale(sim, i, j, col_f.metric, 
                                              noneq_densities);
      nfiltered               +=  1;

      @inbounds rhoij          =   sim.msm.rho[i, j];
      @inbounds uij            =   view(sim.msm.u, :, i, j);
      feq                      =   map(k -> feq_f(sim.lat, msm, uij, i, j, k), 
                                       1:sim.lat.n);
      
      col_f.diss!(sim.lat, i, j, delta, feq); 

    end
  end

  const percent_filtered  = nfiltered / ncollided;
  if percent_filtered > col_f.fltr_thrsh_warn
    warn("More than $(col_f.fltr_thrsh_warn * 100)% of the nodes had their "       *
         "non-equilibrium volume collapsed due to entropic filtering. This " *
         "can produce nonphysical results.");
    info("Percent collapsed: $(percent_filtered * 100)");
  end

end

#! Filtered collision function with positivity rule
type FltrPosCol <: FltrColFunction
  feq_f::LBXFunction;
  inner_col_f!::ColFunction;

  FltrPosCol(inner_col_f!::ColFunction) = new(inner_col_f!.feq_f, inner_col_f!);
end

#! Filtered collision function with positivity rule (call)
function (fpc::FltrPosCol)(sim::AbstractSim, args...)
  const ni, nj = size(sim.msm.rho);
  fpc.inner_col_f!(sim, args...);
  for j=1:nj, i=1:ni
    # find minimum value of f
    @inbounds min_f, min_k = sim.lat.f[1, i, j], 1;
    for k=2:sim.lat.n
      if sim.lat.f[k, i, j] < min_f
        @inbounds min_f = sim.lat.f[k, i, j];
        min_k = k;
      end
    end
    
    # If any collisions results in a negative f, back up until zero
    if min_f < 0
      const feq   = map(k -> fpc.feq_f(sim.lat, sim.msm, 
                                       view(sim.msm.u, :, i, j), i, j, k), 
                                       1:sim.lat.n);

      @inbounds const δ     = min_f / (feq[min_k] - min_f);
      for k=1:sim.lat.n
        @inbounds sim.lat.f[k, i, j] += δ * sim.lat.f[k, i, j] - feq[k];
      end
    end
  end
end

#! Calculate the nonequilibrium energy flux
#!
#! \param   f       Particle distributions
#! \param   f_eq    Equilibrium distributions
#! \param   f_neq   Non-equilibrium distributions
#! \return          Nonequilibrium energy flux
function qx_neq(f::AbstractArray{Float64, 1}, f_eq::AbstractArray{Float64, 1},
                f_neq::AbstractArray{Float64, 1})
  return (-2 * f_neq[1] + 2 * f_neq[3] + f_neq[5] - f_neq[6] - f_neq[7] + 
          f_neq[8]);
end

qx_neq_abs(f::AbstractArray{Float64, 1}, f_eq::AbstractArray{Float64, 1},
           f_neq::AbstractArray{Float64, 1}) = abs(qx_neq(f, f_eq, f_neq));

#! Calculate the nonequilibrium energy flux
#!
#! \param   f       Particle distributions
#! \param   f_eq    Equilibrium distributions
#! \param   f_neq   Non-equilibrium distributions
#! \return          Nonequilibrium energy flux
function qy_neq(f::AbstractArray{Float64, 1}, f_eq::AbstractArray{Float64, 1},
                f_neq::AbstractArray{Float64, 1})
  return (-2 * f_neq[2] + 2 * f_neq[4] + f_neq[5] + f_neq[6] - f_neq[7] -
          f_neq[8]);
end

qy_neq_abs(f::AbstractArray{Float64, 1}, f_eq::AbstractArray{Float64, 1},
           f_neq::AbstractArray{Float64, 1}) = abs(qy_neq(f, f_eq, f_neq));

#! Calculate the nonequilibrium energy flux
#!
#! \param   f       Particle distributions
#! \param   f_eq    Equilibrium distributions
#! \param   f_neq   Non-equilibrium distributions
#! \return          Nonequilibrium energy flux
function qmax_neq(f::AbstractArray{Float64, 1}, f_eq::AbstractArray{Float64, 1},
                  f_neq::AbstractArray{Float64, 1})
  return max(qx_neq_abs(f, f_eq, f_neq), qy_neq_abs(f, f_eq, f_neq));
end
