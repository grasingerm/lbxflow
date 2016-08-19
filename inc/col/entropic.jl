# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import Roots;

const __KAPPA_ERROR_MSG__ = ("Entropic relaxation factor, kappa, must be "  *
                             "greater than or equal to zero and less than " *
                             "or equal to one for the collision operator "  *
                             "to Lyapunov stable."); 

#! Entropically stabilized LBGK collision operator for incompressible flow
#!
#! \param constit_relation_f  Constitutive relationship
#! \param feq_f               Equilibrium particle distribution function
#! \param kappa               Entropic relaxation factor
#! \param eps_ds              Change in entropy threshold for ELBM
#! \param search_entropic_stabiliy Search function for stability factor, alpha
#! \return collision_function!(sim, bounds)
function init_col_entropic_srt(constit_relation_f::Function; 
                               feq_f::Function=feq_incomp_max_entropy,
                               kappa::Real=__KAPPA, eps_ds::Real=__EPS_DS,
                               search_entropic_stability=__SES)
  @assert(0.0 <= kappa && kappa <= 1.0, __KAPPA_ERROR_MSG__);

  return (sim::AbstractSim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    root_not_found = false;

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max

        rhoij       =   msm.rho[i,j];
        uij         =   msm.u[:,i,j];
        feq         =   Array(Float64, lat.n); 
        fneq        =   Array(Float64, lat.n); 
        const f     =   lat.f[:,i,j];

        for k = 1:lat.n 
          feq[k]    = feq_f(lat, rhoij, uij, k);
          fneq[k]   = lat.f[k,i,j] - feq[k];
        end

        const mu      =   constit_relation_f(sim, fneq, i, j);
        const omega   =   @omega(mu, lat.cssq, lat.dt);

        if (entropy_lat_boltzmann(lat, feq) - entropy_lat_boltzmann(lat, f)
            < eps_ds) # do not to enforce entropic stability, use standard BGK
          for k = 1:lat.n
            lat.f[k,i,j] += omega * (feq[k] - lat.f[k,i,j]);
          end
        else
          const alpha    = search_entropic_stability(lat, f, feq);
          if alpha < __ALPHA_EPS
            root_not_found = true;
            for k = 1:lat.n
              lat.f[k,i,j] += omega * (feq[k] - lat.f[k,i,j]);
            end
          else
            for k = 1:lat.n
              lat.f[k,i,j] += kappa * alpha * omega * (feq[k] - lat.f[k,i,j]);
            end
          end
        end

        msm.omega[i,j] = omega;
      end
    end

    if root_not_found
      warn("Root for entropic stabilization not found.");
      warn("Solution may become unstable due to decrease in entropy.");
    end

  end
end

#! Entropically stabilized LBGK collision operator for incompressible flow
#!
#! \param constit_relation_f  Constitutive relationship
#! \param forcing_kf          Forcing functions
#! \param feq_f               Equilibrium particle distribution function
#! \param kappa               Entropic relaxation factor
#! \param eps_ds              Change in entropy threshold for ELBM
#! \param search_entropic_stabiliy Search function for stability factor, alpha
#! \return collision_function!(sim, bounds)
function init_col_entropic_srt(constit_relation_f::Function,
                               forcing_kf::Tuple{Function, Function};
                               feq_f::Function=feq_incomp_max_entropy,
                               kappa::Real=__KAPPA, eps_ds::Real=__EPS_DS,
                               search_entropic_stability=__SES)
  @assert(0.0 <= kappa && kappa <= 1.0, __KAPPA_ERROR_MSG__);

  const uf, colf = forcing_kf;
  return (sim::AbstractSim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    root_not_found = false;

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max

        rhoij       =   msm.rho[i,j];
        uij         =   uf(lat, msm.u[:,i,j]);
        feq         =   Array(Float64, lat.n); 
        fneq        =   Array(Float64, lat.n); 
        const f     =   lat.f[:,i,j];

        for k = 1:lat.n 
          feq[k]  = feq_f(lat, rhoij, uij, k);
          fneq[k] = lat.f[k,i,j] - feq[k];
        end

        const mu      =   constit_relation_f(sim, fneq, i, j);
        const omega   =   @omega(mu, lat.cssq, lat.dt);

        if (entropy_lat_boltzmann(lat, feq) - entropy_lat_boltzmann(lat, f)
            < eps_ds) # do not to enforce entropic stability, use standard BGK
          for k = 1:lat.n
            lat.f[k,i,j] += (omega * (feq[k] - lat.f[k,i,j]) 
                             + colf(lat, omega, uij, k));
          end
        else
          const alpha    = search_entropic_stability(lat, f, feq);
          if alpha < __ALPHA_EPS
            root_not_found = true;
            for k = 1:lat.n
              lat.f[k,i,j] += (omega * (feq[k] - lat.f[k,i,j])
                               + colf(lat, omega, uij, k));
            end
          else
            for k = 1:lat.n
              lat.f[k,i,j] += (kappa * alpha * omega * (feq[k] - lat.f[k,i,j])
                               + colf(lat, omega, uij, k));
            end
          end
        end
        msm.omega[i,j] = omega;
      end
    end

    if root_not_found
      warn("Root for entropic stabilization not found.");
      warn("Solution may become unstable due to decrease in entropy.");
    end

  end
end


#! Entropically stabilized MRT collision operator for incompressible flow
#!
#! \param constit_relation_f  Constitutive relationship
#! \param feq_f               Equilibrium particle distribution function
#! \param S                   Function for building a relaxation matrix
#! \param kappa               Entropic relaxation factor
#! \param eps_ds              Change in entropy threshold for ELBM
#! \param search_entropic_stabiliy Search function for stability factor, alpha
#! \return collision_function!(sim, bounds)
function init_col_entropic_mrt(constit_relation_f::Function;
                               feq_f::Function=feq_incomp_max_entropy,
                               S::Function=__S,
                               kappa::Real=__KAPPA, eps_ds::Real=__EPS_DS,
                               search_entropic_stability=__SES)
  @assert(0.0 <= kappa && kappa <= 1.0, __KAPPA_ERROR_MSG__);
  return (sim::AbstractSim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const M = @DEFAULT_MRT_M();
    const iM = @DEFAULT_MRT_IM();
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    root_not_found = false;

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max

        rhoij     = msm.rho[i,j];
        uij       = msm.u[:,i,j];
        feq       = Array(Float64, lat.n); 
        fneq      = Array(Float64, lat.n); 
        const f   = lat.f[:,i,j];

        for k = 1:lat.n 
          feq[k]  = feq_f(lat, rhoij, uij, k);
          fneq[k] = lat.f[k,i,j] - feq[k];
        end
        const muij = constit_relation_f(sim, fneq, S, M, iM, i, j);
        const Sij  = S(muij, rhoij, lat.cssq, lat.dt);
        if (entropy_lat_boltzmann(lat, feq) - entropy_lat_boltzmann(lat, f)
            < eps_ds) # do not to enforce entropic stability, use standard BGK
          lat.f[:,i,j] = f - iM * Sij * M * fneq; # perform collision
        else
          const alpha    = search_entropic_stability(lat, f, feq);
          if alpha < __ALPHA_EPS
            root_not_found = true;
            lat.f[:,i,j] = f - iM * Sij * M * fneq;
          else
            for a in (8,9); Sij[a] *= kappa * alpha; end;
            lat.f[:,i,j] = f - iM * Sij * M * fneq;
          end
        end

        msm.omega[i,j] = @omega(mu, lat.cssq, lat.dt);;
      end
    end

    if root_not_found
      warn("Root for entropic stabilization not found.");
      warn("Solution may become unstable due to decrease in entropy.");
    end

  end
end

#! Entropically stabilized MRT collision operator for incompressible flow
#!
#! \param constit_relation_f  Constitutive relationship
#! \param forcing_kf          Forcing functions
#! \param feq_f               Equilibrium particle distribution function
#! \param S                   Function for building a relaxation matrix
#! \param kappa               Entropic relaxation factor
#! \param eps_ds              Change in entropy threshold for ELBM
#! \param search_entropic_stabiliy Search function for stability factor, alpha
#! \return collision_function!(sim, bounds)
function init_col_entropic_mrt(constit_relation_f::Function,
                               forcing_kf::Tuple{Function, Function};
                               feq_f::Function=feq_incomp_max_entropy,
                               S::Function=__S,
                               kappa::Real=__KAPPA, eps_ds::Real=__EPS_DS,
                               search_entropic_stability=__SES)
  @assert(0.0 <= kappa && kappa <= 1.0, __KAPPA_ERROR_MSG__);
  
  const uf, colf = forcing_kf;
  return (sim::AbstractSim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const M = @DEFAULT_MRT_M();
    const iM = @DEFAULT_MRT_IM();
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    root_not_found = false;

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max

        rhoij     = msm.rho[i,j];
        uij       = uf(lat, msm.u[:,i,j]);
        feq       = Array(Float64, lat.n); 
        fneq      = Array(Float64, lat.n); 
        const f   = lat.f[:,i,j];

        for k = 1:lat.n 
          feq[k]  = feq_f(lat, rhoij, uij, k);
          fneq[k] = lat.f[k,i,j] - feq[k];
        end

        const muij    = constit_relation_f(sim, fneq, S, M, iM, i, j);
        const Sij     = S(muij, rhoij, lat.cssq, lat.dt);
        const omegaij = @omega(mu, lat.cssq, lat.dt);;
        fdl           = Array(Float64, lat.n);

        for k = 1:lat.n
          fdl[k] = colf(lat, omegaij, uij, k);
        end

        if (entropy_lat_boltzmann(lat, feq) - entropy_lat_boltzmann(lat, f)
            < eps_ds) # do not to enforce entropic stability, use standard BGK
          lat.f[:,i,j]   = f - iM * Sij * M * fneq + fdl;
        else
          const alpha    = search_entropic_stability(lat, f, feq);
          if alpha < __ALPHA_EPS
            root_not_found = true;
            lat.f[:,i,j] = f - iM * Sij * M * fneq + fdl;
          else
            for a in (8,9); Sij[a] *= kappa * alpha; end;
            lat.f[:,i,j] = f - iM * Sij * M * fneq + fdl;
          end
        end

        msm.omega[i,j] = omegaij;
      end
    end

    if root_not_found
      warn("Root for entropic stabilization not found.");
      warn("Solution may become unstable due to decrease in entropy.");
    end

  end
end

#TODO should we allow the user to specify their own H-function?
#! Search for alpha, the limit of over-relaxation for entropic involution
#!
#! \param   lat       Lattice
#! \param   f         Particle distribution vector
#! \param   feq       Equilibrium particle distribution vector
#! \param   sbounds   Search bounds, default (0.0, 2.0)
#! \return            Limit of over-relaxation for entropic involution
function search_alpha_entropic_involution(lat::Lattice, f::Vector{Float64}, 
                                          feq::Vector{Float64};
                                          sbounds::Tuple{Real,Real}=__SBOUNDS)
  const F  = (alpha) -> (entropy_lat_boltzmann(lat, f + alpha * (feq - f) ) 
                        - entropy_lat_boltzmann(lat, f));
  const rs = Roots.fzeros(F, sbounds[1], sbounds[2]);
  return (length(rs) > 0) ? rs[end] : 0.0;
end

#! Bind search range to entropic involution search function
#TODO use FastAnnonymous module here
function init_search_alpha_entropic_involution(sbounds::Tuple{Real,Real})
  return (lat, f, feq) -> search_alpha_entropic_involution(lat, f, feq, 
                                                            sbounds=sbounds);
end

#! Search for alpha, the limit of over-relaxation for entropic involution
#!
#! \param   lat       Lattice
#! \param   f         Particle distribution vector
#! \param   feq       Equilibrium particle distribution vector
#! \param   x_0       Initial guess
#! \param   order     Order of method
#! \return            Limit of over-relaxation for entropic involution
function search_alpha_stef_entropic_involution(lat::Lattice, f::Vector{Float64}, 
                                               feq::Vector{Float64};
                                               x_0::Real=__X_0, 
                                               order::Int=__ORD)
  const F  = (alpha) -> (entropy_lat_boltzmann(lat, f + alpha * (feq - f) ) 
                        - entropy_lat_boltzmann(lat, f));
  return Roots.fzero(F, x_0, order=order);
end

#! Bind search range to entropic involution search function
#TODO use FastAnnonymous module here
function init_search_alpha_stef_entropic_involution(x_0::Real; order::Int=__ORD)
  return (lat, f, feq) -> search_alpha_entropic_involution(lat, f, feq, x_0=x_0,
                                                           order=order);
end

#! Search for alpha, the limit of over-relaxation for entropic involution
#!
#! \param   lat       Lattice
#! \param   f         Particle distribution vector
#! \param   feq       Equilibrium particle distribution vector
#! \return            Limit of over-relaxation for entropic involution
function search_alpha_entropic_involution_db(lat::Lattice, f::Vector{Float64}, 
                                             feq::Vector{Float64})
  const F  = (alpha) -> (entropy_lat_boltzmann(lat, f + alpha * (feq - f) ) 
                        - entropy_lat_boltzmann(lat, f));

  # Determine appropriate search upperbound
  bs = Array{Float64}(lat.n);
  for i=1:lat.n
    bs[i] = abs(f[i] / (f[i] - feq[i]));
  end
  
  const rs = Roots.fzeros(F, 0.0, minimum(bs));
  return (length(rs) > 0) ? rs[end] : 0.0;
end

#! Search for alpha, the limit of over-relaxation for entropic contraction
#!
#! \param   lat       Lattice
#! \param   f         Particle distribution vector
#! \param   feq       Equilibrium particle distribution vector
#! \param   sbounds   Search bounds, default (0.0, 2.0)
#! \param   omega     Entropic contraction factor
#! \return            Limit of over-relaxation for entropic involution
function search_alpha_entropic_contraction(lat::Lattice, f::Vector{Float64}, 
                                           feq::Vector{Float64};
                                           sbounds::Tuple{Real,Real}=__SBOUNDS,
                                           omega::Real=__OMEGA)
  @assert(1.0 <= omega && omega < 2.0, ("Entropic contraction factor must be " *
                                        "greater than or equal to 1.0 and "    *
                                        "less than two"));

  const F = (alpha) -> begin;
    const sf      = entropy_lat_boltzmann(lat, f);
    const sfeq    = entropy_lat_boltzmann(lat, feq);
    const sfp     = entropy_lat_boltzmann(lat, f + alpha * (feq - f));
    return sqrt(omega - 1) * (sf - sfeq) - (sfp - sfeq);
  end

  const rs        = Roots.fzeros(F, sbounds[1], sbounds[2]);
  return (length(rs) > 0) ? rs[end] : 0.0;
end

#! Bind search range and contraction factor to entropic contraction search
#TODO use FastAnnonymous module here
function init_search_alpha_entropic_contraction(sbounds::Tuple{Real,Real},
                                                omega::Real=__OMEGA)
  @assert(1.0 <= omega && omega < 2.0, ("Entropic contraction factor must be " *
                                        "greater than or equal to 1.0 and "    *
                                        "less than two"));
  return (lat, f, feq) -> search_alpha_entropic_contraction(lat, f, feq, 
                                                            sbounds=sbounds,
                                                            omega=omega);
end

#! Search for alpha, the limit of over-relaxation for entropic contraction
#!
#! \param   lat       Lattice
#! \param   f         Particle distribution vector
#! \param   feq       Equilibrium particle distribution vector
#! \param   omega     Entropic contraction factor
#! \return            Limit of over-relaxation for entropic involution
function search_alpha_entropic_contraction_db(lat::Lattice, f::Vector{Float64}, 
                                              feq::Vector{Float64};
                                              omega::Real=__OMEGA)
  @assert(1.0 <= omega && omega < 2.0, ("Entropic contraction factor must be " *
                                        "greater than or equal to 1.0 and "    *
                                        "less than two"));

  const F = (alpha) -> begin;
    const sf      = entropy_lat_boltzmann(lat, f);
    const sfeq    = entropy_lat_boltzmann(lat, feq);
    const sfp     = entropy_lat_boltzmann(lat, f + alpha * (feq - f));
    return sqrt(omega - 1) * (sf - sfeq) - (sfp - sfeq);
  end

  # Determine appropriate search upperbound
  bs = Array{Float64}(lat.n);
  for i=1:lat.n
    bs[i] = abs(f[i] / (f[i] - feq[i]));
  end

  const rs        = Roots.fzeros(F, 0.0, minimum(bs));
  return (length(rs) > 0) ? rs[end] : 0.0;
end

#! Helper functions to be used in Newton and other gradient based methods
const __F_ENTROPY = (lat, x_n, f, feq) -> begin;
  return (entropy_lat_boltzmann(lat, f + x_n * (feq - f) ) 
         - entropy_lat_boltzmann(lat, f));
end;

function __F_PRIME_ENTROPY(lat, x_n, f, feq)
  sum = 0.0;
  for i=1:lat.n
    fneq = feq[i] - f[i];
    sum -= (fneq * ( log( (f[i] + x_n * fneq) / lat.w[i] ) + 1 ));
  end
  return sum;
end

const __MAX_ITERS_NEWTON = convert(Int, 1e6);

#! Search for alpha using Newton's method
#! Details for this algorithm were mostly borrowed from Gorban and Packwood 2014
#!
#! \param   lat         Lattice
#! \param   f           Particle distribution vector
#! \param   feq         Equilibrium particle distribution vector
#! \param   x_0         Initial guess
#! \param   esp_search  Search tolerance
#! \return              Limit of over-relaxation for entropic involution
function search_alpha_newton_entropic_involution(lat::Lattice,
                                                 f::Vector{Float64}, 
                                                 feq::Vector{Float64};
                                                 x_0::Real=__X_0,
                                                 eps_search::Real=__EPS_SEARCH)
  @assert(eps_search >= 0.0, "Search tolerance must be nonnegative");

  const fneq_norm = norm(feq - f, 2);
  k = 0;
  x_n = x_0;

  while true

    if x_n == 0.0 # heuristic so that we don't get stuck at zero
      x_n += rand(0.0:1e-4:2.0);
    end

    fe  =  __F_ENTROPY(lat, x_n, f, feq);
    fep =  __F_PRIME_ENTROPY(lat, x_n, f, feq);

    x_n -= if feq != 0.0
             fe / fep;
           else
             fe;
           end

    #=for i=1:lat.n # heuristic, don't want to fall out of polytope
      if f[i] + x_n * (feq[i]-f[i]) < -2*eps()
        x_n = rand(0.0:1e-4:2.0);
        continue;
      end
    end=#

    #TODO perhaps wrap this in a try catch clause and just bail if shit goes bad
    #if fneq_norm * abs(__F_ENTROPY(lat, x_n, f, feq)) < eps_search 
    if abs(__F_ENTROPY(lat, x_n, f, feq)) < eps_search
      break;
    end

    k += 1;
    @assert(k <= __MAX_ITERS_NEWTON, "Root not found");
  end

  fe = __F_ENTROPY(lat, x_n, f, feq);
  return (fe >= 0) ? x_n : x_n - 2 * fe / __F_PRIME_ENTROPY(lat, x_n, f, feq); 
end

#! Bind initial guess to entropic involution search using Newton's method
#TODO use FastAnnonymous module here
function init_search_alpha_newton_entropic_involution(x_0::Real, 
                                                      eps_search::Real)
  return (lat, f, feq) -> (search_alpha_newton_entropic_involution(lat, f, feq, 
                            x_0=x_0, eps_search=eps_search));
end

#! Helper functions to be used in Newton and other gradient based methods
function __FC_ENTROPY(lat, x_n, f, feq, omega)
  const sfeq = entropy_lat_boltzmann(lat, feq);
  return (entropy_lat_boltzmann(lat, f + x_n * (feq - f) ) 
         - sfeq 
         - sqrt(omega-1) * (entropy_lat_boltzmann(lat, f) - sfeq));
end

#! Search for alpha, the limit of over-relaxation for entropic contraction
#!
#! \param   lat         Lattice
#! \param   f           Particle distribution vector
#! \param   feq         Equilibrium particle distribution vector
#! \param   x_0         Initial guess
#! \param   esp_search  Search tolerance
#! \param   omega       Entropic contraction factor
#! \return              Limit of over-relaxation for entropic contraction
function search_alpha_newton_entropic_contraction(lat::Lattice, 
                                                  f::Vector{Float64}, 
                                                  feq::Vector{Float64};
                                                  x_0::Real=__X_0,
                                                  eps_search::Real=__EPS_SEARCH,
                                                  omega::Real=__OMEGA)
  @assert(1.0 <= omega && omega < 2.0, ("Entropic contraction factor must be " *
                                        "greater than or equal to 1.0 and "    *
                                        "less than two"));
  @assert(eps_search >= 0.0, "Search tolerance must be nonnegative");

  const fneq_norm = norm(feq - f, 2);
  k = 0;
  x_n = x_0;

  while true

    if x_n == 0.0 # heuristic so that we don't get stuck at zero
      x_n += rand(0.0:1e-4:2.0);
    end

    fe  =  __FC_ENTROPY(lat, x_n, f, feq, omega);
    fep =  __F_PRIME_ENTROPY(lat, x_n, f, feq);

    x_n -= if fep != 0.0
             fe / fep;
           else
             fe;
           end

    #=for i=1:lat.n # heuristic, don't want to fall out of polytope
      if f[i] + x_n * (feq[i]-f[i]) < -2*eps()
        x_n = rand(0.0:1e-4:2.0);
        continue;
      end
    end=#

    #TODO perhaps wrap this in a try catch clause and just bail if shit goes bad
    #if fneq_norm * abs(__FC_ENTROPY(lat, x_n, f, feq, omega)) < eps_search 
    if abs(__F_ENTROPY(lat, x_n, f, feq)) < eps_search
      break;
    end

    k += 1;
    @assert(k <= __MAX_ITERS_NEWTON, "Root not found");
  end

  fe = __FC_ENTROPY(lat, x_n, f, feq, omega);
  return (fe >= 0) ? x_n : x_n - 2 * fe / __F_PRIME_ENTROPY(lat, x_n, f, feq);
end

#! Bind search range and contraction factor to entropic contraction search
#TODO use FastAnnonymous module here
function init_search_alpha_newton_entropic_contraction(x_0::Real,
                                                       eps_search::Real=__EPS_SEARCH,
                                                       omega::Real=__OMEGA)

  return (lat, f, feq) -> search_alpha_entropic_contraction(lat, f, feq,
                                                          x_0=x_0,
                                                          esp_search=eps_search,
                                                          omega=omega);
end

#! Test functions
using Base.Test;
import PyPlot;
function __test_search_alpha_newton_entropic_involution(x_0::Real=__X_0, 
                                                        eps_search::Real=__EPS_SEARCH)
  global time_spent_newton  = 0;
  global time_spent_root    = 0;
  const _sanei = init_search_alpha_newton_entropic_involution(x_0, eps_search);
  return (lat, f, feq) -> begin;
    const norm_fneq     = norm(f-feq, 2);

    t                   = time();
    const alpha_newton  = _sanei(lat, f, feq);
    time_spent_newton  += time() - t;
    
    t                   = time();
    const alpha_roots   = search_alpha_stef_entropic_involution(lat, f, feq);
    time_spent_root    += time() - t;
   
    const fe_an         = __F_ENTROPY(lat, alpha_newton, f, feq); 
    const fe_ar         = __F_ENTROPY(lat, alpha_roots, f, feq); 
    @show alpha_newton, alpha_roots;
    if abs(alpha_newton - alpha_roots) > 0.05 && false
      const xs = linspace(-0.2, max(alpha_newton, alpha_roots)+0.2, 1000);
      ys = Array{Float64}(1000);
      for i=1:1000; ys[i] = __F_ENTROPY(lat, xs[i], f, feq); end
      #PyPlot.scatter([alpha_newton; alpha_roots], [fe_an; fe_ar]);
      PyPlot.plot(xs, ys, xs, zeros(1000));
      PyPlot.show();
    end
    @test_approx_eq_eps fe_an 0.0 1e-10; 
    @test_approx_eq_eps fe_ar 0.0 1e-10; 
    @test_approx_eq_eps fe_an fe_ar 1e-9;
    @show time_spent_newton, time_spent_root, time_spent_newton-time_spent_root;
    return alpha_newton;
  end;
end

function __test_search_alpha_newton_entropic_contraction(x_0::Real=__X_0, 
                                                         eps_search::Real=__EPS_SEARCH,
                                                         omega::Real=__OMEGA)
  const _sanec = init_search_alpha_newton_entropic_contraction(x_0, eps_search,
                                                               omega);
  return (lat, f, feq) -> begin;
    const norm_fneq     = norm(f-feq, 2);
    const alpha_newton  = _sanec(lat, f, feq);
    const alpha_roots   = search_alpha_entropic_contraction_db(lat, f, feq,
                                                               omega=omega);
    @show alpha_newton, alpha_roots;
    if alpha_roots != 0.0
      @test_approx_eq_eps alpha_newton alpha_roots 0.05;
    end
    return alpha_newton;
  end;
end

# Default arguments
const __ALPHA_EPS         = 1e-5;   #TODO do something with this heuristic
const __KAPPA             = 0.5;    #TODO do something with this heuristic
const __EPS_DS            = 1e-15;  #TODO do something with this heuristic
const __SES               = search_alpha_newton_entropic_involution;
const __S                 = S_fallah;

#   Search constants
const __X_0               =     1.6;
const __SBOUNDS           =     (0.0, 2.5);
const __EPS_SEARCH        =     1e-10;
const __ORD               =     1;
const __OMEGA             =     1.5;
