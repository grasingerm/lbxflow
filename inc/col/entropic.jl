# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

module Entropic

include("constitutive.jl");
include("forcing.jl");
include("equilibrium.jl");
include(joinpath("..", "lattice.jl"));
include("mrt_matrices.jl");
include(joinpath("..", "multiscale.jl"));
include(joinpath("..", "numerics.jl"));
include(joinpath("..", "sim", "simtypes.jl"));

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
function init_col_srt(constit_relation_f::Function; 
                      feq_f::Function=feq_incomp_max_entropy,
                      kappa=0.5::Real, eps_ds=1e-15::Real,
                      search_entropic_stabiliy=search_alpha_entropic_involution)
  @assert(0.0 <= kappa && kappa <= 1.0, __KAPPA_ERROR_MSG__);

  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        rhoij = msm.rho[i,j];
        uij = msm.u[:,i,j];
        feq = Array(Float64, lat.n); 
        fneq = Array(Float64, lat.n); 
        for k = 1:lat.n 
          feq[k] = feq_f(lat, rhoij, uij, k);
          fneq[k] = lat.f[k,i,j] - feq[k];
        end
        const mu = constit_relation_f(sim, fneq, i, j);
        const omega = @omega(mu, lat.cssq, lat.dt);
        if (entropy_lat_boltzmann(lat, f_eq) - entropy_lat_boltzmann(lat, f)
            < eps_ds) # do not to enforce entropic stability, use standard BGK
          for k = 1:lat.n
            lat.f[k,i,j] += omega * (feq[k] - lat.f[k,i,j]);
          end
        else
          const alpha = search_alpha_entropic_involution(lat, f, f_eq);
          if alpha < 1e-4 #TODO fix this heuristic...
            for k = 1:lat.n
              warn("Root for entropic stabilization not found at node ($i, $j)");
              warn("Solution may become unstable due to decrease in entropy.");
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
function init_col_srt(constit_relation_f::Function,
                      forcing_kf::Tuple{Function, Function};
                      feq_f::Function=feq_incomp_max_entropy,
                      kappa=0.5::Real, eps_ds=1e-15::Real,
                      search_entropic_stabiliy=search_alpha_entropic_involution)
  @assert(0.0 <= kappa && kappa <= 1.0, __KAPPA_ERROR_MSG__);

  const uf, colf = forcing_kf;
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        rhoij = msm.rho[i,j];
        uij = uf(lat, msm.u[:,i,j]);
        feq = Array(Float64, lat.n); 
        fneq = Array(Float64, lat.n); 
        for k = 1:lat.n 
          feq[k] = feq_f(lat, rhoij, uij, k);
          fneq[k] = lat.f[k,i,j] - feq[k];
        end
        const mu = constit_relation_f(sim, fneq, i, j);
        const omega = @omega(mu, lat.cssq, lat.dt);
        if (entropy_lat_boltzmann(lat, f_eq) - entropy_lat_boltzmann(lat, f)
            < eps_ds) # do not to enforce entropic stability, use standard BGK
          for k = 1:lat.n
            lat.f[k,i,j] += (omega * (feq[k] - lat.f[k,i,j]) 
                             + colf(lat, omega, uij, k));
          end
        else
          const alpha = search_alpha_entropic_involution(lat, f, f_eq);
          if alpha < 1e-4 #TODO fix this heuristic...
            for k = 1:lat.n
              warn("Root for entropic stabilization not found at node ($i, $j)");
              warn("Solution may become unstable due to decrease in entropy.");
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
  end
end


#! Entropically stabilized MRT collision operator for incompressible flow
#!
#! \param constit_relation_f  Constitutive relationship
#! \param feq_f               Equilibrium particle distribution function
#! \param kappa               Entropic relaxation factor
#! \param eps_ds              Change in entropy threshold for ELBM
#! \param search_entropic_stabiliy Search function for stability factor, alpha
#! \return collision_function!(sim, bounds)
function init_col_mrt(constit_relation_f::Function;
                      feq_f::Function=feq_incomp_max_entropy,
                      kappa=0.5::Real, eps_ds=1e-15::Real,
                      search_entropic_stabiliy=search_alpha_entropic_involution)
  @assert(0.0 <= kappa && kappa <= 1.0, __KAPPA_ERROR_MSG__);
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const M = @DEFAULT_MRT_M();
    const iM = @DEFAULT_MRT_IM();
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        rhoij     = msm.rho[i,j];
        uij       = msm.u[:,i,j];
        feq       = Array(Float64, lat.n); 
        fneq      = Array(Float64, lat.n); 
        for k = 1:lat.n 
          feq[k]  = feq_f(lat, rhoij, uij, k);
          fneq[k] = lat.f[k,i,j] - feq[k];
        end
        const fij  = lat.f[:,i,j];
        const muij = constit_relation_f(sim, fneq, S, M, iM, i, j);
        const Sij  = S(muij, rhoij, lat.cssq, lat.dt);
        if (entropy_lat_boltzmann(lat, f_eq) - entropy_lat_boltzmann(lat, f)
            < eps_ds) # do not to enforce entropic stability, use standard BGK
          lat.f[:,i,j] = fij - iM * Sij * M * fneq; # perform collision
        else
          const alpha = search_alpha_entropic_involution(lat, f, f_eq);
          if alpha < 1e-4 #TODO fix this heuristic...
            warn("Root for entropic stabilization not found at node ($i, $j)");
            warn("Solution may become unstable due to decrease in entropy.");
            lat.f[:,i,j] = fij - iM * Sij * M * fneq;
          else
            for a in (8,9); Sij[a] *= kappa * alpha; end;
            lat.f[:,i,j] = fij - iM * Sij * M * fneq;
          end
        end

        msm.omega[i,j] = @omega(mu, lat.cssq, lat.dt);;
      end
    end
  end
end

#! Entropically stabilized MRT collision operator for incompressible flow
#!
#! \param constit_relation_f  Constitutive relationship
#! \param forcing_kf          Forcing functions
#! \param feq_f               Equilibrium particle distribution function
#! \param kappa               Entropic relaxation factor
#! \param eps_ds              Change in entropy threshold for ELBM
#! \param search_entropic_stabiliy Search function for stability factor, alpha
#! \return collision_function!(sim, bounds)
function init_col_mrt(constit_relation_f::Function,
                      forcing_kf::Tuple{Function, Function};
                      feq_f::Function=feq_incomp_max_entropy,
                      kappa=0.5::Real, eps_ds=1e-15::Real,
                      search_entropic_stabiliy=search_alpha_entropic_involution)
  @assert(0.0 <= kappa && kappa <= 1.0, __KAPPA_ERROR_MSG__);
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const M = @DEFAULT_MRT_M();
    const iM = @DEFAULT_MRT_IM();
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        rhoij     = msm.rho[i,j];
        uij       = uf(lat, msm.u[:,i,j]);
        feq       = Array(Float64, lat.n); 
        fneq      = Array(Float64, lat.n); 
        for k = 1:lat.n 
          feq[k]  = feq_f(lat, rhoij, uij, k);
          fneq[k] = lat.f[k,i,j] - feq[k];
        end

        const fij     = lat.f[:,i,j];
        const muij    = constit_relation_f(sim, fneq, S, M, iM, i, j);
        const Sij     = S(muij, rhoij, lat.cssq, lat.dt);
        const omegaij = @omega(mu, lat.cssq, lat.dt);;
        fdl           = Array(Float64, lat.n);
        for k = 1:lat.n
          fdl[k] = colf(lat, omegaij, uij, k);
        end

        if (entropy_lat_boltzmann(lat, f_eq) - entropy_lat_boltzmann(lat, f)
            < eps_ds) # do not to enforce entropic stability, use standard BGK
          lat.f[:,i,j]   = fij - iM * Sij * M * fneq + fdl;
        else
          const alpha    = search_alpha_entropic_involution(lat, f, f_eq);
          if alpha < 1e-4 #TODO fix this heuristic...
            warn("Root for entropic stabilization not found at node ($i, $j)");
            warn("Solution may become unstable due to decrease in entropy.");
            lat.f[:,i,j] = fij - iM * Sij * M * fneq + fdl;
          else
            for a in (8,9); Sij[a] *= kappa * alpha; end;
            lat.f[:,i,j] = fij - iM * Sij * M * fneq + fdl;
          end
        end

        msm.omega[i,j] = omegaij;
      end
    end
  end
end

#! Search for alpha, the limit of over-relaxation for entropic involution
#!
#! \param   lat       Lattice
#! \param   f         Particle distribution vector
#! \param   f_eq      Equilibrium particle distribution vector
#! \param   sbounds   Search bounds, default (0.0, 2.0)
#! \return            Limit of over-relaxation for entropic involution
function search_alpha_entropic_involution(lat::Lattice, f::Vector{Float64}, 
                                          f_eq::Vector{Float64};
                                          sbounds=(0.0, 2.0)::Tuple{Real,Real})
  const F  = (alpha) -> (entropy_lat_boltzmann(lat, f + alpha * (f_eq - f) ) 
                        - entropy_lat_boltzmann(lat, f));
  const rs = Roots.fzeros(F, sbounds[1], sbounds[2]);
  return (length(rs) > 0) ? rs[end] : 0.0;
end

#! Bind search range to entropic involution search function
#TODO use FastAnnonymous module here
function init_search_alpha_entropic_involution(sbounds::Tuple{Real,Real})
  return (lat, f, f_eq) -> search_alpha_entropic_involution(lat, f, f_eq, 
                                                            sbounds=sbounds);
end

#! Search for alpha, the limit of over-relaxation for entropic involution
#!
#! \param   lat       Lattice
#! \param   f         Particle distribution vector
#! \param   f_eq      Equilibrium particle distribution vector
#! \param   sbounds   Search bounds, default (0.0, 2.0)
#! \param   omega     Entropic contraction factor
#! \return            Limit of over-relaxation for entropic involution
function search_alpha_entropic_contraction(lat::Lattice, f::Vector{Float64}, 
                                           f_eq::Vector{Float64};
                                           sbounds=(1.0, 2.0)::Tuple{Real,Real},
                                           omega=1.5::Real)
  @assert(1.0 <= omega && omega < 2.0, ("Entropic contraction factor must be " *
                                        "greater than or equal to 1.0 and "    *
                                        "less than two"));

  const F = (alpha) -> begin;
    const sf      = entropy_lat_boltzmann(lat, f);
    const sf_eq   = entropy_lat_boltzmann(lat, f_eq);
    const sfp     = entropy_lat_boltzmann(lat, f + alpha * (f_eq - f));
    return sqrt(omega - 1) * (sf - sf_eq) - (sfp - sf_eq);
  end

  const rs        = Roots.fzeros(F, sbounds[1], sbounds[2]);
  return (length(rs) > 0) ? rs[end] : 0.0;
end

#! Bind search range and contraction factor to entropic contraction search
#TODO use FastAnnonymous module here
function init_search_alpha_entropic_contraction(sbounds::Tuple{Real,Real},
                                                omega=1.5::Real)
  @assert(1.0 <= omega && omega < 2.0, ("Entropic contraction factor must be " *
                                        "greater than or equal to 1.0 and "    *
                                        "less than two"));
  return (lat, f, f_eq) -> search_alpha_entropic_contraction(lat, f, f_eq, 
                                                             sbounds=sbounds,
                                                             omega=omega);
end

end # module Entropic
