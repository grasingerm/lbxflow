# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

function adapt_time_step!(sim::AdaptiveTimeStepSim, col_f!::ColFunction)
  _adapt_time_step!(sim, sim.isim, col_f!::ColFunction);
end

function _adapt_time_step!(sim::AdaptiveTimeStepSim, isim::AbstractSim, 
                           col_f!::ColFunction)
  lat           =   sim.isim.lat;
  msm           =   sim.isim.msm;
  const ni, nj  =   size(sim.isim.msm.rho);
  const nk      =   lat.n;

  u_mag_max = maximum(u_mag(sim.msm));

  if  ((u_mag_max > (1/6 / sim.ξ) && sim.decr) 
       || (u_mag_max < (sim.ξ / 6) && maximum(isim.msm.omega) < 1.99 
           && sim.incr))
    const Δt_o  =   sim.Δt;
    
                    # Decrease time step
    const Δt_n  =   if      u_mag_max > (1/6 / sim.ξ) && sim.decr
                      val = Δt_o * sim.ξ * sim.relax + Δt_o * (1 - sim.relax);
                      info("Decreasing time step size $(Δt_o) => $(val)");
                      val; 
                    # Increase time step
                    elseif  (u_mag_max < (sim.ξ / 6) && sim.incr &&
                             maximum(isim.msm.omega) < 1.95)
                      val = Δt_o / sim.ξ * sim.relax + Δt_o * (1 - sim.relax);
                      info("Increasing time step size $(Δt_o) => $(val)");
                      val;
                    else
                      error("This shouldn't have happened, fam");
                    end;
    
    const st    =   Δt_n / Δt_o;
    ρ_n         =   map(ρ -> st * ρ, msm.rho);
    ω_n         =   map(ω -> 1.0 / ( st*(1.0 / ω - 0.5) + 0.5), msm.omega);
    u_n         =   map(u -> st * u, msm.u);
    for j=1:nj, i=1:ni # rescale mass, fluid fraction, and particle distributions
      s_ω           =   st * msm.omega[i, j] / ω_n[i, j];
      for k=1:nk
        const feq_o    =   col_f!.feq_f(lat, ρ_n[i, j], sub(u_n, :, i, j), k);
        const feq_n    =   col_f!.feq_f(lat, msm.rho[i, j], sub(msm.u, :, i, j), k);
        const s_f      =   feq_n / feq_o; 
        lat.f[k, i, j] =   s_f * (feq_o + s_ω * (lat.f[k, i, j] - feq_o));
      end
    end

    sim.Δt      =  Δt_n; 
    copy!(msm.rho, ρ_n);
    copy!(msm.u, u_n);
    copy!(msm.omega, ω_n);
    rescale!(col_f!, st);
  end
end

function _adapt_time_step!(sim::AdaptiveTimeStepSim, isim::FreeSurfSim, 
                           col_f!::ColFunction)
  lat           =   sim.isim.lat;
  msm           =   sim.isim.msm;
  const ni, nj  =   size(sim.isim.msm.rho);
  const nk      =   lat.n;

  u_mag_max = maximum(u_mag(sim.msm));

  if  ((u_mag_max > (1/6 / sim.ξ) && sim.decr) 
       || (u_mag_max < (sim.ξ / 6) && maximum(isim.msm.omega) < 1.95 
           && sim.incr))
    const Δt_o  =   sim.Δt;
    
                    # Decrease time step
    const Δt_n  =   if      u_mag_max > (1/6 / sim.ξ) && sim.decr
                      val = Δt_o * sim.ξ * sim.relax + Δt_o * (1 - sim.relax);
                      info("Decreasing time step size $(Δt_o) => $(val)");
                      val; 
                    # Increase time step
                    elseif  (u_mag_max < (sim.ξ / 6) && sim.incr &&
                             maximum(isim.msm.omega) < 1.95)
                      val = Δt_o / sim.ξ * sim.relax + Δt_o * (1 - sim.relax);
                      info("Increasing time step size $(Δt_o) => $(val)");
                      val;
                    else
                      error("This shouldn't have happened, fam");
                    end;
    
    const st    =   Δt_n / Δt_o;
    const ρ_med =   sum(sim.isim.tracker.eps) / sum(sim.isim.tracker.M);
    ρ_n         =   map(ρ -> st * (ρ - ρ_med) + ρ_med, msm.rho);
    ω_n         =   map(ω -> 1.0 / ( st*(1.0 / ω - 0.5) + 0.5), msm.omega);
    u_n         =   map(u -> st * u, msm.u);
    for j=1:nj, i=1:ni # rescale mass, fluid fraction, and particle distributions
      @inbounds sim.isim.tracker.M[i, j] *= msm.rho[i, j] / ρ_n[i, j];
      @inbounds sim.isim.tracker.eps[i, j] = sim.isim.tracker.M[i, j] / ρ_n[i, j];
      s_ω           =   st * msm.omega[i, j] / ω_n[i, j];
      for k=1:nk
        const feq_o    =   col_f!.feq_f(lat, ρ_n[i, j], sub(u_n, :, i, j), k);
        const feq_n    =   col_f!.feq_f(lat, msm.rho[i, j], sub(msm.u, :, i, j), k);
        const s_f      =   feq_n / feq_o; 
        lat.f[k, i, j] =   s_f * (feq_o + s_ω * (lat.f[k, i, j] - feq_o));
      end
    end

    sim.Δt      =  Δt_n; 
    copy!(msm.rho, ρ_n);
    copy!(msm.u, u_n);
    copy!(msm.omega, ω_n);
    rescale!(col_f!, st);
  end
end

#! Rescale body forces and if applicable, the constant \mu
function rescale!(col_f!::Union{BGK_F, MRT_F}, st::Real)
  try
    col_f!.forcing_f.c *= st*st;
  catch e
    warn("Forcing functions in adaptive simulations must scalable");
    rethrow(e)
  end
  _rescale_constit_relation_f!(col_f!.constit_relation_f, st);
end

# Rescale constant constitutive relationships, i.e. constitutive relationships
#   that do NOT depend on the local strain rate
function _rescale_constit_relation_f!(cc::_ConstConstit, st::Real)
  cc.μ *= st;  
end
# Default, do nothing
function _rescale_constit_relation_f!(cr::LBXFunction, st::Real)
  return; 
end

#! Rescale body forces
function rescale!(col_f!::FltrColFunction, st::Real)
  # forward call to inner collision function
  rescale!(col_f!.inner_col_f!, st);
end

#! Default should be to do nothing...
function rescale!(col_f!::ColFunction, st::Real)
  return;
end
