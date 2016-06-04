# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

function adapt_time_step!(sim::AdaptiveTimeStepSim, col_f!::ColFunction)
  lat           =   sim.isim.lat;
  msm           =   sim.isim.msm;
  const ni, nj  =   size(sim.isim.msm.rho);
  const nk      =   lat.n;

  u_mag_max = maximum(u_mag(sim.msm));

  # Decrease time step
  if      u_mag_max > (1/6 / sim.ξ) && maximum(isim.msm.omega) < 1.99
    const Δt_o  =   sim.Δt;
    Δt_n        =   Δt_o * sim.ξ;
    const st    =   Δt_n / Δt_o;
    const ρ_med =   sum(sim.isim.tracker.eps) / sum(sim.isim.tracker.M);
    ρ_n         =   map(@anon ρ -> st * (ρ - ρ_med) + ρ_med, msm.rho);
    for j=1:nj, i=1:ni
      sim.isim.tracker.M[i, j] *= msm.rho[i, j] / ρ_n[i, j];
    end
    copy!(msm.rho, ρ_n);

  # Decrease time step
  elseif  u_mag_max < (sim.ξ / 6)

  # Basically, do nothing
  else
  end
end

#! Rescale body forces
function rescale(col_f!::Union{BGK_F, MRT_F}, st::Real)
  try
    col_f!.forcing_f.c *= st*st;
  catch e
    warn("Forcing functions in adaptive simulations must scalable");
    rethrow(e)
  end
end

#! Rescale body forces
function rescale(col_f!::FltrColFunction, st::Real)
  # forward call to inner collision function
  rescale(col_f!.inner_col_f!, st);
end

#! Default should be to do nothing...
function rescale(col_f!::ColFunction, st::Real)
  return;
end
