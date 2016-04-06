# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

function adapt_time_step!(isim::FreeSurfSim, ξ::Real)
  const ni, nj  =   size(isim.tracker.state);
  const nk      =   lat.n;

  u_mag_max = 0.0;
  if (isim.tracker.states[i, j] != GAS && 
      (u_mag_ij = norm(isim.msm.u[:, i, j], 2)) > u_mag_max)
    u_mag_max = u_mag_ij;
  end

  # Decrease time step
  if      u_max_max > (1/6 / ξ) && maximum(isim.msm.omega) < 1.99

  # Decrease time step
  elseif  u_mag_max < (ξ / 6)

  # Basically, do nothing
  else
  end
end
