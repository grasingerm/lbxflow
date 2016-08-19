# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Calculate the p norm distance from equilibrium
macro dist_from_eq_norm(f_neq::Vector{Float64}, p::Int)
  return :(norm(f_neq, p));
end

#! Stabilize mass
function stabilize_mass_callback(sim::FreeSurfSim, k::Real)
  const ni, nj = size(sim.msm.rho);
  for j=1:nj, i=1:ni
    if isnan(sim.tracker.M[i, j])
      mass_sum = 0;
      counter::UInt = 0;
      for k=sim.lat.n-1
        const i_nbr = i + sim.lat.c[1, k];
        const j_nbr = j + sim.lat.c[2, k];
        if !isnan(sim.tracker.M[i_nbr, j_nbr])
          counter += 1;
          mass_sum += sim.tracker.M[i_nbr, j_nbr];
        end
      end
      sim.tracker.M[i, j] = if counter != 0
                              mass_sum / counter;
                            else
                              0.0
                            end;
    end
  end
end
