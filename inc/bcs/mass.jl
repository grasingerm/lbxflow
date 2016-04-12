# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Mass inlet
function mass_inlet!(sim::FreeSurfSim, i::Int, j_range::UnitRange{Int}, m::Real)
  for j in j_range
    sim.tracker.M[i, j] = m;
  end
end

#! Mass inlet
function mass_inlet!(sim::FreeSurfSim, i_range::UnitRange{Int}, j::Int, m::Real)
  for i in i_range
    sim.tracker.M[i, j] = m;
  end
end

#! Mass inlet
function mass_inlet_north!(sim::FreeSurfSim, m::Real)
  const ni, nj = size(sim.msm.rho);
  mass_inlet!(sim, 1:ni, nj, m);
end

#! Mass inlet
function mass_inlet_south!(sim::FreeSurfSim, m::Real)
  const ni = size(sim.msm.rho, 1);
  mass_inlet!(sim, 1:ni, 1, m);
end

#! Mass inlet
function mass_inlet_east!(sim::FreeSurfSim, m::Real)
  const ni, nj = size(sim.msm.rho);
  mass_inlet!(sim, ni, 1:nj, m);
end

#! Mass inlet
function mass_inlet_west!(sim::FreeSurfSim, m::Real)
  const ni, nj = size(sim.msm.rho);
  mass_inlet!(sim, 1, 1:nj, m);
end
