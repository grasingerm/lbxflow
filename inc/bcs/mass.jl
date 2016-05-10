# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Mass inlet
function mass_inlet!(sim::FreeSurfSim, i::Int, j_range::UnitRange{Int}, m::Real)
  for j in j_range
    if sim.tracker.state[i, j] != GAS; sim.tracker.M[i, j] = m; end
  end
end

#! Mass inlet
function mass_inlet!(sim::FreeSurfSim, i_range::UnitRange{Int}, j::Int, m::Real)
  for i in i_range
    if sim.tracker.state[i, j] != GAS; sim.tracker.M[i, j] = m; end
  end
end

#! Mass inlet
function north_mass_inlet!(sim::FreeSurfSim, m::Real)
  const ni, nj = size(sim.msm.rho);
  mass_inlet!(sim, 1:ni, nj, m);
end

#! Mass inlet
function south_mass_inlet!(sim::FreeSurfSim, m::Real)
  const ni = size(sim.msm.rho, 1);
  mass_inlet!(sim, 1:ni, 1, m);
end

#! Mass inlet
function east_mass_inlet!(sim::FreeSurfSim, m::Real)
  const ni, nj = size(sim.msm.rho);
  mass_inlet!(sim, ni, 1:nj, m);
end

#! Mass inlet
function west_mass_inlet!(sim::FreeSurfSim, m::Real)
  const ni, nj = size(sim.msm.rho);
  mass_inlet!(sim, 1, 1:nj, m);
end

#! Mass outlet
function mass_outlet!(sim::FreeSurfSim, i::Int, j_range::UnitRange{Int}, 
                      ks::Vector{Int})
  for j in j_range
    if sim.tracker.state[i, j] == FLUID
      for k in ks
        sim.tracker.M[i, j] -= sim.lat.f[k, i, j];
      end
    end
  end
end

#! Mass outlet
function mass_outlet!(sim::FreeSurfSim, i_range::UnitRange{Int}, j::Int, 
                      ks::Vector{Int})
  for i in i_range
    if sim.tracker.state[i, j] == FLUID
      for k in ks
        sim.tracker.M[i, j] -= sim.lat.f[k, i, j];
      end
    end
  end
end

#TODO use type system to get correct velocity vector indices for each lattice,
#     e.g. mass_outlet!(sim, 1:ni, nj, north_vecs(sim.lat));

#! Mass inlet
function north_mass_outlet!(sim::FreeSurfSim)
  const ni, nj = size(sim.msm.rho);
  mass_outlet!(sim, 1:ni, nj, [6; 2; 5]);
end

#! Mass inlet
function south_mass_outlet!(sim::FreeSurfSim)
  const ni = size(sim.msm.rho, 1);
  mass_outlet!(sim, 1:ni, 1, [7; 4; 8]);
end

#! Mass inlet
function east_mass_outlet!(sim::FreeSurfSim)
  const ni, nj = size(sim.msm.rho);
  mass_outlet!(sim, ni, 1:nj, [5; 1; 8]);
end

#! Mass inlet
function west_mass_outlet!(sim::FreeSurfSim)
  const nj = size(sim.msm.rho, 2);
  mass_outlet!(sim, 1, 1:nj, [6; 3; 7]);
end
