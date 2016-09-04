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

#! Mass inlet
function mass_inlet!(sim::FreeSurfSim, i::Int, j_range::UnitRange{Int})
  for j in j_range
    if sim.tracker.state[i, j] != GAS
      sim.tracker.M[i, j] = sim.msm.rho[i, j];
    end
  end
end

#! Mass inlet
function mass_inlet!(sim::FreeSurfSim, i_range::UnitRange{Int}, j::Int)
  for i in i_range
    if sim.tracker.state[i, j] != GAS
      sim.tracker.M[i, j] = sim.msm.rho[i, j];
    end
  end
end

#! Mass inlet
function north_mass_inlet!(sim::FreeSurfSim)
  const ni, nj = size(sim.msm.rho);
  mass_inlet!(sim, 1:ni, nj);
end

#! Mass inlet
function south_mass_inlet!(sim::FreeSurfSim)
  const ni = size(sim.msm.rho, 1);
  mass_inlet!(sim, 1:ni, 1);
end

#! Mass inlet
function east_mass_inlet!(sim::FreeSurfSim)
  const ni, nj = size(sim.msm.rho);
  mass_inlet!(sim, ni, 1:nj);
end

#! Mass inlet
function west_mass_inlet!(sim::FreeSurfSim)
  const ni, nj = size(sim.msm.rho);
  mass_inlet!(sim, 1, 1:nj);
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

#! Mass outlet
function mass_boutlet!(sim::FreeSurfSim, i::Int, j_range::UnitRange{Int}, 
                       ks::Vector{Int}, inner_bc!::LBXFunction)
  for j in j_range
    if sim.tracker.state[i, j] == FLUID
      inner_bc!(sim.lat, i, start(j_range), endof(j_range));
      for k in ks
        sim.tracker.M[i, j] -= sim.lat.f[k, i, j];
      end
    end
  end
end

#! Mass outlet
function mass_boutlet!(sim::FreeSurfSim, i_range::UnitRange{Int}, j::Int, 
                       ks::Vector{Int}, inner_bc!::LBXFunction)
  for i in i_range
    if sim.tracker.state[i, j] == FLUID
      inner_bc!(sim.lat, start(i_range), endof(i_range), j);
      for k in ks
        sim.tracker.M[i, j] -= sim.lat.f[k, i, j];
      end
    end
  end
end

#TODO use type system to get correct velocity vector indices for each lattice,
#     e.g. mass_outlet!(sim, 1:ni, nj, north_vecs(sim.lat));

#! Mass outlet
function north_mass_outlet!(sim::FreeSurfSim)
  const ni, nj = size(sim.msm.rho);
  mass_outlet!(sim, 1:ni, nj, [6; 2; 5]);
end

#! Mass outlet
function south_mass_outlet!(sim::FreeSurfSim)
  const ni = size(sim.msm.rho, 1);
  mass_outlet!(sim, 1:ni, 1, [7; 4; 8]);
end

#! Mass outlet
function east_mass_outlet!(sim::FreeSurfSim)
  const ni, nj = size(sim.msm.rho);
  mass_outlet!(sim, ni, 1:nj, [5; 1; 8]);
end

#! Mass outlet
function west_mass_outlet!(sim::FreeSurfSim)
  const nj = size(sim.msm.rho, 2);
  mass_outlet!(sim, 1, 1:nj, [6; 3; 7]);
end

#! Mass outlet
function north_mass_open!(sim::FreeSurfSim)
  const ni, nj = size(sim.msm.rho);
  mass_boutlet!(sim, 1:ni, nj, [6; 2; 5], north_open!);
end

#! Mass outlet
function south_mass_open!(sim::FreeSurfSim)
  const ni = size(sim.msm.rho, 1);
  mass_boutlet!(sim, 1:ni, 1, [7; 4; 8], south_open!);
end

#! Mass outlet
function east_mass_open!(sim::FreeSurfSim)
  const ni, nj = size(sim.msm.rho);
  mass_boutlet!(sim, ni, 1:nj, [5; 1; 8], east_open!);
end

#! Mass outelt
function west_mass_open!(sim::FreeSurfSim)
  const nj = size(sim.msm.rho, 2);
  mass_boutlet!(sim, 1, 1:nj, [6; 3; 7], west_open!);
end
