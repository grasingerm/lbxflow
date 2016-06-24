# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Prescribe velocity profile at inlet
function north_vprof!(sim::AbstractSim, i_range::UnitRange{Int}, j::Int, 
                      args...; 
                      analytical_soln::LBXFunction=analytical_poise_newton)
  north_velocity!(sim.lat, analytical_soln(args..., length(i_range)), 
                  i_range, j); 
end

#! Prescribe velocity profile at inlet
function south_vprof!(sim::AbstractSim, i_range::UnitRange{Int}, j::Int, 
                      args...; 
                      analytical_soln::LBXFunction=analytical_poise_newton)
  south_velocity!(sim.lat, analytical_soln(args..., length(i_range)), 
                  i_range, j); 
end

#! Prescribe velocity profile at inlet
function east_vprof!(sim::AbstractSim, i::Int, j_range::UnitRange{Int}, 
                      args...; 
                      analytical_soln::LBXFunction=analytical_poise_newton)
  east_velocity!(sim.lat, analytical_soln(args..., length(j_range)), 
                  i, j_range); 
end

#! Prescribe velocity profile at inlet
function west_vprof!(sim::AbstractSim, i::Int, j_range::UnitRange{Int}, 
                      args...; 
                      analytical_soln::LBXFunction=analytical_poise_newton)
  west_velocity!(sim.lat, analytical_soln(args..., length(j_range)), 
                  i, j_range); 
end
