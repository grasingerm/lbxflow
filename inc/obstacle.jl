# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Add a bounceback obstacle to the domain
#!
#! \param   active_cells    Matrix of active flags
#! \param   bcs!            Vector of boundary conditions
#! \param   i_min           ith index of beginning of obstacle
#! \param   i_max           ith index of end of obstacle
#! \param   j_min           jth index of beginning of obstacle
#! \param   j_max           jth index of end of obstacle
function _add_obstacle_bounce_back!(active_cells::Matrix{Bool}, 
                                    bcs!::Vector{LBXFunction},
                                    i_min::Int, i_max::Int, j_min::Int, 
                                    j_max::Int)
  active_cells[i_min+1:i_max-1, j_min+1:j_max-1] = false;
  push!(bcs!, @anon sim -> north_bounce_back!(sim.lat, i_min, i_max, j_min));
  push!(bcs!, @anon sim -> south_bounce_back!(sim.lat, i_min, i_max, j_max));
  push!(bcs!, @anon sim -> east_bounce_back!(sim.lat, i_min, j_min, j_max));
  push!(bcs!, @anon sim -> west_bounce_back!(sim.lat, i_max, j_min, j_max));
end

#! Add a reflecting obstacle to the domain
#!
#! \param   active_cells    Matrix of active flags
#! \param   bcs!            Vector of boundary conditions
#! \param   i_min           ith index of beginning of obstacle
#! \param   i_max           ith index of end of obstacle
#! \param   j_min           jth index of beginning of obstacle
#! \param   j_max           jth index of end of obstacle
function _add_obstacle_reflector!(active_cells::Matrix{Bool}, 
                                  bcs!::Vector{LBXFunction},
                                  i_min::Int, i_max::Int, j_min::Int, 
                                  j_max::Int)
  active_cells[i_min+1:i_max-1, j_min+1:j_max-1] = false;
  push!(bcs!, @anon sim -> north_reflect!(sim.lat, i_min, i_max, j_min));
  push!(bcs!, @anon sim -> south_reflect!(sim.lat, i_min, i_max, j_max));
  push!(bcs!, @anon sim -> east_reflect!(sim.lat, i_min, j_min, j_max));
  push!(bcs!, @anon sim -> west_reflect!(sim.lat, i_max, j_min, j_max));
end

const _obstacle_hash_table! = Dict([
  (:bounce_back,   _add_obstacle_bounce_back!),
  (:reflector,     _add_obstacle_reflector!) 
  ]);

#! Add an obstacle to the domain
#!
#! \param   active_cells    Matrix of active flags
#! \param   bcs!            Vector of boundary conditions
#! \param   i_min           ith index of beginning of obstacle
#! \param   i_max           ith index of end of obstacle
#! \param   j_min           jth index of beginning of obstacle
#! \param   j_max           jth index of end of obstacle
function add_obstacle!(active_cells::Matrix{Bool}, bcs!::Vector{LBXFunction},
                       i_min::Int, i_max::Int, j_min::Int, j_max::Int,
                       obs_type::Symbol = :bounce_back)
  if haskey(_obstacle_hash_table!, obs_type)
    _obstacle_hash_table![obs_type](active_cells, bcs!, i_min, i_max, j_min, 
                                    j_max);
  else
    error("Obstacle type, $obs_type, is not understood");
  end
end
