# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#=require(abspath(joinpath(__simulate_root__, "..", "col", "collision.jl")));
require(abspath(joinpath(__simulate_root__, "..", "lattice.jl")));
require(abspath(joinpath(__simulate_root__, "..", "multiscale.jl")));
require(abspath(joinpath(__simulate_root__, "simtypes.jl")));
require(abspath(joinpath(__simulate_root__, "tracking.jl")));=#

#! stream particle densities
function stream!(lat::Lattice, temp_f::Array{Float64,3}, bounds::Array{Int64,2})
  const nbounds = size(bounds, 2);
  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max, k = 1:lat.n
      i_new = i + lat.c[1,k];
      j_new = j + lat.c[2,k];

      if i_new > i_max || j_new > j_max || i_new < i_min || j_new < j_min
        continue;
      end

      temp_f[k,i_new,j_new] = lat.f[k,i,j];
    end
  end

  copy!(lat.f, temp_f);
end

#! Stream particle densities in free surface conditions
#!
#! \param lat Lattice to stream on
#! \param temp_f Temp lattice to store streamed particle distributions on
#! \param bounds Boundaries enclosing active streaming regions
#! \param t Mass tracker
function stream!(lat::Lattice, temp_f::Array{Float64,3}, bounds::Array{Int64,2},
                 t::Tracker)
  const nbounds = size(bounds, 2);
  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max
      if t.state[i,j] == GAS; continue; end
      for k = 1:lat.n
        i_new = i + lat.c[1,k];
        j_new = j + lat.c[2,k];

        if (i_new > i_max || j_new > j_max || i_new < i_min || j_new < j_min 
            || t.state[i_new,j_new] == INTERFACE
            || t.state[i_new,j_new] == GAS)
          continue;
        end
        temp_f[k,i_new,j_new] = lat.f[k,i,j];
      end
    end
  end

  copy!(lat.f, temp_f);
end

#! Simulate a single step
function sim_step!(sim::Sim, temp_f::Array{Float64,3},
                   sbounds::Array{Int64,2}, collision_f!::Function, 
                   cbounds::Array{Int64,2}, bcs!::Array{Function})
  lat = sim.lat;
  msm = sim.msm;

  collision_f!(sim, cbounds);
  stream!(lat, temp_f, sbounds);

  for bc! in bcs!
    bc!(lat);
  end

  map_to_macro!(lat, msm);
end

#! Simulate a single step free surface flow step
function sim_step!(sim::FreeSurfSim, temp_f::Array{Float64,3},
                   sbounds::Array{Int64,2}, collision_f!::Function, 
                   cbounds::Array{Int64,2}, bcs!::Array{Function})
  lat = sim.lat;
  msm = sim.msm;
  t = sim.tracker;

  # Reconstruct missing distribution functions at the interface
  for node in t.interfacels
    f_reconst!(sim, t, node.val, sbounds, sim.rhog);
  end
  stream!(lat, temp_f, sbounds, t);
  masstransfer!(sim, sbounds); # Calculate mass transfer across interface
  collision_f!(sim, cbounds);
  update!(sim, sbounds); # Update the state of cells

  for bc! in bcs!
    bc!(lat);
  end

  map_to_macro!(lat, msm);
end

#! Run simulation
function simulate!(sim::AbstractSim, sbounds::Array{Int64,2},
                   collision_f!::Function, cbounds::Array{Int64,2},
                   bcs!::Array{Function}, n_steps::Int, test_for_term::Function,
                   callbacks!::Array{Function}, k::Int = 0)

  temp_f = copy(sim.lat.f);

  sim_step!(sim, temp_f, sbounds, collision_f!, cbounds, bcs!);

  for c! in callbacks!
    c!(sim, k+1);
  end

  prev_msm = MultiscaleMap(sim.msm);

  for i = k+2:n_steps
    try

      sim_step!(sim, temp_f, sbounds, collision_f!, cbounds, bcs!);

      for c! in callbacks!
        c!(sim, i);
      end

      # if returns true, terminate simulation
      if test_for_term(sim.msm, prev_msm)
        return i;
      end

      copy!(prev_msm.omega, sim.msm.omega);
      copy!(prev_msm.rho, sim.msm.rho);
      copy!(prev_msm.u, sim.msm.u);
    
    catch e

      showerror(STDERR, e);
      println();
      Base.show_backtrace(STDERR, catch_backtrace()); # display callstack
      warn("Simulation interrupted at step $i !");
      return i;

    end
  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::AbstractSim, sbounds::Array{Int64,2},
                   collision_f!::Function, cbounds::Array{Int64,2},
                   bcs!::Array{Function}, n_steps::Int,
                   callbacks!::Array{Function}, k::Int = 0)

  temp_f = copy(sim.lat.f);
  for i = k+1:n_steps
    try

      sim_step!(sim, temp_f, sbounds, collision_f!, cbounds, bcs!);

      for c! in callbacks!
        c!(sim, i);
      end

    catch e

      showerror(STDERR, e);
      println();
      Base.show_backtrace(STDERR, catch_backtrace()); # display callstack
      warn("Simulation interrupted at step $i !");
      return i;

    end
  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::AbstractSim, sbounds::Array{Int64,2},
                   collision_f!::Function, cbounds::Array{Int64,2},
                   bcs!::Array{Function}, n_steps::Int, k::Int = 0)

  temp_f = copy(sim.lat.f);
  for i = k+1:n_step
    try; sim_step!(sim, temp_f, sbounds, collision_f!, cbounds, bcs!);
    catch e
      showerror(STDERR, e); println();
      warn("Simulation interrupted at step $i !");
      return i;
    end
  end

  return n_steps;

end

