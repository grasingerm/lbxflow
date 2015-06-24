const __simulate_root__ = dirname(@__FILE__);
require(abspath(joinpath(__simulate_root__, "collision.jl")));
require(abspath(joinpath(__simulate_root__, "lattice.jl")));
require(abspath(joinpath(__simulate_root__, "multiscale.jl")));

abstract AbstractSim

#! Simulation object
#! lat Lattice
#! msm Multiscale map
immutable Sim <: AbstractSim
  lat::Lattice
  msm::MultiscaleMap

  Sim(lat::Lattice, msm::MultiscaleMap) = new(lat, msm);
end

#! Free surface flow simulation object
#! lat Lattice
#! msm Multiscal map
#! tracker Mass tracker
immutable FreeSurfSim <: AbstractSim
  lat::Lattice
  msm::MultiscaleMap
  tracker::T

  function FreeSurfSim(lat::Lattice, msm::MultiscaleMap, tracker::Tracker)
    return new(lat, msm, tracker);
  end
end

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

#! stream particle densities
function stream!(lat::Lattice, temp_f::Array{Float64,3}, bounds::Array{Int64,2},
                 msm::MultiscaleMap, mass_tracker::Tracker)
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

      # mass tracking
      neighbor_state = mass_tracker.state[i,j];
      if neighbor_state != GAS
        opk = opp_lat_vec_D2Q9(k);
        if neighbor_state == LIQUID
          dM = lat.f[opk,i_new,j_new] - lat.f[k,i,j];
        elseif neighbor_state == INTERFACE
          epsx = msm.rho[i,j] / mass_tracker.M[i,j];
          epsxe = msm.rho[i_new,j_new] / mass_tracker.M[i_new,j_new];
          dM = (epsx + epsxe)/2 * (lat.f[opk,i_new,j_new] - lat.f[k,i,j]);
        else
          error("State $neighbor_state is not understood!");
        end
        M[i,j] += dM;
      end

      temp_f[k,i_new,j_new] = lat.f[k,i,j];
    end
  end

  copy!(lat.f, temp_f);
end

#! Simulate a single step
function sim_step!(lat::Lattice, temp_f::Array{Float64,3}, msm::MultiscaleMap,
                   sbounds::Array{Int64,2}, collision_f!::Function, 
                   cbounds::Array{Int64,2}, bcs!::Array{Function})

  collision_f!(lat, msm, cbounds);
  stream!(lat, temp_f, sbounds);

  for bc! in bcs!
    bc!(lat);
  end

  map_to_macro!(lat, msm);

end

#! Run simulation
function simulate!(sim::Sim, sbounds::Array{Int64,2}, collision_f!::Function,
                   cbounds::Array{Int64,2}, bcs!::Array{Function}, n_steps::Int,
                   test_for_term::Function, callbacks!::Array{Function}, k = 0)

  temp_f = copy(sim.lat.f);

  sim_step!(sim.lat, temp_f, sim.msm, sbounds, collision_f!, cbounds, bcs!);

  for c! in callbacks!
    c!(sim, k+1);
  end

  prev_msm = MultiscaleMap(sim.msm);

  for i = k+2:n_steps
    try

      sim_step!(sim.lat, temp_f, sim.msm, sbounds, collision_f!, cbounds, bcs!);

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
function simulate!(sim::Sim, sbounds::Array{Int64,2}, collision_f!::Function,
                   cbounds::Array{Int64,2}, bcs!::Array{Function}, n_steps::Int,
                   callbacks!::Array{Function}, k::Int = 0)

  temp_f = copy(sim.lat.f);
  for i = k+1:n_steps
    try

      sim_step!(sim.lat, temp_f, sim.msm, sbounds, collision_f!, cbounds, bcs!);

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
function simulate!(sim::Sim, sbounds::Array{Int64,2}, collision_f!::Function,
  cbounds::Array{Int64,2}, bcs!::Array{Function}, n_steps::Int, k::Int = 0)

  temp_f = copy(sim.lat.f);
  for i = k+1:n_steps
    try; sim_step!(sim.lat, temp_f, sim.msm, sbounds, collision_f!, cbounds, bcs!);
    catch e
      showerror(STDERR, e); println(); warn("Simulation interrupted at step $i !");
      return i;
    end
  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::FreeSurfSim, sbounds::Array{Int64,2},
                   collision_f!::Function, cbounds::Array{Int64,2}, 
                   bcs!::Array{Function}, n_steps::Int, test_for_term::Function,
                   callbacks!::Array{Function}, k = 0)
  error("not yet implemented");
  temp_f = copy(sim.lat.f);

  sim_step!(sim.lat, temp_f, sim.msm, sbounds, collision_f!, cbounds, bcs!);

  for c! in callbacks!
    c!(sim, k+1);
  end

  prev_msm = MultiscaleMap(sim.msm);

  for i = k+2:n_steps
    try

      sim_step!(sim.lat, temp_f, sim.msm, sbounds, collision_f!, cbounds, bcs!);

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
function simulate!(sim::FreeSurfSim, sbounds::Array{Int64,2},
                   collision_f!::Function, cbounds::Array{Int64,2},
                   bcs!::Array{Function}, n_steps::Int,
                   callbacks!::Array{Function}, k::Int = 0)
  error("not yet implemented");
  temp_f = copy(sim.lat.f);
  for i = k+1:n_steps
    try

      sim_step!(sim.lat, temp_f, sim.msm, sbounds, collision_f!, cbounds, bcs!);

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
function simulate!(sim::FreeSurfSim, sbounds::Array{Int64,2},
                   collision_f!::Function, cbounds::Array{Int64,2},
                   bcs!::Array{Function}, n_steps::Int, k::Int = 0)
  error("not yet implemented");
  temp_f = copy(sim.lat.f);
  for i = k+1:n_steps
    try; sim_step!(sim.lat, temp_f, sim.msm, sbounds, collision_f!, cbounds, bcs!);
    catch e
      showerror(STDERR, e); println(); warn("Simulation interrupted at step $i !");
      return i;
    end
  end

  return n_steps;

end
