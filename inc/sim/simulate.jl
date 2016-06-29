# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

include("recolor.jl");

#! Stream particle densities
#!
#! \param lat Lattice to stream on
#! \param temp_f Temp lattice to store streamed particle distributions on
#! \param bounds Boundaries enclosing active streaming regions
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
      if t.state[i,j] != GAS
        for k = 1:lat.n
          i_new = i + lat.c[1,k];
          j_new = j + lat.c[2,k];

          if (i_new > i_max || j_new > j_max || i_new < i_min || j_new < j_min 
              || t.state[i_new,j_new] == GAS)
            continue;
          end
          temp_f[k,i_new,j_new] = lat.f[k,i,j];
        end
      end
    end
  end

  copy!(lat.f, temp_f);
end

#! Stream particle densities around obstacles
#!
#! \param   lat             Lattice to stream on
#! \param   temp_f          Temp lattice to store streamed particle dists
#! \param   active_cells    Matrix of active flags
function stream!(lat::Lattice, temp_f::Array{Float64,3}, 
                 active_cells::Matrix{Bool})
  const ni, nj = size(lat.f, 2), size(lat.f, 3);

  #! Stream
  for j = 1:nj, i = 1:ni, k = 1:lat.n
    i_new = i + lat.c[1,k];
    j_new = j + lat.c[2,k];

    if (i_new <= ni && j_new <= nj && i_new >= 1 && j_new >= 1 
        && active_cells[i_new, j_new])
      temp_f[k,i_new,j_new] = lat.f[k,i,j];
    end
  end

  copy!(lat.f, temp_f);
end

# TODO use same `inbounds` function or macro throughout code
#! Stream particle densities in free surface conditions
#!
#! \param   lat             Lattice to stream on
#! \param   temp_f          Temp lattice to store streamed particle dists
#! \param   active_cells    Matrix of active flags
#! \param   t               Mass tracker
function stream!(lat::Lattice, temp_f::Array{Float64,3}, 
                 active_cells::Matrix{Bool}, t::Tracker)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);

  #! Stream
  for j = 1:nj, i = 1:ni, k = 1:lat.n
    i_new = i + lat.c[1, k];
    j_new = j + lat.c[2, k];

    if (i_new <= ni && j_new <= nj && i_new >= 1 && j_new >= 1 
        && active_cells[i_new, j_new] && t.state[i_new, j_new] != GAS)
      temp_f[k,i_new,j_new] = lat.f[k,i,j];
    end
  end

  copy!(lat.f, temp_f);
end

#! Simulate a single step
function sim_step!(sim::Sim,
                   sbounds::Matrix{Int64}, collision_f!::LBXFunction, 
                   cbounds::Matrix{Int64}, bcs!::Vector{LBXFunction})
  lat = sim.lat;
  temp_f = copy(lat.f);
  msm = sim.msm;

  collision_f!(sim, cbounds);
  stream!(lat, temp_f, sbounds);

  for bc! in bcs!
    bc!(sim);
  end

  map_to_macro!(lat, msm);
end

#! Simulate a step for immiscible two-phase flow
function sim_step!(sim::M2PhaseSim,
                   sbounds::Matrix{Int64}, col_f!::M2PhaseColFunction, 
                   cbounds::Matrix{Int64}, bcs!::Vector{LBXFunction})
 
  # forward collision function calls
  println("collision: ");
  col_f!(sim, cbounds);
  println("Is there an NaN?", true in isnan(sim.simr.lat.f));
  readline(STDIN);

  println("recolor: ");
  recolor!(sim, sbounds, col_f!.col_fr!.feq_f, col_f!.col_fb!.feq_f);
  println("Is there an NaN?", true in isnan(sim.simr.lat.f));
  readline(STDIN);

  println("stream: ");
  temp_f = copy(sim.simr.lat.f);
  stream!(sim.simr.lat, temp_f, sbounds);
  stream!(sim.simb.lat, temp_f, sbounds);
  println("Is there an NaN?", true in isnan(sim.simr.lat.f));
  readline(STDIN);

  println("bcs: ");
  for bc! in bcs!
    bc!(sim.simr);
    bc!(sim.simb);
  end
  println("Is there an NaN?", true in isnan(sim.simr.lat.f));
  readline(STDIN);

  println("m2m: ");
  map_to_macro!(sim.simr.lat, sim.simr.msm);
  map_to_macro!(sim.simb.lat, sim.simb.msm);
  println("Is there an NaN?", true in isnan(sim.simr.lat.f));
  readline(STDIN);
end

#! Simulate a single step free surface flow step
function sim_step!(sim::FreeSurfSim,
                   sbounds::Matrix{Int64}, collision_f!::ColFunction, 
                   cbounds::Matrix{Int64}, 
                   bcs!::Vector{LBXFunction})
  lat               =   sim.lat;
  temp_f            =   copy(lat.f);
  msm               =   sim.msm;
  t                 =   sim.tracker;
  unorms            =   Dict{Tuple{Int, Int}, Vector{Float64}}();

  @_checkdebug_mass_cons("whole step", t.M, begin
  # Algorithm should be:
  # 1.  mass transfer
  @_checkdebug_mass_cons("masstransfer!", t.M, masstransfer!(sim, sbounds), 1e-9);

  # 2.  stream
  @_checkdebug_mass_cons("stream!", t.M, stream!(lat, temp_f, sbounds, t), 1e-9);

  # 3.  reconstruct distribution functions from empty cells
  # 4.  reconstruct distribution functions along interface normal
  @_checkdebug_mass_cons("f_reconst!", t.M, for (i, j) in t.interfacels
    unorms[(i, j)] = f_reconst!(sim, t, (i, j), collision_f!.feq_f, sim.rho_g);
  end, 1e-9);

  # 5.  particle collisions
  @_checkdebug_mass_cons("collision_f!", t.M, collision_f!(sim, cbounds), 1e-9);
  
  # 6.  enforce boundary conditions
  @_checkdebug_mass_cons("bcs!", t.M, for bc! in bcs!
    bc!(sim);
  end, 1e-9);

  # 7.  calculate macroscopic variables
  @_checkdebug_mass_cons("map_to_macro!", t.M, map_to_macro!(lat, msm), 1e-9);

  # 8.  update fluid fractions
  # 9.  update cell states
  @_checkdebug_mass_cons("update!", t.M, 
    update!(sim, collision_f!.feq_f, unorms), 1e-9);
  end, 1e-9);
end

# Adaptive time step simulation step
function sim_step!(sim::AdaptiveTimeStepSim,
                   sbounds::Matrix{Int64}, collision_f!::ColFunction, args...)
  sim_step!(sim.isim, sbounds, collision_f!, args...); # forward sim_step!
  adapt_time_step!(sim, collision_f!);
end

#! Simulate a single step
function sim_step!(sim::Sim,
                   collision_f!::LBXFunction, active_cells::Matrix{Bool}, 
                   bcs!::Vector{LBXFunction})
  lat = sim.lat;
  temp_f = copy(lat.f);
  msm = sim.msm;

  collision_f!(sim, active_cells);
  stream!(lat, temp_f, active_cells);

  for bc! in bcs!
    bc!(sim);
  end

  map_to_macro!(lat, msm);
end

#! Simulate a step for immiscible two-phase flow
function sim_step!(sim::M2PhaseSim,
                   col_f!::M2PhaseColFunction, 
                   active_cells::Matrix{Bool}, bcs!::Vector{LBXFunction})
 
  # forward collision function calls
  col_f!(sim, active_cells);

  recolor!(sim, active_cells, col_fr!.feq_f, col_fb!.feq_f);

  temp_f = copy(sim.simr.lat.f);
  stream!(sim.simr.lat, temp_f, active_cells);
  stream!(sim.simb.lat, temp_f, active_cells);

  for bc! in bcs!
    bc!(sim.simr);
    bc!(sim.simb);
  end

  map_to_macro!(sim.simr.lat, sim.simr.msm);
  map_to_macro!(sim.simb.lat, sim.simb.msm);
end

#! Simulate a single step free surface flow step
function sim_step!(sim::FreeSurfSim,
                   collision_f!::ColFunction, active_cells::Matrix{Bool}, 
                   bcs!::Vector{LBXFunction})
  lat               =   sim.lat;
  temp_f            =   copy(lat.f);
  msm               =   sim.msm;
  t                 =   sim.tracker;
  unorms            =   Dict{Tuple{Int, Int}, Vector{Float64}}();

  @_checkdebug_mass_cons("whole step", t.M, begin
  # Algorithm should be:
  # 1.  mass transfer
  @_checkdebug_mass_cons("masstransfer!", t.M, masstransfer!(sim, active_cells), 1e-9);

  # 2.  stream
  @_checkdebug_mass_cons("stream!", t.M, stream!(lat, temp_f, active_cells, t), 1e-9);

  # 3.  reconstruct distribution functions from empty cells
  # 4.  reconstruct distribution functions along interface normal
  @_checkdebug_mass_cons("f_reconst!", t.M,
  for (i, j) in t.interfacels #TODO maybe abstract out interface list...
    unorms[(i, j)] = f_reconst!(sim, t, (i, j), collision_f!.feq_f, sim.rho_g);
  end, 1e-9);

  # 5.  particle collisions
  @_checkdebug_mass_cons("collision_f!", t.M, collision_f!(sim, active_cells), 1e-9);
  
  # 6.  enforce boundary conditions
  @_checkdebug_mass_cons("bcs!", t.M, for bc! in bcs!
    bc!(sim);
  end, 1e-9);

  # 7.  calculate macroscopic variables
  @_checkdebug_mass_cons("map_to_macro!", t.M, map_to_macro!(lat, msm), 1e-9);

  # 8.  update fluid fractions
  # 9.  update cell states
  @_checkdebug_mass_cons("update!", t.M, 
    update!(sim, collision_f!.feq_f, unorms), 1e-9);
  end, 1e-9);
end

# Adaptive time step simulation step
function sim_step!(sim::AdaptiveTimeStepSim,
                   collision_f!::ColFunction, args...)
  sim_step!(sim.isim, collision_f!, args...); # forward sim_step!
  adapt_time_step!(sim, collision_f!);
end

#! Run simulation
function simulate!(sim::AbstractSim, sbounds::Matrix{Int64},
                   collision_f!::LBXFunction, cbounds::Matrix{Int64},
                   bcs!::Vector{LBXFunction}, n_steps::Real, 
                   test_for_term::LBXFunction,
                   callbacks!::Vector{LBXFunction}, k::Real = 0)

  sim_step!(sim, sbounds, collision_f!, cbounds, bcs!);

  for c! in callbacks!
    c!(sim, (k+=sim.Δt));
  end

  prev_msm = MultiscaleMap(sim.msm);

  while (k+=sim.Δt) <= n_steps
    try

      sim_step!(sim, sbounds, collision_f!, cbounds, bcs!);

      for c! in callbacks!
        c!(sim, k);
      end

      # if returns true, terminate simulation
      if test_for_term(sim.msm, prev_msm)
        return k;
      end

      copy!(prev_msm.omega, sim.msm.omega);
      copy!(prev_msm.rho, sim.msm.rho);
      copy!(prev_msm.u, sim.msm.u);
    
    catch e

      @_report_and_exit(e, k);

    end
  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::AbstractSim, sbounds::Matrix{Int64},
                   collision_f!::LBXFunction, cbounds::Matrix{Int64},
                   bcs!::Vector{LBXFunction}, n_steps::Real, 
                   test_for_term::LBXFunction,
                   steps_for_term::Int, callbacks!::Vector{LBXFunction}, 
                   k::Real = 0)

  sim_step!(sim, sbounds, collision_f!, cbounds, bcs!);

  for c! in callbacks!
    c!(sim, (k+=sim.Δt));
  end

  prev_msms = Vector{MultiscaleMap}(steps_for_term);
  for i=1:steps_for_term; prev_msms[i] = MultiscaleMap(sim.msm); end;

  idx = 0;
  while (k+=sim.Δt) <= n_steps
    try

      sim_step!(sim, sbounds, collision_f!, cbounds, bcs!);

      for c! in callbacks!
        c!(sim, k);
      end

      # if returns true, terminate simulation
      if test_for_term(sim.msm, prev_msms)
        return k;
      end
  
      if (idx += 1) > steps_for_term; idx = 1; end
      copy!(prev_msms[idx].omega, sim.msm.omega);
      copy!(prev_msms[idx].rho,   sim.msm.rho);
      copy!(prev_msms[idx].u,     sim.msm.u);
    
    
    catch e

      @_report_and_exit(e, k);

    end

  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::AbstractSim, sbounds::Matrix{Int64},
                   collision_f!::LBXFunction, cbounds::Matrix{Int64},
                   bcs!::Vector{LBXFunction}, n_steps::Real,
                   callbacks!::Vector{LBXFunction}, k::Real = 0)

  while (k+=sim.Δt) <= n_steps
    try

      sim_step!(sim, sbounds, collision_f!, cbounds, bcs!);

      for c! in callbacks!
        c!(sim, k);
      end

    
    catch e

      @_report_and_exit(e, k);

    end

  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::AbstractSim, sbounds::Matrix{Int64},
                   collision_f!::LBXFunction, cbounds::Matrix{Int64},
                   bcs!::Vector{LBXFunction}, n_steps::Real, k::Real = 0)

  while (k+=sim.Δt) <= n_steps
    try; sim_step!(sim, sbounds, collision_f!, cbounds, bcs!);
    catch e

      @_report_and_exit(e, k);

    end
  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::AbstractSim, collision_f!::LBXFunction,
                   active_cells::Matrix{Bool},
                   bcs!::Vector{LBXFunction}, n_steps::Real, 
                   test_for_term::LBXFunction,
                   callbacks!::Vector{LBXFunction}, k::Real = 0)

  sim_step!(sim, collision_f!, active_cells, bcs!);

  for c! in callbacks!
    c!(sim, (k+=sim.Δt));
  end

  prev_msm = MultiscaleMap(sim.msm);

  while (k+=sim.Δt) <= n_steps
    try

      sim_step!(sim, collision_f!, active_cells, bcs!);

      for c! in callbacks!
        c!(sim, k);
      end

      # if returns true, terminate simulation
      if test_for_term(sim.msm, prev_msm)
        return k;
      end

      copy!(prev_msm.omega, sim.msm.omega);
      copy!(prev_msm.rho, sim.msm.rho);
      copy!(prev_msm.u, sim.msm.u);
    
    catch e

      @_report_and_exit(e, k);

    end
  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::AbstractSim,
                   collision_f!::LBXFunction, active_cells::Matrix{Bool},
                   bcs!::Vector{LBXFunction}, n_steps::Real, 
                   test_for_term::LBXFunction,
                   steps_for_term::Int, callbacks!::Vector{LBXFunction}, 
                   k::Real = 0)

  sim_step!(sim, collision_f!, active_cells, bcs!);

  for c! in callbacks!
    c!(sim, (k+=sim.Δt));
  end

  prev_msms = Vector{MultiscaleMap}(steps_for_term);
  for i=1:steps_for_term; prev_msms[i] = MultiscaleMap(sim.msm); end;

  idx = 0;
  while (k+=sim.Δt) <= n_steps
    try

      sim_step!(sim, collision_f!, active_cells, bcs!);

      for c! in callbacks!
        c!(sim, k);
      end

      # if returns true, terminate simulation
      if test_for_term(sim.msm, prev_msms)
        return k;
      end
  
      if (idx += 1) > steps_for_term; idx = 1; end
      copy!(prev_msms[idx].omega, sim.msm.omega);
      copy!(prev_msms[idx].rho,   sim.msm.rho);
      copy!(prev_msms[idx].u,     sim.msm.u);
    
    catch e

      @_report_and_exit(e, k);

    end
  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::AbstractSim,
                   collision_f!::LBXFunction, active_cells::Matrix{Bool},
                   bcs!::Vector{LBXFunction}, n_steps::Real,
                   callbacks!::Vector{LBXFunction}, k::Real = 0)

  while (k+=sim.Δt) <= n_steps
    try

      sim_step!(sim, collision_f!, active_cells, bcs!);

      for c! in callbacks!
        c!(sim, k);
      end

    catch e

      @_report_and_exit(e, k);

    end
  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::AbstractSim,
                   collision_f!::LBXFunction, active_cells::Matrix{Bool},
                   bcs!::Vector{LBXFunction}, n_steps::Real, k::Real = 0)

  while (k+=sim.Δt) <= n_steps
    try; sim_step!(sim, collision_f!, active_cells, bcs!);
    catch e

      @_report_and_exit(e, k);

    end
  end

  return n_steps;

end
