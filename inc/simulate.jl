const __simulate_root__ = dirname(@__FILE__);
require(abspath(joinpath(__simulate_root__, "collision.jl")));
require(abspath(joinpath(__simulate_root__, "lattice.jl")));
require(abspath(joinpath(__simulate_root__, "multiscale.jl")));

#! Simulation object
#! lat Lattice
#! msm Multiscale map
immutable Sim
  lat::Lattice
  msm::MultiscaleMap

  Sim(lat::Lattice, msm::MultiscaleMap) = new(lat, msm);
end

#! stream particle densities
function stream!(lat::Lattice, temp_f::Array{Float64,3}, bounds::Array{Int64,2})

  const ni, nj = size(temp_f);
  const nbounds, = size(bounds);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];
    for i = i_min:i_max, j = j_min:j_max, k = 1:9
      i_new = i + lat.c[k,1];
      j_new = j + lat.c[k,2];

      if i_new > i_max || j_new > j_max || i_new < i_min || j_new < j_min
        continue;
      end

      temp_f[i_new,j_new,k] = lat.f[i,j,k];
    end
  end

  copy!(lat.f, temp_f);
end

#! Simulate a single step
function sim_step!(lat::Lattice, temp_f::Array{Float64,3}, msm::MultiscaleMap,
  sbounds::Array{Int64,2}, collision_f!::Function, cbounds::Array{Int64,2},
  bcs!::Array{Function})

  collision_f!(lat, msm);
  stream!(lat, temp_f);

  for bc! in bcs!
    bc!(lat);
  end

  map_to_macro!(lat, msm);

end

#! Run simulation
function simulate!(sim::Sim, sbounds::Array{Int64,2}, collision_f!::Function,
  cbounds::Array{Int64,2}, bcs!::Array{Function}, n_steps::Int,
  tests_for_term::Array{Function}, callbacks!::Array{Function}, k = 1)

  temp_f = copy(sim.lat.f);

  sim_step!(sim.lat, temp_f, sim.msm, sbounds, collision_f!, cbounds, bcs!);

  for c! in callbacks!
    c!(sim, k);
  end

  prev_msm = MultiscaleMap(sim.msm);

  for i = k+1:n_steps
    sim_step!(sim.lat, temp_f, sim.msm, collision_f!, bcs!);

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
  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::Sim, sbounds::Array{Int64,2}, collision_f!::Function,
  cbounds::Array{Int64,2}, bcs!::Array{Function}, n_steps::Int,
  callbacks!::Array{Function}, k = 1)

  temp_f = copy(sim.lat.f);
  for i = k:n_steps
    sim_step!(sim.lat, temp_f, sim.msm, sbounds, collision_f!, cbounds, bcs!);

    for c! in callbacks!
      c!(sim, i);
    end
  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::Sim, sbounds::Array{Int64,2}, collision_f!::Function,
  cbounds::Array{Int64,2}, bcs!::Array{Function}, n_steps::Int, k = 1)

  temp_f = copy(sim.lat.f);
  for i = k:n_steps
    sim_step!(sim.lat, temp_f, sim.msm, sbounds, collision_f!, cbounds, bcs!);
  end

  return n_steps;

end

# ------------------------- Deprecated functions in v0.2 ---------------------

# TODO: deprecate this
#! stream particle densities
function stream!(lat::Lattice, temp_f::Array{Float64,3})

  global VERSION;
  if VERSION > 0.1
    warn("This stream function is deprecated as of v0.2.");
  end

  const ni, nj = size(temp_f);

  #! Stream
  for i = 1:ni, j = 1:nj, k = 1:9
    i_new = i + lat.c[k,1];
    j_new = j + lat.c[k,2];

    if i_new > ni || j_new > nj || i_new < 1 || j_new < 1
      continue;
    end

    temp_f[i_new,j_new,k] = lat.f[i,j,k];
  end

  copy!(lat.f, temp_f);
end

# TODO: deprecate this
#! stream particle densities
function stream!(lat::Lattice, temp_f::Array{Float64,3}, is::UnitRange{Int64},
  js::UnitRange{Int64})

  global VERSION;
  if VERSION > 0.1
    warn("This stream function is deprecated.");
  end

  const ni, nj = size(temp_f);

  #! Stream
  for i = is, j = js, k = 1:9
    i_new = i + lat.c[k,1];
    j_new = j + lat.c[k,2];

    if i_new > ni || j_new > nj || i_new < 1 || j_new < 1
      continue;
    end

    temp_f[i_new,j_new,k] = lat.f[i,j,k];
  end

  copy!(lat.f, temp_f);
end

#! Simulate a single step
function sim_step!(lat::Lattice, temp_f::Array{Float64,3}, msm::MultiscaleMap,
  collision_f!::Function, bcs!::Array{Function}, stream_f!::Function)

  global VERSION;
  if VERSION > 0.1
    warn("This stream function is deprecated.");
  end

  collision_f!(lat, msm);
  stream_f!(lat, temp_f);

  for bc! in bcs!
    bc!(lat);
  end

  map_to_macro!(lat, msm);
end

#! Run simulation
function simulate!(sim::Sim, collision_f!::Function,
  bcs!::Array{Function}, n_steps::Int, test_for_term::Function,
  callbacks!::Array{Function}, stream_f!::Function)

  global VERSION;
  if VERSION > 0.1
    warn("This stream function is deprecated.");
  end

  temp_f = copy(sim.lat.f);

  sim_step!(sim.lat, temp_f, sim.msm, collision_f!, bcs!, stream_f!);

  for c! in callbacks!
    c!(sim, 1);
  end

  prev_msm = MultiscaleMap(sim.msm);

  for i = 2:n_steps
    sim_step!(sim.lat, temp_f, sim.msm, collision_f!, bcs!);

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
  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::Sim, collision_f!::Function,
  bcs!::Array{Function}, n_steps::Int, callbacks!::Array{Function},
  stream_f!::Function)

  global VERSION;
  if VERSION > 0.1
    warn("This stream function is deprecated.");
  end

  temp_f = copy(sim.lat.f);
  for i = 1:n_steps
    sim_step!(sim.lat, temp_f, sim.msm, collision_f!, bcs!, stream_f!);

    for c! in callbacks!
      c!(sim, i);
    end
  end

  return n_steps;

end

#! Run simulation
function simulate!(sim::Sim, collision_f!::Function,
  bcs!::Array{Function}, n_steps::Int, stream_f!::Function)

  global VERSION;
  if VERSION > 0.1
    warn("This stream function is deprecated.");
  end

  temp_f = copy(sim.lat.f);
  for i = 1:n_steps
    sim_step!(sim.lat, temp_f, sim.msm, collision_f!, bcs!, stream_f!);
  end

  return n_steps;

end
