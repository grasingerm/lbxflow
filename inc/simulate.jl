__simulate_root__ = dirname(@__FILE__);
require(abspath(joinpath(__simulate_root__, "collision.jl")));
require(abspath(joinpath(__simulate_root__, "lattice.jl")));
require(abspath(joinpath(__simulate_root__, "multiscale.jl")));

#! stream particle densities
function stream!(lat::Lattice, temp_f::Array{Float64,3})
  ni, nj = size(temp_f);

  #! Stream
  for i = 1:ni, j = 1:nj, k = 1:9
    i_new = i + lat.c[k,1];
    j_new = j + lat.c[k,2];

    if i_new > ni || j_new > nj || i_new < 1 || j_new < 1
      continue;
    end

    temp_f[i_new,j_new,k] = lat.f[i,j,k];
  end

  lat.f = copy(temp_f);
end

#! Run simulation
function simulate!(lat::Lattice, msm::MultiscaleMap, collision_f!::Function,
  bcs!::Array{Function}, n_steps::Int, callbacks!::Array{Function})

  temp_f = copy(lat.f);
  for i = 1:n_steps
    sim_step!(lat, temp_f, msm, collision_f!, bcs!);

    for c! in callbacks!
      c!(msm, i);
    end
  end

  return n_steps;

end

#! Run simulation
function simulate!(lat::Lattice, msm::MultiscaleMap, collision_f!::Function,
  bcs!::Array{Function}, n_steps::Int)

  temp_f = copy(lat.f);
  for i = 1:n_steps
    sim_step!(lat, temp_f, msm, collision_f!, bcs!);
  end

  return n_steps;

end

#! Simulate a single step
function sim_step!(lat::Lattice, temp_f::Array{Float64,3}, msm::MultiscaleMap,
  collision_f!::Function, bcs!::Array{Function})

  collision_f!(lat, msm);
  stream!(lat, temp_f);

  for bc! in bcs!
    bc!(lat);
  end

  map_to_macro!(lat, msm);

end
