const __dump_bak_root__ = dirname(@__FILE__);
require(abspath(joinpath(__dump_bak_root__, "..", "inc", "lattice.jl")));
require(abspath(joinpath(__dump_bak_root__, "..", "inc", "multiscale.jl")));
require(abspath(joinpath(__dump_bak_root__, "..", "inc", "sim", "simtypes.jl")));

using HDF5, JLD;

if length(ARGS) != 1
  info("Usage: dump_bak.jl path_to_bak.jld");
  exit(1);
end

fname = ARGS[1];
data = jldopen(fname,"r") do file
  read(file);
end
println(data);
