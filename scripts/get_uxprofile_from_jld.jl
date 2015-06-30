const __get_uxprofile_from_jld_root__ = dirname(@__FILE__);
require(abspath(joinpath(__get_uxprofile_from_jld_root__, "..", "inc", "lattice.jl")));
require(abspath(joinpath(__get_uxprofile_from_jld_root__, "..", "inc", "multiscale.jl")));
require(abspath(joinpath(__get_uxprofile_from_jld_root__, "..", "inc", "sim", "simtypes.jl")));

using HDF5, JLD;

if length(ARGS) != 2
  info("Usage: get_uxprofile_from_jld.jl fname i_idx");
  exit(1);
end

fname, i_idx = ARGS;
data = jldopen(fname,"r") do file
  read(file);
end
println(data);
