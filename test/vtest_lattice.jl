const __vtl_root__ = dirname(@__FILE__);
require(abspath(joinpath(__vtl_root__, "..", "inc", "lattice.jl")));

println(Lattice.c);
println(Lattice.w);
