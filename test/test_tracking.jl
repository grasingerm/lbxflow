using Base.Test
const __test_tracking_root__ = dirname(@__FILE__);
require(abspath(joinpath(__test_tracking_root__, "..", "inc", "sim", 
                         "tracking.jl")));
require(abspath(joinpath(__test_tracking_root__, "..", "inc", "sim", 
                         "simtypes.jl")));

const abounds = [1 1 1 1;]'
const bbounds = [1 10 2 10;
                 2 10 3 10]';
const cbounds = [-12 3 5 55;
                 -10 4 11 34;
                 -4  9 9  36;]'

println("Testing `inbounds`");

@test inbounds(1, 1, abounds) == true
@test inbounds(1, 2, abounds) == false
@test inbounds(2, 1, abounds) == false
@test inbounds(-1, 1, abounds) == false
@test inbounds(1, -1, abounds) == false
@test inbounds(3, 4, bbounds) == true
@test inbounds(2, 10, bbounds) == true
@test inbounds(5, 8, bbounds) == true
@test inbounds(8, 5, bbounds) == true
@test inbounds(2, 3, bbounds) == true
@test inbounds(10, 10, bbounds) == true
@test inbounds(11, 10, bbounds) == false
@test inbounds(10, 11, bbounds) == false
@test inbounds(21, 10, bbounds) == false
@test inbounds(-1, 10, bbounds) == false
@test inbounds(1, -10, bbounds) == false
@test inbounds(-3, 11, cbounds) == true
@test inbounds(-4, 34, cbounds) == true
@test inbounds(-10, 12, cbounds) == false
@test inbounds(-1, 12, cbounds) == true
@test inbounds(3, 15, cbounds) == true
@test inbounds(-2, 5, cbounds) == false
@test inbounds(-100, 400, cbounds) == false
@test inbounds(0, 12, cbounds) == true
@test inbounds(0, 0, cbounds) == false

println("Tests passed.");

#=
println("Testing `masstransfer!`");
lat = D2Q9Lattice(
=#
