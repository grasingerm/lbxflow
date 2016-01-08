module LBXFlow

include("boundary.jl");
include(joinpath("col","collision.jl"));
include("convergence.jl");
include("lattice.jl");
include("lbxio.jl");
include("multiscale.jl");
include("profile.jl");
include(joinpath("sim","simulate.jl"));

using PyCall;
@pyimport yaml;

export simulate!, sim_step!;
export MultiscaleMap, map_to_macro!, strain_rate_tensor;
export north_bounce_back!, south_bounce_back!, east_bounce_back!,           ...
       west_bounce_back!, north_half_bounce_back!, south_half_bounce_back!, ...
       east_half_bounce_back!, 

end
