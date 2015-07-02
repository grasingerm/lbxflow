const __collision_root__ = dirname(@__FILE__);
require(abspath(joinpath(__collision_root__, "freecol.jl"))); # free surface collisions
require(abspath(joinpath(__collision_root__, "modcol.jl")));  # "modular" collisions
require(abspath(joinpath(__collision_root__, "pmodcol.jl"))); # parallel mod collisions
require(abspath(joinpath(__collision_root__, "stdcol.jl")));  # "standard" collisions
