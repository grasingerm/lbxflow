# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

const __collision_root__ = dirname(@__FILE__);
require(abspath(joinpath(__collision_root__, "freecol.jl"))); # free surface collisions
require(abspath(joinpath(__collision_root__, "modcol.jl")));  # "modular" collisions
require(abspath(joinpath(__collision_root__, "pmodcol.jl"))); # parallel mod collisions
require(abspath(joinpath(__collision_root__, "stdcol.jl")));  # "standard" collisions
