# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

include("forcing.jl");
include("equilibrium.jl");
include("constitutive.jl");
include("mrt_matrices.jl");
include("freecol.jl"); # free surface collisions
include("modcol.jl");  # "modular" collisions
include("pmodcol.jl"); # parallel mod collisions
include("bgk.jl");
include("mrt.jl");
include("entropic.jl");
include("filtering.jl");
