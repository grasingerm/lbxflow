# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

module LBXFlow

include("debug.jl");
include("lattice.jl");
include("multiscale.jl");
include(joinpath("sim","simtypes.jl"));
include(joinpath("sim","tracking.jl"));
include(joinpath("sim","simulate.jl"));
include("entropy.jl");
include(joinpath("col","collision.jl"));
include("boundary.jl");
include("obstacle.jl");
include("convergence.jl");
include("lbxio.jl");
include("profile.jl");
include("stability.jl");
include("api.jl");
include("analytical.jl");

end
