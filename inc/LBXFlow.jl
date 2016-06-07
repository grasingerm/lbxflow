# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

module LBXFlow

#TODO clean up simulate! code with some kernal functions...
macro _report_and_exit(e, i)
  return quote
    const bt = catch_backtrace(); 
    showerror(STDERR, $e, bt);
    println();
    println("Showing backtrace:");
    Base.show_backtrace(STDERR, backtrace()); # display callstack
    println();
    warn("Simulation interrupted at step ", $i, "!");
    return $i;
  end
end

include("debug.jl");
include("numerics.jl");
include("lattice.jl");
include("multiscale.jl");
include(joinpath("sim","simtypes.jl"));
include(joinpath("sim","tracking.jl"));
include(joinpath("col","collision.jl"));
include(joinpath("sim","adapt.jl"));
include(joinpath("sim","simulate.jl"));
include("accessors.jl");
include("entropy.jl");
include("boundary.jl");
include("obstacle.jl");
include("convergence.jl");
include("lbxio.jl");
include("profile.jl");
include("stability.jl");
include("api.jl");
include("analytical.jl");

end
