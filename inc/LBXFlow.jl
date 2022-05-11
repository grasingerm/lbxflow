# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Logging;
using LinearAlgebra;
using DelimitedFiles;

module LBXFlow

#TODO clean up simulate! code with some kernal functions...
macro _report_and_exit(err, i)
  return quote
    bt = catch_backtrace(); 
    showerror(stderr, $err, bt);
    @warn("Showing backtrace:");
    Base.show_backtrace(stderr, backtrace()); # display callstack
    @warn("Simulation interrupted at step ", $i, "!");
    return $i;
  end
end

@inline function _report_and_exit(err, i)
  bt = catch_backtrace(); 
  showerror(stderr, err, bt);
  @warn("Showing backtrace:");
  Base.show_backtrace(stderr, backtrace()); # display callstack
  @warn("Simulation interrupted at step ", i, "!");
  rethrow(err);
end

include("debug.jl");
include("numerics.jl");
include("lattice.jl");
include("multiscale.jl");
include(joinpath("sim","simtypes.jl"));
include(joinpath("sim","tracking.jl"));
include("entropy.jl");
include(joinpath("col","collision.jl"));
include(joinpath("sim","adapt.jl"));
include(joinpath("sim","simulate.jl"));
include("accessors.jl");
include("boundary.jl");
include("obstacle.jl");
include("convergence.jl");
include("lbxio.jl");
include("profile.jl");
include("stability.jl");
include("api.jl");
include("analytical.jl");

export vel_prof_acsr, vbar_prof_acsr, vel_mag_acsr, vel_field_acsr;
export density_acsr, pressure_acsr, mass_acsr, ff_acsr, streamlines_acsr;
export m2phase_acsr, fluid_frac_acsr;
export dumpsim_jld, loadsim_jld, dumpsim, latest_backup_dir, load_backup_dir;
export load_latest_backup, load_latest_jld, rrm, write_jld_file_callback;
export write_backup_file_callback, write_datafile_callback, take_snapshot_callback;
export sample_u_uniform_sqr, sample_v_uniform_sqr, sample_ρ_uniform_sqr;
export pause_sim_callback, print_step_callback, gc_callback;
export pyplot_callback, pycontour_callback, pyquiver_callback;
export pystream_callback;

end
