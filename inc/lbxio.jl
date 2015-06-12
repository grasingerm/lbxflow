const __lbxio_root__ = dirname(@__FILE__);
require(abspath(joinpath(__lbxio_root__, "multiscale.jl")));
require(abspath(joinpath(__lbxio_root__, "simulate.jl")));

function dumpsim(datadir::String, sim::Sim, step::Int)
  const dumpdir = joinpath(datadir, "bak_$step");
  if !isdir(dumpdir)
    mkdir(dumpdir);
  end
 
  # backup particle distributions and lattice config
  ni, nj, nk = size(sim.lat.f);
  w = open(joinpath(dumpdir, "lat.dat"), "w");
  write(w, "ni: $ni\n");
  write(w, "nj: $nj\n");
  write(w, "nk: $nk\n");
  write(w, string("dx: ", sim.lat.dx, "\n"));
  write(w, string("dt: ", sim.lat.dt, "\n"));
  close(w);

  for k=1:nk
    writedlm(joinpath(dumpdir, "f$k.dat"), sim.lat.f[:,:,k]);
  end

  # backup multiscale map
  w = open(joinpath(dumpdir, "msm.dat"), "w");
  write(w, "ni: $ni\n");
  write(w, "nj: $nj\n");
  write(w, string("dx: ", sim.msm.dx, "\n"));
  write(w, string("dt: ", sim.msm.dt, "\n"));
  write(w, string("rho_0: ", sim.msm.rho_0), "\n");
  close(w);

  for d=1:2
    writedlm(joinpath(dumpdir, "u$d.dat"), sim.msm.u[:,:,d]);
  end
  writedlm(joinpath(dumpdir, "rho.dat"), sim.msm.rho);
  writedlm(joinpath(dumpdir, "omega.dat"), sim.msm.omega);
end

# load simulation from backup dir
function loadsim(datadir)
  error("Function not yet implemented.");
end

#! Create a callback function for writing a backup file
function write_backup_file_callback(datadir::String)
  return (sim::Sim, k::Int) -> begin
    dumpsim(datadir, sim, k)
  end;
end

#! Create a callback function for writing a backup file
function write_backup_file_callback(datadir::String, stepout::Int)
  return (sim::Sim, k::Int) -> begin
    if k % stepout == 0
      dumpsim(datadir, sim, k)
    end
  end;
end


#! Create a callback function for writing to a delimited file
#!
#! \param pre Prefix of filename
#! \param stepout Number of steps in between writing
#! \param A Function for extracting a 2D array from the sim
#! \param delim Delimiter to separate values with
function write_datafile_callback (pre::String, stepout::Int, A::Function,
  dir=".", delim=',')

  return (sim::Sim, k::Int) -> begin
    if k % stepout == 0
      writedlm(joinpath(dir, pre*"_step-$k.dsv"), A(sim), delim);
    end
  end;

end

#! Extract velocity profile cut parallel to y-axis
function extract_prof_callback(i::Int)

  return (sim::Sim) -> begin
    const nj = size(sim.msm.u)[2];
    x = Array(Float64, (nj, 3));

    for j=1:nj
      x[j,:] = [j, sim.msm.u[i,j,1], sim.msm.u[i,j,2]];
    end

    return x;
  end;

end

#! Extract ux profile cut parallel to y-axis
function extract_ux_prof_callback(i::Int)

  return (sim::Sim) -> begin
    const ni, nj = size(sim.msm.u);
    const u = vec(sim.msm.u[i,:,1]);

    x = linspace(-0.5, 0.5, nj);
    y = u;

    return [x y];
  end;

end

#! Extract u_bar profile cut parallel to y-axis
function extract_ubar_prof_callback(i::Int)

  return (sim::Sim) -> begin
    const ni, nj = size(sim.msm.u);
    const u = vec(sim.msm.u[i,:,1]);

    x = linspace(-0.5, 0.5, nj);
    y = u / maximum(u);

    return [x y];
  end;

end

#! Create callback for reporting step
function print_step_callback(step::Int)
  return (sim::Sim, k::Int) -> begin
    if k % step == 0
      println("step $k");
    end
  end
end

#! Create callback for reporting step
function print_step_callback(step::Int, name::String)
  return (sim::Sim, k::Int) -> begin
    if k % step == 0
      println(name * ":\tstep $k");
    end
  end
end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int,
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u)[2];

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[i,:,1]);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux (lat / sec)");
      sleep(pause);
    end
  end

end

#! Plot nondimensional x-component of velocity profile cut parallel to y-axis
function plot_ubar_profile_callback(i::Int, iters_per_frame::Int,
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u)[2];
      const u = vec(sim.msm.u[i,:,1]);

      x = linspace(-0.5, 0.5, nj);
      y = u / maximum(u);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux / u_max");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int,
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(u_mag(sim.msm)));
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int,
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(transpose(sim.msm.u[:,:,1]), transpose(sim.msm.u[:,:,2]));
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int,
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = linspace(0.0, 1.0, ni);
      y = linspace(0.0, 1.0, nj);

      clf();
      streamplot(x, y, transpose(sim.msm.u[:,:,1]), transpose(sim.msm.u[:,:,2]));
      ylim(0.0, 1.0);
      xlim(0.0, 1.0);
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int, fname::String,
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u)[2];

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[i,:,1]);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux (lat / sec)");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot nondimensional x-component of velocity profile cut parallel to y-axis
function plot_ubar_profile_callback(i::Int, iters_per_frame::Int, fname::String,
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u)[2];
      const u = vec(sim.msm.u[i,:,1]);

      x = linspace(-0.5, 0.5, nj);
      y = u / maximum(u);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux / u_max");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int, fname::String,
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(u_mag(sim.msm)));
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int, fname::String,
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(transpose(sim.msm.u[:,:,1]), transpose(sim.msm.u[:,:,2]));
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int, fname::String,
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = linspace(0.0, 1.0, ni);
      y = linspace(0.0, 1.0, nj);

      clf();
      streamplot(x, y, transpose(sim.msm.u[:,:,1]), transpose(sim.msm.u[:,:,2]));
      ylim(0.0, 1.0);
      xlim(0.0, 1.0);
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int,
  xy::(Number,Number), pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u)[2];

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[i,:,1]);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux (lat / sec)");
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot nondimensional x-component of velocity profile cut parallel to y-axis
function plot_ubar_profile_callback(i::Int, iters_per_frame::Int,
  xy::(Number,Number), pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u)[2];
      const u = vec(sim.msm.u[i,:,1]);

      x = linspace(-0.5, 0.5, nj);
      y = u / maximum(u);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux / u_max");
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int, xy::(Number,Number),
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(u_mag(sim.msm)));
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int, xy::(Number,Number),
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(transpose(sim.msm.u[:,:,1]), transpose(sim.msm.u[:,:,2]));
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int, xy::(Number,Number),
  pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = linspace(0.0, 1.0, ni);
      y = linspace(0.0, 1.0, nj);

      clf();
      streamplot(x, y,transpose(sim.msm.u[:,:,1]), transpose(sim.msm.u[:,:,2]));
      ylim(0.0, 1.0);
      xlim(0.0, 1.0);
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int,
  xy::(Number,Number), fname::String, pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u)[2];

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[i,:,1]);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux (lat / sec)");
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot nondimensional x-component of velocity profile cut parallel to y-axis
function plot_ubar_profile_callback(i::Int, iters_per_frame::Int,
  xy::(Number,Number), fname::String, pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u)[2];
      const u = vec(sim.msm.u[i,:,1]);

      x = linspace(-0.5, 0.5, nj);
      y = u / maximum(u);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux / u_max");
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int, xy::(Number,Number),
  fname::String, pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(u_mag(sim.msm)));
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int, xy::(Number,Number),
  fname::String, pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(transpose(sim.msm.u[:,:,1]), transpose(sim.msm.u[:,:,2]));
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int, xy::(Number,Number),
  fname::String, pause::FloatingPoint = 0.1)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = linspace(0.0, 1.0, ni);
      y = linspace(0.0, 1.0, nj);

      clf();
      streamplot(x, y,transpose(sim.msm.u[:,:,1]), transpose(sim.msm.u[:,:,2]));
      ylim(0.0, 1.0);
      xlim(0.0, 1.0);
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#=
#! Animate x-component of velocity profile cut parallel to y-axis
function animate_ux_profile_callback(i::Int, iters_per_frame::Int,
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(msm.u)[2];

      x = linspace(-0.5, 0.5, nj);
      y = vec(msm.u[i,:,1]);

      plot(x,y);
      xlabel("x / width");
      ylabel("ux (lat / sec)");
      sleep(pause);
    end
  end

end

#! Animate nondimensional x-component of velocity profile cut parallel to y-axis
function animate_ubar_profile_callback(i::Int, iters_per_frame::Int,
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(msm.u)[2];

      x = linspace(-0.5, 0.5, nj);
      y = vec(msm.u[i,:,1]) / max(msm.u[i,:,1]);

      plot(x,y);
      xlabel("x / width");
      ylabel("ux / u_max");
      sleep(pause);
    end
  end

end

#! Animate x-component of velocity profile cut parallel to y-axis
function animate_umag_contour_callback(iters_per_frame::Int,
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      contour(u_mag(msm));
      sleep(pause);
    end
  end

end
=#
