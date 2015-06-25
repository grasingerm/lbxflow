const __lbxio_root__ = dirname(@__FILE__);
require(abspath(joinpath(__lbxio_root__, "multiscale.jl")));
require(abspath(joinpath(__lbxio_root__, "sim", "simulate.jl")));

import HDF5, JLD

# ===========================================================================
# ================================= Data output and loading =================
# ===========================================================================

#! Dump simulation to JLD file
function dumpsim_jld(datadir::String, sim::Sim, step::Int,
  append_step_to_name::Bool = false)
  
  const jldpath = (append_step_to_name) ? joinpath(datadir, "bak_$step.jld") :
                                          joinpath(datadir, "bak.jld");
  JLD.jldopen(jldpath, "w") do file
    write(file, "sim", sim);
    write(file, "step", step);
  end
end

#! Load simulation from JLD file
function loadsim_jld(path::String)
  d = JLD.jldopen(path, "r") do file
    read(file);
  end
  return d["step"], d["sim"];
end

#! Dump simulation to a text file
function dumpsim(datadir::String, sim::Sim, step::Int)
  const dumpdir = joinpath(datadir, "bak_$step");
  if !isdir(dumpdir)
    mkdir(dumpdir);
  end
 
  # backup particle distributions and lattice config
  nk, ni, nj = size(sim.lat.f);
  w = open(joinpath(dumpdir, "lat.dat"), "w");
  write(w, "ni: $ni\n");
  write(w, "nj: $nj\n");
  write(w, "nk: $nk\n");
  write(w, string("dx: ", sim.lat.dx, "\n"));
  write(w, string("dt: ", sim.lat.dt, "\n"));
  close(w);

  for k=1:nk
    writedlm(joinpath(dumpdir, "f$k.dat"), sim.lat.f[k,:,:]);
  end

  # backup multiscale map
  w = open(joinpath(dumpdir, "msm.dat"), "w");
  write(w, "ni: $ni\n");
  write(w, "nj: $nj\n");
  write(w, string("rho_0: ", sim.msm.rho_0), "\n");
  close(w);

  for d=1:2
    writedlm(joinpath(dumpdir, "u$d.dat"), sim.msm.u[d,:,:]);
  end
  writedlm(joinpath(dumpdir, "rho.dat"), sim.msm.rho);
  writedlm(joinpath(dumpdir, "omega.dat"), sim.msm.omega);
end

#! Find the latest backup directory 
function latest_backup_dir(datadir::String)
  step = 0;
  dir = "";

  dircontents = readdir(datadir);
  for f in dircontents
    if isdir(joinpath(datadir, f)) && beginswith(f, "bak")
      this_step = parse(split(f,"_")[2]);
      if this_step > step
        step = this_step;
        dir = f;
      end
    end
  end

  return step, joinpath(datadir, dir);
end

#! Load backup directory
function load_backup_dir(backup_dir::String)
  latdefs = Dict();
  msmdefs = Dict();

  for (defs, fname) in zip((latdefs, msmdefs),("lat.dat", "msm.dat"))
    try
      dat = open(joinpath(backup_dir, fname));
      datlines = readlines(dat);
      close(dat);
      for line in datlines
        k, v = [strip(s) for s in split(line, ":")];
        defs[k] = parse(v);
      end
    catch e
      return 0, nothing;
    end
  end

  lat = Lattice(latdefs["dx"], latdefs["dt"], latdefs["ni"], latdefs["nj"]);
  msm = MultiscaleMap(0.0, lat, msmdefs["rho_0"]);

  for k=1:lat.n
    try
      lat.f[k,:,:] = readdlm(joinpath(backup_dir,"f$k.dat"));
    catch e
      return 0, nothing;
    end
  end

  for d=1:2
    try
      msm.u[d,:,:] = readdlm(joinpath(backup_dir,"u$d.dat"));
    catch e
      return 0, nothing;
    end
  end

  try
    msm.rho = readdlm(joinpath(backup_dir, "rho.dat"));
  catch e
    return 0, nothing;
  end

  try
    msm.omega = readdlm(joinpath(backup_dir, "omega.dat"));
  catch e
    return 0, nothing;
  end

  return parse(split(backup_dir, "_")[end]), Sim(lat, msm);
end

#! Search data directory for latest backup information
function load_latest_backup(datadir::String)
  if isfile(joinpath(datadir, "bak.jld"))
    return loadsim_jld(joinpath(datadir, "bak.jld"))
  end

  step = 0;
  latest_path = "";
  paths = filter((p) -> beginswith(p, "bak"), readdir(datadir));
  for path in paths
    if contains(path, "_")
      this_step = parse(split(split(path, "_")[2], ".")[1]);
      if this_step > step
        step = this_step;
        latest_path = path;
      end
    end
  end


  const latest_fpath = joinpath(datadir, latest_path);
  if latest_path == ""; return 0, nothing;
  elseif isdir(latest_fpath)
    return load_backup_dir(latest_fpath);
  elseif isfile(latest_fpath) && endswith(latest_fpath, ".jld")
    return loadsim_jld(latest_fpath);
  end

  return 0, nothing; # no backup files found
end

#! Search data directory for latest backup information
function load_latest_jld(datadir::String)
  if isfile(joinpath(datadir, "bak.jld"))
    return loadsim_jld(joinpath(datadir, "bak.jld"))
  end

  step = 0;
  latest_path = "";
  paths = filter((p) -> beginswith(p, "bak") && endswith(".jld"),
                 readdir(datadir));
  for path in paths
    this_step = split(split(path, "_")[2], ".")[1];
    if this_step > step
      step = this_step;
      latest_path = path;
    end
  end

  if step == 0; return 0, nothing; end;
  return loadsim_jld(latest_path);
end

#! Recursively remove files and directories
function rrm(path::String)
  if isdir(path)
    for f in readdir(path); rrm(joinpath(path, f)); end
  end
  rm(path);
end
# ===========================================================================
# ===================== Sim callbacks: IO and visualization =================
# ===========================================================================

#! Create a callback function for writing a jld backup file
function write_jld_file_callback(datadir::String,
                                 append_step_to_name::Bool = false)
 
  return (sim::Sim, k::Int) -> begin
    dumpsim_jld(datadir, sim, k, append_step_to_name);
  end;
end

#! Create a callback function for writing a jld backup file
function write_jld_file_callback(datadir::String, stepout::Int,
                                 append_step_to_name::Bool = false)

  return (sim::Sim, k::Int) -> begin
    if k % stepout == 0
      dumpsim_jld(datadir, sim, k, append_step_to_name);
    end
  end;
end
#! Create a callback function for writing a backup file
function write_backup_file_callback(datadir::String)
  return (sim::Sim, k::Int) -> begin
    dumpsim(datadir, sim, k);
  end;
end

#! Create a callback function for writing a backup file
function write_backup_file_callback(datadir::String, stepout::Int)
  return (sim::Sim, k::Int) -> begin
    if k % stepout == 0
      dumpsim(datadir, sim, k);
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
    const nj = size(sim.msm.u, 3);
    x = Array(Float64, (nj, 3));

    for j=1:nj
      x[j,:] = [j, sim.msm.u[1,i,j], sim.msm.u[2,i,j]];
    end

    return x;
  end;

end

#! Extract ux profile cut parallel to y-axis
function extract_ux_prof_callback(i::Int)

  return (sim::Sim) -> begin
    const ni, nj = size(sim.msm.u, 2), size(sim.msm.u, 3);
    const u = vec(sim.msm.u[1,i,:]);

    x = linspace(-0.5, 0.5, nj);
    y = u;

    return [x y];
  end;

end

#! Extract u_bar profile cut parallel to y-axis
function extract_ubar_prof_callback(i::Int)

  return (sim::Sim) -> begin
    const ni, nj = size(sim.msm.u, 2), size(sim.msm.u, 3);
    const u = vec(sim.msm.u[1,i,:]);

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
                                  pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[1,i,:]);

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
                                    pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);
      const u = vec(sim.msm.u[1,i,:]);

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
                                    pause::FloatingPoint = 0.025)

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
                             pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int,
                                   pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho, 2), size(sim.msm.rho, 3);
      x = linspace(0.0, 1.0, ni);
      y = linspace(0.0, 1.0, nj);

      clf();
      streamplot(x, y, transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      ylim(0.0, 1.0);
      xlim(0.0, 1.0);
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int, fname::String,
                                  pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[1,i,:]);

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
                                    pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);
      const u = vec(sim.msm.u[1,i,:]);

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
                                    pause::FloatingPoint = 0.025)

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
                             pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int, fname::String,
                                   pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho, 2), size(sim.msm.rho, 3);
      x = linspace(0.0, 1.0, ni);
      y = linspace(0.0, 1.0, nj);

      clf();
      streamplot(x, y, transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      ylim(0.0, 1.0);
      xlim(0.0, 1.0);
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int,
                                  xy::(Number,Number),
                                  pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[1,i,:]);

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
                                    xy::(Number,Number),
                                    pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);
      const u = vec(sim.msm.u[1,i,:]);

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
                                    pause::FloatingPoint = 0.025)

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
                             pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int, xy::(Number,Number),
                                   pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = linspace(0.0, 1.0, ni);
      y = linspace(0.0, 1.0, nj);

      clf();
      streamplot(x, y,transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      ylim(0.0, 1.0);
      xlim(0.0, 1.0);
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int,
                                  xy::(Number,Number), fname::String,
                                  pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[1,i,:]);

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
                                    xy::(Number,Number), fname::String,
                                    pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);
      const u = vec(sim.msm.u[1,i,:]);

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
                                    fname::String, pause::FloatingPoint = 0.025)

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
                             fname::String, pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int, xy::(Number,Number),
                                   fname::String, pause::FloatingPoint = 0.025)

  return (sim::Sim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho, 2), size(sim.msm.rho, 3);
      x = linspace(0.0, 1.0, ni);
      y = linspace(0.0, 1.0, nj);

      clf();
      streamplot(x, y,transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      ylim(0.0, 1.0);
      xlim(0.0, 1.0);
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end
