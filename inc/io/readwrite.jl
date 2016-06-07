# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Dump simulation to JLD file
function dumpsim_jld(datadir::AbstractString, sim::AbstractSim, step::Real,
                     append_step_to_name::Bool = false)
  
  const jldpath = (append_step_to_name) ? joinpath(datadir, "bak_$step.jld") :
                                          joinpath(datadir, "bak.jld");
  JLD.jldopen(jldpath, "w") do file
    write(file, "sim", sim);
    write(file, "step", step);
  end
end

#! Load simulation from JLD file
function loadsim_jld(path::AbstractString)
  d = JLD.jldopen(path, "r") do file
    read(file);
  end
  return d["step"], d["sim"];
end

#! Dump simulation to a text file
function dumpsim(datadir::AbstractString, sim::Sim, step::Real)
  const dumpdir = joinpath(datadir, "bak_$step");
  if !isdir(dumpdir)
    mkpath(dumpdir);
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
function latest_backup_dir(datadir::AbstractString)
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
function load_backup_dir(backup_dir::AbstractString)
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
function load_latest_backup(datadir::AbstractString)
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
function load_latest_jld(datadir::AbstractString)
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
function rrm(path::AbstractString)
  if isdir(path)
    for f in readdir(path); rrm(joinpath(path, f)); end
  end
  rm(path);
end
# ===========================================================================
# =============================== Sim callbacks =============================
# ===========================================================================

#! Create a callback function for writing a jld backup file
function write_jld_file_callback(datadir::AbstractString,
                                 append_step_to_name::Bool = false)
 
  return (sim::AbstractSim, k::Real) -> begin
    dumpsim_jld(datadir, sim, k, append_step_to_name);
  end;
end

#! Create a callback function for writing a jld backup file
function write_jld_file_callback(datadir::AbstractString, stepout::Real,
                                 append_step_to_name::Bool = false)

  return @create_callback(stepout,
      dumpsim_jld(datadir, sim, k, append_step_to_name);
  );
end
#! Create a callback function for writing a backup file
function write_backup_file_callback(datadir::AbstractString)
  return (sim::Sim, k::Real) -> begin
    dumpsim(datadir, sim, k);
  end;
end

#! Create a callback function for writing a backup file
function write_backup_file_callback(datadir::AbstractString, stepout::Real)
  return @create_callback(stepout,
      dumpsim(datadir, sim, k);
  );
end

#! Create a callback function for writing to a delimited file
#!
#! \param pre Prefix of filename
#! \param stepout Number of steps in between writing
#! \param A Function for extracting a 2D array from the sim
#! \param delim Delimiter to separate values with
function write_datafile_callback(pre::AbstractString, stepout::Real, 
                                 A::LBXFunction; dir=".", delim=',')

  return @create_callback(stepout,
      writedlm(joinpath(dir, pre*"_step-$k.dsv"), A(sim), delim);
  );

end

#! Record a snapshot in a text file
#!
#! \param   fname     File name
#! \param   stepout   Number of steps in between writing
#! \param   A         Function for extracting data points
#! \param   delim     Delimiter for file
#! \return            Callback function
function take_snapshot_callback(fname::AbstractString, stepout::Real, 
                                A::LBXFunction; dir=".", delim::Char=',')
  if !isdir(dir); mkpath(dir); end
  fhandle = open(joinpath(dir, fname), "a");
  f = @create_callback(stepout, begin
                        write(fhandle, join(A(sim), delim));
                        write(fhandle, '\n');
                       end);
  # finalizer(f, () -> close(fhandle));
  return f;
end

include("sample.jl");
