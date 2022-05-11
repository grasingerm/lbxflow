# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ccall(:jl_exit_on_sigint, Void, (Cint,), 0); # Allows Ctrl+C to be caught
const LBX_VERSION = v"1.1.0";

# load dependencies
using ArgParse;
using Logging;
using DelimitedFiles;

s = ArgParseSettings();
@add_arg_table s begin
  "--file", "-f"
    help = "path to input file"
  "--dir", "-d"
    help = "directory with input file(s)"
  "--recursive", "-r"
    help = "recursively search directories for input file(s)"
    action = :store_true
  "--ext", "-x"
    help = "file extension of input file(s) (for directory search)"
    default = "yaml"
  "--verbose", "-v"
    help = "verbosity level: 0-nothing, 1-errors, 2-warnings, 3-@info"
    arg_type = Int
    default = 3
  "--clean", "-c"
    help = "clean out `datadir` defined in input file(s)"
    action = :store_true
  "--resume", "-R"
    help = "search for backup files and use if applicable"
    action = :store_true
  "--post-process", "-G"
    help = "ONLY post process, don't simulate"
    action = :store_true
  "--ndebug", "-N"
    help = "turn off debugging macros"
    action = :store_true
  "--ndebug-mass-cons"
    help = "turn off conservation of mass debugging macros"
    action = :store_true
  "--profile", "-P"
    help = "profile `simulate` function"
    action = :store_true
  "--profile-file", "-p"
    help = "file for profiler to print to"
  "--profile-view", "-V"
    help = "load and use the ProfileView package"
    action = :store_true
  "--profile-delay", "-D"
    help = "time delay between calls to sampler"
    arg_type = Real
    default = 0.001
  "--profile-cols"
    help = "number of columns in Profile.print"
    arg_type = Int
    default = 40
  "--version"
    help = "display @information about the program"
    action = :store_true
end

pa = parse_args(s);
pa["LBX_VERSION"] = LBX_VERSION; # add version to arguments

if pa["verbose"] == 3
  global_logger(ConsoleLogger(stderr, Logging.Info));
elseif pa["verbose"] == 2
  global_logger(ConsoleLogger(stderr, Logging.Warn));
elseif pa["verbose"] == 1
  global_logger(ConsoleLogger(stderr, Logging.Error));
else
  global_logger(Logging.NullLogger());
end

pa["verbose"] = pa["verbose"] >= 2;

if pa["version"]
  println();
  println(read(abspath(joinpath(dirname(@__FILE__), "banner.txt")), String));
  println("version:  $LBX_VERSION");
  println();
  exit(0);
end

# organize profiling args
if pa["profile-file"] != nothing
  pa["profile"] = true # e have an output file location, we should be profiling
  pa["profile-io"] = open((pa["profile-file"]*".txt"), "w");
elseif pa["profile"]
  pa["profile-io"] = STDOUT;
end

if pa["profile-view"]; 
  using ProfileView; 
end

# Load lattice Boltzmann method simulation module
push!(LOAD_PATH, "inc");
import LBXFlow;

if pa["ndebug"];            LBXFlow.turn_off_debugging();           end;
if pa["ndebug-mass-cons"];  LBXFlow.turn_off_mass_cons_debugging(); end;

# run input file
if pa["file"] != nothing
  LBXFlow.parse_and_run(pa["file"], pa);
end

# recursively search directory for input files
function recursively_search_for_input_files(dir::AbstractString, 
                                            ext::AbstractString)
  files = String[];
  function recursively_add_input_files!(files::Array{String}, 
                                        dir::AbstractString, 
                                        ext::AbstractString)
    for f in readdir(dir)
      if isdir(joinpath(dir, f)); recursively_add_input_files!(files, joinpath(dir, f), ext); end
      if isfile(joinpath(dir, f)) && endswith(f, ext); push!(files, joinpath(dir, f)); end
    end
  end

  recursively_add_input_files!(files, dir, ext);
  return files;
end

# run input files from directory in parallel
if pa["dir"] != nothing

  if pa["recursive"]
    files = recursively_search_for_input_files(pa["dir"], pa["ext"]);
  else
    files = filter((x)->endswith(x, ".$(pa["ext"])"), readdir(pa["dir"]));
    map!((x)->joinpath(pa["dir"],x), files);
  end
  
  for file in files
    LBXFlow.parse_and_run(file, pa);
  end

end
