# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

ccall(:jl_exit_on_sigint, Void, (Cint,), 0); # Allows Ctrl+C to be caught
const LBX_VERSION = v"0.2.4";
const root = dirname(@__FILE__);

# load dependencies
require(abspath(joinpath(root, "inc", "api.jl")));
using ArgParse;

const term_rows, term_cols = Base.tty_size();

s = ArgParseSettings();
@add_arg_table s begin
  "--file", "-f"
    help = "path to input file"
    arg_type = AbstractString
  "--dir", "-d"
    help = "directory with input file(s)"
    arg_type = AbstractString
  "--recursive", "-r"
    help = "recursively search directories for input file(s)"
    action = :store_true
  "--ext", "-x"
    help = "file extension of input file(s) (for directory search)"
    default = "yaml"
  "--verbose", "-v"
    help = "print extra information about execution status"
    action = :store_true
  "--clean", "-c"
    help = "clean out `datadir` defined in input file(s)"
    action = :store_true
  "--resume", "-R"
    help = "search for backup files and use if applicable"
    action = :store_true
  "--noexe", "-X"
    help = "do not execute input file(s)"
    action = :store_true
  "--post-process", "-G"
    help = "ONLY post process, don't simulate"
    action = :store_true
  "--profile", "-P"
    help = "profile `simulate` function"
    action = :store_true
  "--profile-file", "-p"
    help = "file for profiler to print to"
    arg_type = AbstractString
  "--profile-view", "-V"
    help = "load and use the ProfileView package"
    action = :store_true
  "--profile-delay", "-D"
    help = "time delay between calls to sampler"
    arg_type = Number
    default = 0.001
  "--profile-cols"
    help = "number of columns in Profile.print"
    arg_type = Int
    default = term_cols
  "--version"
    help = "display information about the program"
    action = :store_true
  "--debug"
    help = "turns on debugging mode"
    action = :store_true
end

pa = parse_args(s);

if pa["version"]
  println();
  println(readall(abspath(joinpath(root, "banner.txt"))));
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

if pa["profile-view"]; using ProfileView; end

# run input file
if pa["file"] != nothing
  parse_and_run(pa["file"], pa);
end

# recursively search directory for input files
function recursively_search_for_input_files(dir::AbstractString, 
                                            ext::AbstractString)
  files = Array(ASCIIString, 0);
  function recursively_add_input_files!(files::Array{ASCIIAbstractString}, 
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
    parse_and_run(file, pa);
  end

end
