ccall(:jl_exit_on_sigint, Void, (Cint,), 0); # Allows Ctrl+C to be caught
const VERSION = 0.2;
const root = dirname(@__FILE__);

println();
println(readall(abspath(joinpath(root, "banner.txt"))));
println("version: $VERSION");
println();

# load dependencies
require(abspath(joinpath(root, "inc", "api.jl")));
require("argparse");
using ArgParse;

s = ArgParseSettings();
@add_arg_table s begin
  "--file", "-f"
    help = "path to input file"
    arg_type = String
  "--dir", "-d"
    help = "directory with input files"
    arg_type = String
  "--ext", "-x"
    help = "file extension of input files (for directory search)"
    default = "lbx"
  "--verbose", "-v"
    help = "print extra information about execution status"
    action = :store_true
  "--resume"
    help = "search for backup files and use if applicable"
    action = :store_true
  "--log"
    help = "create a logfile"
    action = :store_true
  "--debug"
    help = "turns on debugging mode"
    action = :store_true
end

pa = parse_args(s);

if pa["file"] != nothing
  parse_and_run(pa["file"], pa);
end

if pa["dir"] != nothing
  @parallel for file in readdir(pa["dir"])
    if endswith(file, "." * pa["ext"])
      parse_and_run(file, pa);
    end
  end
end
