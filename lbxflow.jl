const root = dirname(@__FILE__);

println(readall(abspath(joinpath(root, "banner.txt"))));

require(abspath(joinpath(root, "inc", "api.jl")));
require(abspath(joinpath(root, "inc", "boundary.jl")));
require(abspath(joinpath(root, "inc", "collision.jl")));
require(abspath(joinpath(root, "inc", "convergence.jl")));
require(abspath(joinpath(root, "inc", "lattice.jl")));
require(abspath(joinpath(root, "inc", "lbxio.jl")));
require(abspath(joinpath(root, "inc", "multiscale.jl")));
require(abspath(joinpath(root, "inc", "simulate.jl")));

function main(inputfile)
  println("Load simulation defintions from inputfile $inputfile");
  const def = load_sim_definitions(inputfile);

  println("Definitions: ", def);

  lat = Lattice(def["dx"], def["dt"], def["ni"], def["nj"], def["rhoo"]);
  msm = MultiscaleMap(def["nu"], lat, def["rhoo"]);

  # if data directory does not exist, create it
  if !isdir(def["datadir"])
    println(def["datadir"], " does not exist. Creating now...");
    mkdir(def["datadir"]);
  end

  if in("presim", keys(def))
    def["presim"](msm);
  end

  println("Starting simulation...");
  println();

  if in("test_for_term", keys(def))
    const n = simulate!(lat, msm, def["col_f"], def["bcs"], def["nsteps"],
      def["test_for_term"], def["callbacks"]);
  else
    const n = simulate!(lat, msm, def["col_f"], def["bcs"], def["nsteps"],
      def["callbacks"]);
  end

  if in("postsim", keys(def))
    def["postsim"](msm);
  end

  println("\nSteps simulated: $n");

end

# parsing arguments and user interface
if length(ARGS) != 1
  println("usage: julia lbxflow.jl inputfile.js");
  exit(1);
end

if !isfile(ARGS[1])
  println(ARGS[1], " not found. Please use valid path to input file.");
end

main(ARGS[1]);
