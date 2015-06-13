const api_root = dirname(@__FILE__);

require(abspath(joinpath(api_root, "boundary.jl")));
require(abspath(joinpath(api_root, "collision.jl")));
require(abspath(joinpath(api_root, "convergence.jl")));
require(abspath(joinpath(api_root, "lattice.jl")));
require(abspath(joinpath(api_root, "lbxio.jl")));
require(abspath(joinpath(api_root, "multiscale.jl")));
require(abspath(joinpath(api_root, "simulate.jl")));

#import YAML
using PyCall
@pyimport yaml # more mature parser than YAML.jl

function parse_and_run(infile::String, args::Dict)

  if !isfile(infile)
    warn(infile * " not found. Please use valid path to input file.");
    return nothing;
  end

  const DEF_EXPR_ATTRS = {
    "col_f"         =>  { :store => true,  :array => false, :type => Function },
    "bcs"           =>  { :store => true,  :array => true,  :type => Function },
    "callbacks"     =>  { :store => true,  :array => true,  :type => Function },
    "finally"       =>  { :store => true,  :array => true,  :type => Function },
    "test_for_term" =>  { :store => true,  :array => false, :type => Function }
  };

  if args["verbose"]; info("parsing $infile from yaml..."); end
  ins = yaml.load(readall(infile));
  defs = Dict();

  if haskey(ins, "preamble")
    eval(parse(pop!(ins, "preamble")));
  end

  for (k, v) in ins

    if haskey(DEF_EXPR_ATTRS, k)

      if !DEF_EXPR_ATTRS[k][:store]
        eval(parse(v));
      else
        if DEF_EXPR_ATTRS[k][:array]
          n = length(v);
          # TODO: this syntax is soon deprecated
          defs[k] = Array(DEF_EXPR_ATTRS[k][:type], (n));
          for i=1:n
            defs[k][i] = eval(parse(v[i]));
          end
        else
          defs[k] = convert(DEF_EXPR_ATTRS[k][:type], eval(parse(v)));
        end
      end

    else
      if typeof(v) <: Dict
        if v["expr"]
          defs[k] = eval(parse(v["value"]));
        else
          defs[k] = v["value"];
        end
      else
        defs[k] = v;
      end

    end

  end

  const DEFAULTS = {
    "datadir" =>  (defs::Dict) -> begin; global datadir; return datadir; end,
    "rho_0"   =>  (defs::Dict) -> error("`rho_0` is a required parameter."),
    "nu"      =>  (defs::Dict) -> error("`nu` is a required parameter."),
    "dx"      =>  (defs::Dict) -> return 1.0,
    "dt"      =>  (defs::Dict) -> return 1.0,
    "ni"      =>  (defs::Dict) -> error("`ni` is a required parameter."),
    "nj"      =>  (defs::Dict) -> error("`nj` is a required parameter."),
    "nsteps"  =>  (defs::Dict) -> error("`nsteps` is a required parameter."),
    "col_f"   =>  (defs::Dict) -> error("`col_f` is a required parameter."),
    "sbounds" =>  (defs::Dict) -> [[1, defs["ni"], 1, defs["nj"]]],
    "cbounds" =>  (defs::Dict) -> [[1, defs["ni"], 1, defs["nj"]]],
    "bcs"     =>  (defs::Dict) -> Array(Function, 0),
    "callbacks" =>  (defs::Dict) -> Array(Function, 0),
    "finally" =>  (defs::Dict) -> Array(Function, 0)
  };

  if args["verbose"]; info("setting defaults."); end
  # set defaults
  for (k, f) in DEFAULTS
    if !haskey(defs, k)
      defs[k] = f(defs);
    end
  end

  # syntactic sugar for backing up simulation on program exit
  if haskey(defs, "backup_on_exit") && defs["backup_on_exit"]
    defs["finally"][end] = write_backup_file_callback(defs["datadir"], 1);
  end

  if args["verbose"]; println("$infile definitions:"); println(defs); end

  # if datadir does not exist, create it
  if !isdir(defs["datadir"])
    info(defs["datadir"] * " does not exist. Creating now...");
    mkpath(defs["datadir"]); # makes all directories in a given path
  end

  if args["resume"] && isfile(joinpath(defs["datadir"], "sim.bak"))
    k, sim = load_backup_file(joinpath(defs["datadir"], "sim.bak"));
  else
    # construct objects
    k = 0; # this is so every simulation can start from "k+1"
    lat = Lattice(defs["dx"], defs["dt"], defs["ni"], defs["nj"], defs["rho_0"]);
    msm = MultiscaleMap(defs["nu"], lat, defs["rho_0"]);
    sim = Sim(lat, msm);
  end

  if !args["debug"]

    try
      if !haskey(defs, "test_for_term")
        # this simulate should be more memory and computationally efficient
        nsim = simulate!(sim, defs["sbounds"], defs["col_f"], defs["cbounds"], 
          defs["bcs"], defs["nsteps"], defs["callbacks"], k);
      else
        nsim = simulate!(sim, defs["sbounds"], defs["col_f"], defs["cbounds"], 
          defs["bcs"], defs["nsteps"], defs["test_for_term"], defs["callbacks"],
          k);
      end

    catch e
      showerror(STDERR, e);
      println();
      Base.show_backtrace(STDERR, catch_backtrace()); # display callstack
      warn("$infile: not completed successfully.");

    finally
      for fin in defs["finally"]
        fin(sim, nsim);
      end

      println("$infile:\tSteps simulated: $nsim");

    end

  else # debugging on, run simulation without exception catching

    if !haskey(defs, "test_for_term")
      # this simulate should be more memory and computationally efficient
      nsim = simulate!(sim, defs["sbounds"], defs["col_f"], defs["cbounds"], 
        defs["bcs"], defs["nsteps"], defs["callbacks"], k);
    else
      nsim = simulate!(sim, defs["sbounds"], defs["col_f"], defs["cbounds"], 
        defs["bcs"], defs["nsteps"], defs["test_for_term"], defs["callbacks"],
        k);
    end
    
    println("$infile:\tSteps simulated: $nsim");

  end

end
