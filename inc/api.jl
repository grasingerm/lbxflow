# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

const __api_root__ = dirname(@__FILE__);

require(abspath(joinpath(__api_root__, "boundary.jl")));
require(abspath(joinpath(__api_root__, "col", "collision.jl")));
require(abspath(joinpath(__api_root__, "convergence.jl")));
require(abspath(joinpath(__api_root__, "lattice.jl")));
require(abspath(joinpath(__api_root__, "lbxio.jl")));
require(abspath(joinpath(__api_root__, "multiscale.jl")));
require(abspath(joinpath(__api_root__, "profile.jl")));
require(abspath(joinpath(__api_root__, "sim", "simulate.jl")));

using PyCall
@pyimport yaml

function parse_and_run(infile::AbstractString, args::Dict)

  if !isfile(infile)
    warn(infile * " not found. Please use valid path to input file.");
    return nothing;
  end

  const DEF_EXPR_ATTRS = Dict{AbstractString, Dict{Symbol, Any}}(
    "col_f"         =>  Dict{Symbol, Any}( :store => true,  :array => false,
                                           :type => Function ),
    "bcs"           =>  Dict{Symbol, Any}( :store => true,  :array => true,  
                                           :type => Function ),
    "callbacks"     =>  Dict{Symbol, Any}( :store => true,  :array => true,  
                                           :type => Function ),
    "finally"       =>  Dict{Symbol, Any}( :store => true,  :array => true,  
                                           :type => Function ),
    "test_for_term" =>  Dict{Symbol, Any}( :store => true,  :array => false, 
                                           :type => Function )
  );

  if args["verbose"]; info("parsing $infile from yaml..."); end
  ins = yaml.load(readall(infile));

  if haskey(ins, "version") && eval(parse("v\"$(ins["version"])\"")) > LBX_VERSION
    warn("$infile recommends v$(ins["version"]), consider updating.");
  end

  if haskey(ins, "preamble")
    if args["verbose"]; info("evaluating preamble..."); end;
    pre = pop!(ins, "preamble");
    if start(search(pre, "#")) != 0
      warn("`#` character in the preamble can lead to undefined behavior");
    end
    eval(parse(pre));
  end

  defs = Dict();
  for (k, v) in ins

    if haskey(DEF_EXPR_ATTRS, k)

      if !DEF_EXPR_ATTRS[k][:store]
        if args["verbose"]; info("evalling... $v"); end
        eval(parse(v));
      else
        if DEF_EXPR_ATTRS[k][:array]
          n = length(v);
          # TODO: this syntax is soon deprecated; since when??
          defs[k] = Array(DEF_EXPR_ATTRS[k][:type], (n));
          for i=1:n
            if args["verbose"]; info("evalling... $(v[i])"); end
            defs[k][i] = eval(parse(v[i]));
          end
        else
          if args["verbose"]; info("evalling... $v"); end
          defs[k] = convert(DEF_EXPR_ATTRS[k][:type], eval(parse(v)));
        end
      end

    else
      if typeof(v) <: Dict
        if v["expr"]
          if args["verbose"]; info("evalling... $v"); end
          defs[k] = eval(parse(v["value"]));
        else
          defs[k] = v["value"];
        end
      else
        defs[k] = v;
      end

    end

  end

  const DEF_DEFAULTS = Dict{AbstractString, Any}(
    "simtype" =>  (defs::Dict) -> begin; return "default"; end,
    "rhog"    =>  (defs::Dict) -> begin; return 1.0; end,
    "datadir" =>  (defs::Dict) -> begin; global datadir; return datadir; end,
    "rho_0"   =>  (defs::Dict) -> begin; error("`rho_0` is a usingd parameter."); end,
    "nu"      =>  (defs::Dict) -> begin; error("`nu` is a usingd parameter."); end,
    "dx"      =>  (defs::Dict) -> begin; return 1.0; end,
    "dt"      =>  (defs::Dict) -> begin; return 1.0; end,
    "ni"      =>  (defs::Dict) -> begin; error("`ni` is a usingd parameter."); end,
    "nj"      =>  (defs::Dict) -> begin; error("`nj` is a usingd parameter."); end,
    "nsteps"  =>  (defs::Dict) -> begin; error("`nsteps` is a usingd parameter."); end,
    "col_f"   =>  (defs::Dict) -> begin; error("`col_f` is a usingd parameter."); end,
    "sbounds" =>  (defs::Dict) -> begin; [1 defs["ni"] 1 defs["nj"];]'; end,
    "cbounds" =>  (defs::Dict) -> begin; [1 defs["ni"] 1 defs["nj"];]'; end,
    "bcs"     =>  (defs::Dict) -> begin; Array(Function, 0); end,
    "callbacks" =>  (defs::Dict) -> begin; Array(Function, 0); end,
    "finally" =>  (defs::Dict) -> begin; Array(Function, 0); end
  );

  if args["verbose"]; info("setting defaults."); end
  # set defaults
  for (k, f) in DEF_DEFAULTS
    if !haskey(defs, k)
      defs[k] = f(defs);
    end
  end

  # syntactic sugar for backing up simulation on program exit
  if haskey(defs, "backup_on_exit") && defs["backup_on_exit"]
    defs["finally"][end] = write_jld_file_callback(defs["datadir"], 1);
  end

  if args["verbose"]; println("$infile definitions:"); println(defs); end

  # if datadir does not exist, create it
  if !isdir(defs["datadir"])
    info(defs["datadir"] * " does not exist. Creating now...");
    mkpath(defs["datadir"]); # makes all directories in a given path
  end

  # recursively remove data from previous runs
  if args["clean"]
    for f in readdir(defs["datadir"]); rrm(joinpath(defs["datadir"], f)); end
    if args["noexe"]; rm(defs["datadir"]); end
  end

  # check for `noexe` switch
  if args["noexe"]
    info("`noexe` switch was passed. $infile will not be simulated.");
    return nothing;
  end

  is_init = false;
  if args["resume"]
    k, sim = load_latest_backup(defs["datadir"]);
    if k == 0 || sim == nothing
      if args["verbose"]; info("No backup files found."); end
    else
      if args["verbose"]; info("Loaded previous simulation data from step $k."); end
      is_init = true;
    end
  end

  if !is_init
    # construct objects
    k = 0; # this is so every simulation can start from "k+1"
    lat = LatticeD2Q9(defs["dx"], defs["dt"], defs["ni"], defs["nj"], defs["rho_0"]);
    msm = MultiscaleMap(defs["nu"], lat, defs["rho_0"]);
    if defs["simtype"] == "default"; sim = Sim(lat, msm)
    elseif defs["simtype"] == "free_surface"
      if haskey(defs, "states")
        sim = FreeSurfSim(lat, msm, Tracker(msm, defs["states"]), defs["rhog"]);
      else
        error("No `states` matrix provided. Cannot initialize free surface flow");
      end
    else
      error("`simtype` $(defs["simtype"]) is not understood");
    end
  end
 
  if args["profile"]
    Profile.clear();
    Profile.init(delay=args["profile-delay"]);
  end

  # Do not simulate, only post process
  if args["post-process"]
    defs["test_for_term"] = (msm, prev_msm) -> true;
  end

  nsim = 0;
  if !args["debug"]

    try
      tic();

      if !haskey(defs, "test_for_term")
        # this simulate should be more memory and computationally efficient
        @profif(args["profile"], begin; nsim = simulate!(sim,
                defs["sbounds"], defs["col_f"], defs["cbounds"], defs["bcs"],
                defs["nsteps"], defs["callbacks"], k); end);
      else
        @profif(args["profile"], begin; nsim = simulate!(sim,
                defs["sbounds"], defs["col_f"], defs["cbounds"], defs["bcs"],
                defs["nsteps"], defs["test_for_term"], defs["callbacks"], k); 
                end);
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
      toc();

      if args["profile"]
        Profile.print(args["profile-io"], cols=args["profile-cols"]);
      end

      if args["profile-file"] != nothing; 
        close(args["profile-io"]);
        bt, lidict = Profile.retrieve();
        JLD.jldopen(args["profile-file"]*".jld", "w") do file
          write(file, "bt", bt);
          write(file, "lidict", lidict);
        end
      elseif args["profile"] && args["profile-view"]
        ProfileView.view();
        println("Press enter to continue...");
        readline(STDIN);
      end
    end
  else
    error("Debugging mode is not yet available");
  end

  #=
  # TODO: is `debugging` mode actually useful?
  else # debugging on, run simulation without exception catching
    display_help = true
    println("======================");
    println("DEBUGGING MODE... REPL");
    println("======================");
    warn("debugging mode is still experimental and not yet working smoothly");
    warn("... having trouble with variable and function scoping");
    simulate!(nsteps) = simulate!(sim, defs["sbounds"], defs["col_f"],
                                  defs["cbounds"], defs["bcs"], nsteps,
                                  defs["callbacks"], k);
    while true
      # bring in everything from the global scope
      # TODO: find a better solution...
      global simulate!
      global display_help
      global sim
      global msm
      if display_help
        println("Binding input file definitions to `simulate!` function...");
        println("In interactive debugging mode...");
        println("`simulate!(<nsteps>)` will run the simulation with input ",
                "file defintions for <nsteps> number of time steps.");
        println("`exit()` will end the REPL");
        println("to suppress this message, `display_help = false`");
      end
      print(">> ");
      command = readline(STDIN);
      try
        # bring in everything from the global scope
        # TODO: find a better solution...
        global simulate!
        global display_help
        global sim
        global msm
        eval(parse(command));
      catch e
        showerror(STDERR, e);
        Base.show_backtrace(STDERR, catch_backtrace()); # display callstack
        println();
      end
    end

    if args["profile"]; Profile.print(args["profile-io"], cols=args["profile-cols"]); end
      if args["profile-file"] != nothing; 
        close(args["profile-io"]);
        bt, lidict = Profile.retrieve();
        JLD.jldopen(args["profile-file"]*".jld", "w") do file
          write(file, "bt", bt);
          write(file, "lidict", lidict);
        end
      elseif args["profile"] && args["profile-view"]
        ProfileView.view();
        println("Press enter to continue...");
        readline(STDIN);
      end
    end
    =#
end
