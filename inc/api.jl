# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# TODO this entire mess needs rewritten... If you're reading this it's too late

using PyCall;
@pyimport yaml;

function parse_and_run(infile::String, args::Dict)

  if !isfile(infile)
    warn(infile * " not found. Please use valid path to input file.");
    return nothing;
  end

  DEF_EXPR_ATTRS = Dict{String, Dict{Symbol, Any}}(
    "col_f"         =>  Dict{Symbol, Any}( :store => true,  :array => false,
                                           :type => LBXFunction ),
    "bcs"           =>  Dict{Symbol, Any}( :store => true,  :array => true,  
                                           :type => LBXFunction ),
    "callbacks"     =>  Dict{Symbol, Any}( :store => true,  :array => true,  
                                           :type => LBXFunction ),
    "finally"       =>  Dict{Symbol, Any}( :store => true,  :array => true,  
                                           :type => LBXFunction ),
    "test_for_term" =>  Dict{Symbol, Any}( :store => true,  :array => false, 
                                           :type => LBXFunction )
  );

  if args["verbose"]; info("parsing $infile from yaml..."); end
  ins = yaml.load(readall(infile));

  # Check version is consistent with input file
  if haskey(ins, "version")

    @assert(haskey(args, "LBX_VERSION"), "LBX_VERSION should be passed from " *
                                         "lbxflow.jl");
    ins["version"] = eval(parse("v\"$(ins["version"])\""));
    @assert(ins["version"].major == args["LBX_VERSION"].major,
            "Major version specified in $infile, $(ins["version"]), does not " *
            "match major version of LBXFlow, $(args["LBX_VERSION"])."); 
    if ins["version"] > args["LBX_VERSION"]
      warn("$infile recommends v$(ins["version"]), consider updating.");
    end

  else

    warn("$infile does not contain versioning information.");

  end

  if haskey(ins, "preamble")
    if args["verbose"]; info("evaluating preamble..."); end;
    pre = pop!(ins, "preamble");
    if start(search(pre, "#")) != 0
      warn("`#` character in the preamble can lead to undefined behavior");
    end
    eval(parse(pre));
  end

  global defs = Dict();
  for (k, v) in ins

    if haskey(DEF_EXPR_ATTRS, k)

      if !DEF_EXPR_ATTRS[k][:store]
        if args["verbose"]; info("evalling... $v"); end
        eval(parse(v));
      else
        if DEF_EXPR_ATTRS[k][:array]
          n = length(v);
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

  DEF_DEFAULTS = Dict{String, Any}(
    "simtype" =>  (defs::Dict) -> begin; return "default"; end,
    "rho_g"   =>  (defs::Dict) -> begin; return 1.0; end,
    "datadir" =>  (defs::Dict) -> begin; global datadir; return datadir; end,
    "dx"      =>  (defs::Dict) -> begin; return 1.0; end,
    "dt"      =>  (defs::Dict) -> begin; return 1.0; end,
    "ni"      =>  (defs::Dict) -> begin; error("`ni` is a required parameter."); end,
    "nj"      =>  (defs::Dict) -> begin; error("`nj` is a required parameter."); end,
    "nsteps"  =>  (defs::Dict) -> begin; error("`nsteps` is a required parameter."); end,
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

  # syntactic sugar for adding obstacles
  if haskey(defs, "obstacles")
    defs["active_cells"]  =   fill(true, defs["ni"], defs["nj"]);
    for obstacle_array in defs["obstacles"]
      obs_type      =   parse(obstacle_array["type"]);
      obs_locs      =   eval(parse(obstacle_array["coords"]));

      @assert(size(obs_locs, 1) == 4, "Obstacle coordinates should be "      *
              "organized in columns, i.e. coords[:, ...] = i_min, i_max, "   *
              "j_min, j_max. Obstacle coordinate matrix should have exactly "*
              "4 rows. $(size(obs_locs, 1)) != 4");

      for j = 1:size(obs_locs, 2)
        add_obstacle!(defs["active_cells"], defs["bcs"], obs_locs[1, j], 
                      obs_locs[2, j], obs_locs[3, j], obs_locs[4, j], obs_type); 
      end
    end
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
    k = 0.0; # this is so every simulation can start from "k+1"
    simtypes = map(s -> strip(s), split(defs["simtype"], ','));
    if "m2phase" in simtypes
      for key in (["Ar", "Ab", "αr", "αb", "β", "nu_r", "nu_b", "rho_0r", 
                   "rho_0b"])
        @assert(haskey(defs, key), "`$key` is a required input");
      end
      fill_r = (haskey(defs, "fill_r")) ? defs["fill_r"] : (0.0, 1.0, 0.0, 1.0);
      fill_b = (haskey(defs, "fill_b")) ? defs["fill_b"] : (0.0, 1.0, 0.0, 1.0);
      isim = M2PhaseSim(defs["nu_r"], defs["nu_b"], defs["rho_0r"], 
                        defs["rho_0b"], defs["ni"], defs["nj"], defs["Ar"], 
                        defs["Ab"], defs["αr"], defs["αb"], defs["β"];
                        fill_r=fill_r, fill_b=fill_b);
    else
      @assert(haskey(defs, "nu"), "`nu` is a required input");
      @assert(haskey(defs, "rho_0"), "`rho_0` is a required input");
      lat = LatticeD2Q9(defs["dx"], defs["dt"], defs["ni"], defs["nj"], defs["rho_0"]);
      msm = MultiscaleMap(defs["nu"], lat, defs["rho_0"]);
      if "default" in simtypes
        @assert(!("free_surface" in simtypes), "`simtype` cannot be both `default` and `free_surface`");
        @assert(!("m2phase" in simtypes), "`simtype` cannot be both `default` and `m2phase`");
        @assert(haskey(defs, "col_f"), "`col_f` must be provided in input file.");
        isim = Sim(lat, msm);
      elseif "free_surface" in simtypes
        @assert(!("m2phase" in simtypes), "`simtype` cannot be both `free_surface` and `m2phase`");
        @assert(haskey(defs, "col_f"), "`col_f` must be provided in input file.");
        if haskey(defs, "states")
          if haskey(defs, "fill_x") || haskey(defs, "fill_y")
            warn("'States' matrix was already provided. 'fill_\$D' variables " *
                 "will be ignored");
          end
          isim = FreeSurfSim(lat, msm, Tracker(msm, defs["states"]), defs["rho_g"]);
        elseif haskey(defs, "fill_x") && haskey(defs, "fill_y")
          isim = FreeSurfSim(lat, msm, defs["rho_0"], defs["rho_g"], 
                            defs["fill_x"], defs["fill_y"]);
        else
          error("No `states` matrix provided. Cannot initialize free surface flow");
        end
      else
        error("`simtype` $(defs["simtype"]) is not understood");
      end
    end
    if "adaptive" in simtypes # adaptive time stepping?
      warn("Adaptive time stepping is highly experimental. It has not yet been validated");
      incr = (!haskey(defs, "incr") || defs["incr"]) ? true : false;
      decr = (!haskey(defs, "decr") || defs["decr"]) ? true : false;
      relax = (haskey(defs, "relax")) ? defs["relax"] : 1.0;
      if haskey(defs, "xi")
        sim = AdaptiveTimeStepSim(isim, defs["xi"]; incr=incr, decr=decr, relax=relax);
      else
        sim = AdaptiveTimeStepSim(isim; incr=incr, decr=decr, relax=relax);
      end
    else
      sim = isim;
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

  mass_matrix = (if (haskey(defs, "simtype") 
                     && defs["simtype"] == "free_surface")
                    sim.tracker.M
                  else 
                    zeros(1);
                  end);

  @_checkdebug_mass_cons("simulation", mass_matrix, try
    tic();

    if !haskey(defs, "obstacles")
      if !haskey(defs, "test_for_term")
        # this simulate should be more memory and computationally efficient
        @profif(args["profile"],
          begin;
            global nsim;
            nsim = simulate!(sim, defs["sbounds"], defs["col_f"], defs["cbounds"],
                             defs["bcs"], defs["nsteps"], defs["callbacks"], k); 
          end;
        );
      else
        if haskey(defs, "test_for_term_steps")
          @assert(defs["test_for_term_steps"] > 1, 
                  "'test_for_term_steps', i.e. the number of steps to average " *
                  "over when checking for steady-state conditions, should be "  *
                  "greater than 1 (or not specified).");
          @profif(args["profile"],
            begin;
              global nsim;
              nsim = simulate!(sim, defs["sbounds"], defs["col_f"], defs["cbounds"],
                               defs["bcs"], defs["nsteps"], defs["test_for_term"], 
                               defs["test_for_term_steps"], defs["callbacks"], k); 
            end;
          );
        else
          @profif(args["profile"],
            begin;
              global nsim;
              nsim = simulate!(sim, defs["sbounds"], defs["col_f"], defs["cbounds"],
                               defs["bcs"], defs["nsteps"], defs["test_for_term"], 
                               defs["callbacks"], k); 
            end;
          );
        end
      end
    else
      if !haskey(defs, "test_for_term")
        # this simulate should be more memory and computationally efficient
        @profif(args["profile"],
          begin;
            global nsim;
            nsim = simulate!(sim, defs["col_f"], defs["active_cells"],
                             defs["bcs"], defs["nsteps"], defs["callbacks"], k); 
          end;
        );
      else
        if haskey(defs, "test_for_term_steps")
          @assert(defs["test_for_term_steps"] > 1, 
                  "'test_for_term_steps', i.e. the number of steps to average " *
                  "over when checking for steady-state conditions, should be "  *
                  "greater than 1 (or not specified).");
          @profif(args["profile"],
            begin;
              global nsim;
              nsim = simulate!(sim, defs["col_f"], defs["active_cells"],
                               defs["bcs"], defs["nsteps"], defs["test_for_term"], 
                               defs["test_for_term_steps"], defs["callbacks"], k); 
            end;
          );
        else
          @profif(args["profile"],
            begin;
              global nsim;
              nsim = simulate!(sim, defs["col_f"], defs["active_cells"],
                               defs["bcs"], defs["nsteps"], defs["test_for_term"], 
                               defs["callbacks"], k); 
            end;
          );
        end
      end
    end

  catch e
    @_report_and_exit(e, 0.0);

  finally
    for fin in defs["finally"]
      fin(sim, nsim);
    end

    if args["verbose"]
      println("$infile:\tSteps simulated: $nsim");
      toc();
    end

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
  end, 1e-12);

end
