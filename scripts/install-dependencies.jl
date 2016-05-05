function confirm(prompt::AbstractString="Would you like to continue, (Y/N)?: ")
  print(prompt);
  answer = chomp(lowercase(readline(STDIN)));

  while answer != "y" && answer != "n"
    println("Do not understand: $answer");
    print(prompt);
    answer = chomp(lowercase(readline(STDIN)));
  end

  return answer == "y";
end

if confirm("Would you like to update Julia packages (Y/N)?: ")
  Pkg.update();

  for pkgname in (
    "ArgParse",
    "PyCall",
    "PyPlot",
    "YAML",
    "HDF5",
    "JLD",
    "ArrayViews",
    "ProfileView",
    "Roots",
    "FastAnonymous"
    )

    println("Checking for package $pkgname");
    Pkg.add(pkgname);
  end

  try
    Pkg.clone("https://github.com/grasingerm/LazyViews.jl.git");
  catch e
    warn("unable to clone LazyViews.jl");
  end
end

if confirm("Would you like to add LBXFlow.jl module to your load path (Y/N)? (NOT recommended): ")
  function add_lbxflow_to_load_path(fp::AbstractString)
    const defhome = homedir();
    println("What is your home directory? (Enter for default: $defhome) ");
    home = readline(STDIN);

    home = (home == "\n") ? defhome : home;
    jrc = open(joinpath(home, ".juliarc.jl"), "a");
    write(jrc, "
  # This line was added by the lbxflow installation in order that LBXFlow
  # module might be found in the LOAD_PATH
  push!(LOAD_PATH, \"$fp\");
  ");
    close(jrc);
  end

  if isdir("inc") && isfile(joinpath("inc", "LBXFlow.jl"))
    const fp = abspath("inc");
    if !(fp in LOAD_PATH || joinpath(fp, "LBXFlow.jl") in LOAD_PATH)
      info("$fp is not contained in LOAD_PATH, adding now...");
      add_lbxflow_to_load_path(fp);
    end
  elseif (isdir(joinpath("..", "inc")) && 
          isfile(joinpath("..", "inc", "LBXFlow.jl")))
    const fp = abspath("..", "inc");
    if !(fp in LOAD_PATH || joinpath(fp, "LBXFlow.jl") in LOAD_PATH)
      info("$fp is not contained in LOAD_PATH, adding now...");
      add_lbxflow_to_load_path(fp);
    end
  else
    warn("Cannot find directory where LBXFlow module is contained.");
    error("Try rerunning the install script from lbxflow root!");
  end
end

if confirm("Would you like to check python dependencies (Y/N)?: ")
  warn("For python dependencies, pip must be installed!");
  try
    println("Attempting to install PyYaml...");
    run(`sudo pip install pyyaml`);
    println("Attempting to install numpy...");
    run(`sudo pip install numpy`);
    println("Attempting to install matplotlib...");
    run(`sudo pip install matplotlib`);
  catch e
    warn("unable to install python dependencies");
  end
end
