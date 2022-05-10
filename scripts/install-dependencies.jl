import Pkg;
Pkg.update();

for pkgname in (
  "ArgParse",
  "PyCall",
  "PyPlot",
  "YAML",
  "HDF5",
  "JLD",
  "ProfileView",
  "Roots"
  )

  println("Checking for package $pkgname");
  Pkg.add(pkgname);
end

  #=
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
=#
