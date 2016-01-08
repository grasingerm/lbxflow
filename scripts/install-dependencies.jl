Pkg.update();

for pkgname in (
  "ArgParse",
  "PyCall",
  "PyPlot",
  "YAML",
  "HDF5",
  "JLD",
  "ArrayViews",
  "ProfileView"
  )

  Pkg.add(pkgname);
end

try
  Pkg.clone("https://github.com/grasingerm/List.jl.git");
catch e
  warn("unable to clone List.jl");
end

warn("For python dependencies, pip must be installed!");
try
  println("Attempting to install PyYaml...");
  run(`pip install pyyaml`);
  println("Attempting to install numpy...");
  run(`pip install numpy`);
  println("Attempting to install matplotlib...");
  run(`pip install matplotlib`);
catch e
  warn("unable to install python dependencies");
end
