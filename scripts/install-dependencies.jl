Pkg.update();

for pkgname in (
  "PyCall",
  "PyPlot",
  "YAML",
  "HDF5",
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

println("Attempting to install PyYaml...");
run(`pip install pyyaml`);
println("Attempting to install numpy...");
run(`pip install numpy`);
println("Attempting to install matplotlib...");
run(`pip install matplotlib`);
