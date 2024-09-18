import Pkg;
Pkg.update();

for pkgname in (
  "ArgParse",
  "PythonCall",
  "PythonPlot",
  "YAML",
  "HDF5",
  "JLD",
  "ProfileView",
  "Roots"
  )

  println("Checking for package $pkgname");
  Pkg.add(pkgname);
end
