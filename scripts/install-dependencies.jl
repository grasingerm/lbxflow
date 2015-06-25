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

Pkg.clone("https://github.com/grasingerm/List.git");
