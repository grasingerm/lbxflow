Pkg.update();

for pkgname in (
  "PyCall",
  "PyPlot",
  "YAML",
  "HDF5",
  "ArrayViews"
  )

  Pkg.add(pkgname);
end

Pkg.clone("https://github.com/grasingerm/List.git");
