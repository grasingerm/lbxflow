println("Testing streaming...");

require("collision.jl");
require("defaults.jl");
require("lattice.jl");
require("multiscale.jl");
require("simulate.jl");

lat = Lattice(1.0, 1.0, 10, 10);
msm = MultiscaleMap(0.02, lat);

for k=1:length(lat.w)
  lat.f[5,5,k] = 1.0;
end

map_to_macro!(lat, msm);
println("initial conditions\n", msm.rho, "\n");

temp_f = zeros(lat.f);

for i=1:3
  stream!(lat, temp_f);
  map_to_macro!(lat, msm);
  println("step i=", i);
  println(msm.rho, "\n");
end
