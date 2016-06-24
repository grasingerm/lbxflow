println("Testing srt_col_f...");

require("collision.jl");
require("defaults.jl");
require("lattice.jl");
require("multiscale.jl");
require("simulate.jl");

lat = Lattice(1.0, 1.0, 10, 10, 9.0);
msm = MultiscaleMap(0.1, lat);

# initialize lattice
for k=1:length(lat.w)
  lat.f[5,5,k] = 10.0 * lat.w[k];
end
map_to_macro!(lat, msm);

println("initial conditions");
println("rho=")
println(msm.rho);
println("u=");
println(msm.u[:,:,1]);
println("v=");
println(msm.u[:,:,2]);

temp_f = zeros(lat.f);

for i=1:4
  stream!(lat, temp_f);
  srt_col_f!(lat, msm);
  map_to_macro!(lat, msm);
  println("step i=", i);
  println("rho=")
  println(msm.rho);
  println("u=");
  println(msm.u[:,:,1]);
  println("v=");
  println(msm.u[:,:,2]);
end
