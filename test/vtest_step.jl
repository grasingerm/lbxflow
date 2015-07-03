require("boundary.jl");
require("collision.jl");
require("defaults.jl")
require("lattice.jl");
require("lbxio.jl");
require("multiscale.jl");
require("simulate.jl");

# get default values from defaults.jl
def = lbx_defaults();

lat = Lattice(1.0, 1.0, 20, 20);
msm = MultiscaleMap(0.02, lat, 5.0);

#! Zero y-component of velocity at outlet
zero_v_at_outlet!(msm::MultiscaleMap, k::Int) = begin
  ni, nj = size(msm.u) 
  for j=1:nj
    msm.u[ni,j,2] = 0.0;
  end
end

temp_f = copy(lat.f);

for i=1:3
  srt_col_f!(lat, msm);
  println("f\n", lat.f);

  stream!(lat, temp_f);
  println(lat.f);

  for bc! in def["bcs"]
    bc!(lat);
  end
  println("f\n", lat.f);

  map_to_macro!(lat, msm);
  println("f\n", lat.f);
  println("rho\n", msm.rho);
  println("u\n", msm.u);
end

println("tau, ", msm.tau);
omega = 1.0/(3.0*def["nu"]+0.5);
println("omega, ", omega);
println("1/omega, ", 1/omega);
