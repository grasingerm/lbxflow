println("Flow through a 2D channel");

require("boundary.jl");
require("collision.jl");
require("defaults.jl");
require("lattice.jl");
require("lbxio.jl");
require("multiscale.jl");
require("simulate.jl");

# get default values from defaults.jl
def = lbx_defaults();

lat = Lattice(def["dx"], def["dt"], def["ni"], def["nj"]);
msm = MultiscaleMap(def["nu"], lat, def["rhoo"]);

#! Initialize velocity at inlet
for j=1:def["nj"]
  msm.u[1,j,1] = def["u_inlet"];
end

#! Zero y-component of velocity at outlet
zero_v_at_outlet!(msm::MultiscaleMap, k::Int) = begin
  ni, nj = size(msm.u);
  for j=1:nj
    msm.u[ni,j,2] = 0.0;
  end
end

extract_prof_f(i::Int) = begin

  return (msm::MultiscaleMap) -> begin
    nj = size(msm.u)[2];
    x = Array(Float64, (nj, 3));

    for j=1:nj
      x[j,:] = [j, msm.u[i,j,1], msm.u[i,j,2]];
    end

    return x;
  end;

end

#! Collect callback functions for end of each iteration
callbacks! = [
  zero_v_at_outlet!,
  write_datafile_callback("u", def["stepout"],
    ((msm::MultiscaleMap) -> return msm.u[:,:,1]), "data"),
  write_datafile_callback("v", def["stepout"],
    ((msm::MultiscaleMap) -> return msm.u[:,:,2]), "data"),
  write_datafile_callback("u_mag", def["stepout"], u_mag, "data"),
  write_datafile_callback("prof-at-100", def["stepout"],
    extract_prof_f(100), "data"),
  write_datafile_callback("prof-at-250", def["stepout"],
    extract_prof_f(250), "data"),
  write_datafile_callback("prof-at-500", def["stepout"],
    extract_prof_f(500), "data"),
  (msm::MultiscaleMap, k::Int) -> begin
    if k % 10 == 0
      println("step $k");
    end
  end
]

@profile simulate!(lat, msm, def["col_f"], def["bcs"], def["nsteps"],
  callbacks!);

s = open("prof/main.prof","w")
Profile.print(s,cols = 500)
close(s)
