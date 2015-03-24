require("boundary.jl");
require("collision.jl");

# Defaults for lbx simulation
function lbx_defaults()
  u_inlet = 0.02;

  curry_west_bc!(lat) = west_inlet!(lat, u_inlet)

  return [
    "rhoo" => 5.0,
    "dx" => 1.0,
    "dt" => 1.0,
    "ni" => 1000,
    "nj" => 40,
    "nu" => 0.02,
    "nsteps" => 25000,
    "stepout" => 1000,
    "col_f" => srt_col_f!,
    "bcs" => [
      north_bounce_back!,
      south_bounce_back!,
      curry_west_bc!,
      east_open!
    ],
    "u_inlet" => u_inlet
  ];
end
