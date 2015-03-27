const __defaults_root__ = dirname(@__FILE__);
require(abspath(joinpath(__defaults_root__, "boundary.jl")));
require(abspath(joinpath(__defaults_root__, "collision.jl")));

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
    "nsteps" => 5000,
    "stepout" => 500,
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
