{
  "datadir": "data/newt-srt_ss",
  "rhoo": 1.0,
  "dx": 1.0,
  "dt": 1.0,
  "ni": 250,
  "nj": 25,
  "nu": 0.02,
  "nsteps": 5000,
  "stepout": 500,
  "col_f": "srt_col_f!",
  "bcs": [
      "north_bounce_back!",
      "south_bounce_back!",
      "begin;
        const rho_in = 3 * 0.1;
        bind_west_pinlet!(lat) = west_pressure_inlet!(lat, rho_in);
        bind_west_pinlet!;
      end",
      "begin;
        const rho_out = 3 * 0.05;
        bind_east_poutlet!(lat) = east_pressure_outlet!(lat, rho_out);
        bind_east_poutlet!;
      end"
  ]
}
