{
  "datadir": "data/newt-srt_5000_500_50_002",
  "rhoo": 1.0,
  "dx": 1.0,
  "dt": 1.0,
  "ni": 250,
  "nj": 25,
  "nu": 0.02,
  "nsteps": 25000,
  "test_for_term": "is_stateadystate",
  "stepout": 500,
  "col_f": "srt_col_f!",
  "bcs": [
      "north_bounce_back!",
      "south_bounce_back!",
      "begin;
        const rho_in = @rho(0.1, @c_ssq(1.0, 1.0));
        bind_west_pinlet!(lat) = west_pressure_inlet!(lat, rho_in);
        bind_west_pinlet!;
      end",
      "begin;
        const rho_out = @rho(0.05, @c_ssq(1.0, 1.0));
        bind_east_poutlet!(lat) = east_pressure_outlet!(lat, rho_out);
        bind_east_poutlet!;
      end"
  ]
}
