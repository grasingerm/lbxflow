{
  "datadir": "data/newt-mrt_5000_500_50_005",
  "rhoo": 5.0,
  "dx": 1.0,
  "dt": 1.0,
  "ni": 500,
  "nj": 50,
  "nu": 0.05,
  "nsteps": 5000,
  "stepout": 500,
  "col_f": "begin;
              const itau = 1.0 / @relax_t(0.02*5.0, 5.0, 1.0, 1.0);
              const S = spdiagm([0.0, itau, itau, 0.0, itau, 0.0, itau,
                itau, itau]);
              curry_mrt_newton_col_f!(lat, msm) = mrt_col_f!(lat, msm, S);
              return curry_mrt_newton_col_f!;
            end",
  "bcs": [
      "north_bounce_back!",
      "south_bounce_back!",
      "begin;
        curry_west_bc!(lat) = west_inlet!(lat, 0.04);
        curry_west_bc!;
      end",
      "east_open!"
  ],
  "u_inlet": 0.04
}
