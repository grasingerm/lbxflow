{
  "datadir": "data/bingham-mrt_5000_500_20_012_000064",
  "dx": 1.0,
  "dt": 1.0,
  "ni": 500,
  "nj": 20,
  "rhoo": 1.0,
  "nu": 0.12,
  "mu_p": 0.12,
  "tau_y": 64.0e-6,
  "m": 1000000,
  "max_iters": 20,
  "tol": 0.005,
  "nsteps": 5000,
  "stepout": 500,
  "col_f": "begin;
              curry_mrt_bingham_col_f!(lat, msm) = mrt_bingham_col_f!(lat, msm,
                vikhansky_relax_matrix, 0.12, 64.0e-6, 1.0e6, 20, 5.0e-3, 1.0e-9,
                0.501, 8.0);
              return curry_mrt_bingham_col_f!;
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
