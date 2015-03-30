{
  "dx": 1.0,
  "dt": 1.0,
  "ni": 500,
  "nj": 20,
  "rhoo": 1.0,
  "nu": 0.015,
  "mu_p": 0.015,
  "tau_y": 0.03,
  "m": 1000000,
  "max_iters": 20,
  "tol": 0.005,
  "nsteps": 5000,
  "stepout": 500,
  "col_f": "begin;
              curry_mrt_bingham_col_f!(lat, msm) = mrt_bingham_col_f!(lat, msm, 
                vikhansky_relax_matrix, 0.015, 0.03, 1.0e6, 20, 5.0e-3)
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
