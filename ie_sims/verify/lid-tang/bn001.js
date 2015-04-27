{
  "datadir": "data/lid-tang/bn001",
  "dx": 1.0,
  "dt": 1.0,
  "ni": 200,
  "nj": 200,
  "rhoo": 1.0,
  "nu": 2.0,
  "tau_y": 1.0e-4,
  "m": 1.0e8,
  "max_iters": 20,
  "tol": 0.005,
  "nsteps": 5000,
  "stepout": 500,
  "col_f": "begin;
              curry_mrt_bingham_col_f!(lat, msm) = mrt_bingham_col_f!(lat, msm,
                vikhansky_relax_matrix, 2.0, 1.0e-4, 1.0e8, 20, 5.0e-3, 1.0e-11
                );
              return curry_mrt_bingham_col_f!;
            end",
  "bcs": [
    "begin; 
      curry_lid_driven!(lat) = lid_driven!(lat, 1.0);
      return curry_lid_driven!;
    end;",
    "south_bounce_back!",
    "east_bounce_back!",
    "west_bounce_back!"
  ]
}