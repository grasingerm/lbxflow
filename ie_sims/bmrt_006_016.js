{
  "preamble": "const datadir = \"data/bmrt_006_008-2\";",
  "datadir": "data/bmrt_006_008-2",
  "dx": 1.0,
  "dt": 1.0,
  "ni": 500,
  "nj": 25,
  "rhoo": 1.0,
  "nu": 0.06,
  "mu_p": 0.06,
  "tau_y": 16.0e-3,
  "m": 10000000,
  "nsteps": 15000,
  "col_f": "begin;
              curry_mrt_bingham_col_f!(lat, msm) = mrt_fallah_bingham_col_f!(
                lat, msm, vikhansky_relax_matrix, 0.06, 16.0e-3, 1.0e8, 1.0e-11);
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
  "u_inlet": 0.04,
  "callbacks": [
    "print_step_callback(100)"
  ],
  "postsim": "(msm::MultiscaleMap) -> begin
                writedlm(joinpath(datadir, \"ux_profile.dsv\"),
                  extract_prof_callback(250)(msm), \",\");
              end",
  "test_for_term": "is_steadystate_x"
}
