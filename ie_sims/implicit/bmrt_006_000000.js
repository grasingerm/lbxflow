{
  "preamble": "const pgrad = -5.2e-6; const f = [-pgrad; 0.0]; const datadir = \"data/bmrt_006_000000\";",
  "datadir": "data/bmrt_006_000000",
  "dx": 1.0,
  "dt": 1.0,
  "ni": 128,
  "nj": 64,
  "rhoo": 1.0,
  "nu": 0.2,
  "mu_p": 0.2,
  "tau_y": 0.0e-5,
  "m": 100000,
  "nsteps": 15000,
  "col_f": "begin;
              curry_mrt_bingham_col_f!(lat, msm) = mrt_bingham_col_f!(lat, msm,
                vikhansky_relax_matrix, 0.2, 0.0, 1.0e8, 7, 1e-5,
                1e-9, f, 0.75);
              return curry_mrt_bingham_col_f!;
            end",
  "bcs": [
      "north_bounce_back!",
      "south_bounce_back!",
      "periodic_east_to_west!"
  ],
  "callbacks": [
    "print_step_callback(100)"
  ],
  "postsim": "(msm::MultiscaleMap) -> begin
                writedlm(joinpath(datadir, \"ux_profile.dsv\"),
                  extract_ux_prof_callback(20)(msm), \",\");
                writedlm(joinpath(datadir, \"ubar_profile.dsv\"),
                  extract_ubar_prof_callback(20)(msm), \",\");
              end",
  "test_for_term": "is_steadystate_x"
}
