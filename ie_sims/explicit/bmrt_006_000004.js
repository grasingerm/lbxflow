{
  "preamble": "const pgrad = -5.2e-6; const f = [-pgrad; 0.0]; const datadir = \"data/bmrt_006_000004\";",
  "datadir": "data/bmrt_006_000004",
  "dx": 1.0,
  "dt": 1.0,
  "ni": 40,
  "nj": 25,
  "rhoo": 1.0,
  "nu": 0.2,
  "mu_p": 0.2,
  "tau_y": 4.0e-5,
  "m": 100000,
  "nsteps": 15000,
  "col_f": "begin;
              curry_mrt_bingham_col_f!(lat, msm) = mrt_fallah_bingham_col_f!(
                lat, msm, vikhansky_relax_matrix, 0.2, 4.0e-5, 1.0e6, 1.0e-11,
                f);
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
              end"
}
