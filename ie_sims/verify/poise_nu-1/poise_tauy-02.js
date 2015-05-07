{
  "preamble": "const pgrad = -1.0; const f = [-pgrad; 0.0]; const datadir = \"data/vpoise_tauy-02\"; const mu_p = 1.0; const tau_y = 0.2; const m = 1.0e8; const max_iters = 50; const tol = 1e-6;",
  "datadir": "data/vpoise_tauy-02",
  "dx": 1.0,
  "dt": 1.0,
  "ni": 25,
  "nj": 21,
  "rhoo": 1.0,
  "nu": 1.0,
  "nsteps": 10000,
  "col_f": "begin;
              curry_mrt_bingham_col_f!(lat, msm) = mrt_bingham_col_f!(lat, msm,
                vikhansky_relax_matrix, mu_p, tau_y, m, max_iters, tol,
                1.0e-11, f);
              return curry_mrt_bingham_col_f!;
            end",
  "bcs": [
    "north_bounce_back!",
    "south_bounce_back!",
    "periodic_east_to_west!"
  ],
  "callbacks": [
    "plot_ux_profile_callback(13,25)",
    "print_step_callback(25)"
  ],
  "postsim": "(msm::MultiscaleMap) -> begin
                writedlm(joinpath(datadir, \"ubar_profile.dsv\"),
                  extract_ubar_prof_callback(13)(msm), \",\");
              end",
  "test_for_term": "is_steadystate_x"
}
