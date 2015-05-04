{
  "preamble": "using PyPlot; const pgrad = -5.2e-6; const f = [-pgrad; 0.0]; const datadir = \"data/poise_tauy-000012\"; const mu_p = 0.2; const tau_y = 0.00012; const m = 1.0e8; const max_iters = 150; const tol = 1e-6;",
  "datadir": "data/poise_tauy-000012",
  "dx": 1.0,
  "dt": 1.0,
  "ni": 20,
  "nj": 20,
  "rhoo": 1.0,
  "nu": 0.2,
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
    "plot_ux_profile_callback(10, 25, 0.1)",
    "print_step_callback(25)"
  ],
  "postsim": "(msm::MultiscaleMap) -> begin
                writedlm(joinpath(datadir, \"ubar_profile.dsv\"),
                  extract_ubar_prof_callback(10)(msm), \",\");
              end"
}
