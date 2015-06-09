{
  "preamble": "using PyPlot; const pgrad = -1.102e-3; const f = [-pgrad; 0.0]; const datadir = \"data/sukop\";",
  "datadir": "data/sukop",
  "rhoo": 1.0,
  "dx": 1.0,
  "dt": 1.0,
  "ni": 20,
  "nj": 12,
  "nu": 0.166666666,
  "col_f":
    "begin;
      bind_srt_col_f!(lat, msm) = srt_guo_col_f!(lat, msm, f);
      return bind_srt_col_f!;
    end;",
  "nsteps": 10000,
  "bcs": [
    "north_halfa_bounce_back!",
    "south_halfa_bounce_back!",
    "periodic_east_to_west!"
  ],
  "callbacks": [
    "plot_ux_profile_callback(10, 10, 0.0)",
    "print_step_callback(25)"
  ],
  "postsim": "(msm::MultiscaleMap) -> begin
                writedlm(joinpath(datadir, \"ux_profile.dsv\"),
                  extract_prof_callback(10)(msm), \",\");
              end",
  "test_for_term": "is_steadystate_x"
}
