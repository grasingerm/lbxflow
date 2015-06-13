{
  "preamble": "const pgrad = -5.2e-6; const f = [-pgrad; 0.0]; const datadir = \"data/newt_006_000000\";",
  "datadir": "data/newt_006_000000",
  "dx": 1.0,
  "dt": 1.0,
  "ni": 40,
  "nj": 25,
  "rhoo": 1.0,
  "nu": 0.2,
  "nsteps": 15000,
  "col_f": "begin;
      bind_srt_col_f!(lat, msm) = srt_guo_col_f!(lat, msm, f);
      return bind_srt_col_f!;
    end;",
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
