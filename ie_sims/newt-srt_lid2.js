{
  "preamble": "using PyPlot; const datadir = \"data/newt-srt_lid2\";",
  "datadir": "data/newt-srt_lid2",
  "rhoo": 5.0,
  "dx": 1.0,
  "dt": 1.0,
  "ni": 100,
  "nj": 100,
  "nu": 0.01,
  "col_f": "srt_col_f!",
  "nsteps": 50000,
  "bcs": [
    "begin;
      curry_lid_driven!(lat) = lid_driven!(lat, 0.1);
      return curry_lid_driven!;
    end;",
    "south_bounce_back!",
    "east_bounce_back!",
    "west_bounce_back!"
  ],
  "callbacks": [
    "plot_streamlines_callback(10, (0.45, -0.075), 0.1)",
    "print_step_callback(25)"
  ],
  "postsim": "(msm::MultiscaleMap) -> begin
                writedlm(joinpath(datadir, \"ux.dsv\"), msm.u[:,:,1], \",\");
                writedlm(joinpath(datadir, \"uy.dsv\"), msm.u[:,:,2], \",\");
              end",
  "test_for_term": "is_steadystate"
}
