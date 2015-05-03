{
  "preamble": "using PyPlot; const rho_in = 0.9; const rho_out = 0.885; const pgrad = -0.0125; const f = [-pgrad; 0.0]; const datadir = \"data/newt-srt_poiseuille_pbc\"; const stepout = 1000;",
  "datadir": "data/newt-srt_poiseuille_pbc",
  "rhoo": 1.0,
  "dx": 1.0,
  "dt": 1.0,
  "ni": 20,
  "nj": 20,
  "nu": 0.1,
  "col_f":
    "begin;
      bind_srt_col_f!(lat, msm) = srt_col_f!(lat, msm, f);
      return bind_srt_col_f!;
    end;",
  "nsteps": 10000,
  "bcs": [
    "north_bounce_back!",
    "south_bounce_back!",
    "periodic_east_to_west!"
  ],
  "callbacks": [
    "plot_umag_contour_callback(1, 0.1)",
    "print_step_callback(25)"
  ],
  "test_for_term": "is_steadystate_x"
}
