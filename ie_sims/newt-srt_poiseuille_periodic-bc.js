{
  "preamble": "const rho_in = 0.9; const rho_out = 0.885; const pgrad = (rho_out - rho_in) / 3.0; const f = [-pgrad; 0.0]; const datadir = \"data/newt-srt_poiseuille_pbc\"; const stepout = 1000;",
  "datadir": "data/newt-srt_poiseuille_pbc",
  "rhoo": 5.0,
  "dx": 1.0,
  "dt": 1.0,
  "ni": 20,
  "nj": 20,
  "nu": 0.02,
  "col_f":
    "begin;
      bind_srt_col_f!(lat, msm) = srt_col_f!(lat, msm, f);
      return bind_srt_col_f!;
    end;",
  "nsteps": 50000,
  "bcs": [
    "north_bounce_back!",
    "south_bounce_back!",
    "periodic_east_to_west!"
  ],
  "callbacks": [
    "write_datafile_callback(\"u\", stepout,
      ((msm::MultiscaleMap) -> return msm.u[:,:,1]), datadir)",
    "write_datafile_callback(\"v\", stepout,
      ((msm::MultiscaleMap) -> return msm.u[:,:,2]), datadir)",
    "write_datafile_callback(\"u_mag\", stepout, u_mag,
      datadir)",
    "write_datafile_callback(\"prof-midchan\", stepout,
      extract_prof_callback(10), datadir)",
    "print_step_callback(25)"
  ],
  "test_for_term": "is_steadystate_x"
}
