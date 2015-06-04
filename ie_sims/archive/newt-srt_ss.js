{
  "preamble": "const rho_in = 0.9; const rho_out = 0.885; const datadir = \"data/newt-srt_ss\"; const stepout = 500;",
  "datadir": "data/newt-srt_ss",
  "rhoo": 1.0,
  "dx": 1.0,
  "dt": 1.0,
  "ni": 250,
  "nj": 25,
  "nu": 0.02,
  "col_f": "srt_col_f!",
  "nsteps": 5000,
  "bcs": [
      "north_bounce_back!",
      "south_bounce_back!",
      "begin;
        bind_west_pinlet!(lat) = west_pressure_inlet!(lat, rho_in);
        return bind_west_pinlet!;
      end",
      "begin;
        bind_east_poutlet!(lat) = east_pressure_outlet!(lat, rho_out);
        return bind_east_poutlet!;
      end"
  ],
  "callbacks": [
    "write_datafile_callback(\"u\", stepout,
      ((msm::MultiscaleMap) -> return msm.u[:,:,1]), datadir)",
    "write_datafile_callback(\"v\", stepout,
      ((msm::MultiscaleMap) -> return msm.u[:,:,2]), datadir)",
    "write_datafile_callback(\"u_mag\", stepout, u_mag, 
      datadir)",
    "write_datafile_callback(\"prof-midchan\", stepout,
      extract_prof_callback(125), datadir)",
    "print_step_callback(5)"
  ]
}
