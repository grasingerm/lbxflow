{
  "preamble": "const datadir = \"data/newt-srt_lid\"; const stepout = 1000;",
  "datadir": "data/newt-srt_lid",
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
    "write_datafile_callback(\"u\", stepout,
      ((msm::MultiscaleMap) -> return msm.u[:,:,1]), datadir)",
    "write_datafile_callback(\"v\", stepout,
      ((msm::MultiscaleMap) -> return msm.u[:,:,2]), datadir)",
    "write_datafile_callback(\"u_mag\", stepout, u_mag,
      datadir)",
    "write_datafile_callback(\"rho\", stepout,
      ((msm::MultiscaleMap) -> return msm.rho), datadir)",
    "print_step_callback(50)"
  ],
  "test_for_term": "is_steadystate"
}
