{
  "preamble": "const datadir = \"data/lid-tang/bn001\"; const mu_p = 0.2; const tau_y = 0.00000; const m = 1.0e8; const max_iters = 20; const tol = 0.001;",
  "datadir": "data//lid-tang/bn001",
  "dx": 1.0,
  "dt": 1.0,
  "ni": 100,
  "nj": 100,
  "rhoo": 1.0,
  "nu": 0.2,
  "nsteps": 50000,
  "col_f": "begin;
              curry_mrt_bingham_col_f!(lat, msm) = mrt_bingham_col_f!(lat, msm,
                vikhansky_relax_matrix, mu_p, tau_y, m, max_iters, tol,
                1.0e-11);
              return curry_mrt_bingham_col_f!;
            end",
  "bcs": [
    "begin;
      curry_lid_driven!(lat) = lid_driven!(lat, 1.0);
      return curry_lid_driven!;
    end;",
    "south_bounce_back!",
    "east_bounce_back!",
    "west_bounce_back!"
  ]
  "callbacks": [
    "print_step_callback(25)"
  ],
  "postsim": "(msm::MultiscaleMap) -> begin
                const ni, nj = size(msm.u);
                const xs = linspace(0, 1.0, ni);
                const ys = linspace(0, 1.0, nj);
                writedlm(joinpath(datadir, \"u_table.dsv\"),
                                  [[x[i] y[j] msm.u[i,j,1] msm.u[i,j,2] for i=1:ni, j=1:nj], \",\");
                writedlm(joinpath(datadir, \"u.dsv\"), transpose(msm.u[:,:,1]), \",\");
                writedlm(joinpath(datadir, \"v.dsv\"), transpose(msm.u[:,:,2]), \",\");
              end"
}
