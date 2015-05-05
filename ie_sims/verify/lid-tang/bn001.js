{
  "preamble": "const datadir = \"data/lid-tang/bn001\"; const U = 0.01; const mu_p = 0.0075; const tau_y = 0.00000001; const m = 1.0e8; const max_iters = 50; const tol = 0.001; const Re = U * 75 / mu_p; const Bn = tau_y * 75 / (mu_p * U); info(\"Re = $Re, Bn = $Bn\"); using PyPlot;",
  "datadir": "data/lid-tang/bn001",
  "dx": 1.0,
  "dt": 1.0,
  "ni": 75,
  "nj": 75,
  "rhoo": 1.0,
  "nu": 0.0075,
  "nsteps": 50000,
  "col_f": "begin;
              curry_mrt_bingham_col_f!(lat, msm) = mrt_bingham_col_f!(lat, msm,
                vikhansky_relax_matrix, mu_p, tau_y, m, max_iters, tol,
                1.0e-11);
              return curry_mrt_bingham_col_f!;
            end",
  "bcs": [
    "begin;
      curry_lid_driven!(lat) = lid_driven!(lat, U);
      return curry_lid_driven!;
    end;",
    "south_bounce_back!",
    "east_bounce_back!",
    "west_bounce_back!"
  ],
  "callbacks": [
  "plot_streamlines_callback(10, (0.45, -0.075),
      joinpath(datadir, \"lid-tang_sl\"), 0.025)",
    "print_step_callback(25)"
  ],
  "postsim": "(msm::MultiscaleMap) -> begin
                const ni, nj = size(msm.u);
                const xs = linspace(0, 1.0, ni);
                const ys = linspace(0, 1.0, nj);

                writedlm(joinpath(datadir, \"u.dsv\"), transpose(msm.u[:,:,1]), \",\");
                writedlm(joinpath(datadir, \"v.dsv\"), transpose(msm.u[:,:,2]), \",\");
                writedlm(joinpath(datadir, \"u_midcav.dsv\"), [vec(msm.u[round(ni/2),:,1]) y], \",\");
                writedlm(joinpath(datadir, \"v_midcav.dsv\"), [x vec(msm.u[:,round(nj/2),2])], \",\");
              end"
}
