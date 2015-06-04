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
    "north_bounce_back!",
    "south_bounce_back!",
    "periodic_east_to_west!"
  ],
  "callbacks": [
    "(msm::MultiscaleMap, k::Int) -> begin
      if k % 10 == 0
        const nj = size(msm.u)[2];

        x = linspace(-0.5, 0.5, nj);
        y1 = vec(msm.u[1,:,1]);
        y2 = vec(msm.u[10,:,1]);
        y3 = vec(msm.u[20,:,1]);

        clf();
        plot(x,y1,x,y2,x,y3);
        legend([\"1\",\"10\",\"20\"]);
        xlabel(\"x / width\");
        ylabel(\"ux (lat / sec)\");
      end
    end",
    "print_step_callback(25)"
  ],
  "postsim": "(msm::MultiscaleMap) -> begin
                writedlm(joinpath(datadir, \"ux_profile.dsv\"),
                  extract_prof_callback(10)(msm), \",\");
              end",
  "test_for_term": "is_steadystate_x"
}
