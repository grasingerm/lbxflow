{
  "preamble": "using PyPlot; const datadir = \"data/newt_channel_open\";",
  "datadir": "data/newt_channel_open",
  "rhoo": 1.0,
  "dx": 1.0,
  "dt": 1.0,
  "ni": 150,
  "nj": 20,
  "nu": 2.0,
  "col_f": "srt_col_f!",
  "nsteps": 10000,
  "bcs": [
    "north_bounce_back!",
    "south_bounce_back!",
    "begin;
        curry_west_bc!(lat) = west_inlet!(lat, 0.004);
        curry_west_bc!;
      end",
    "curry_east_pressure!(lat) = east_pressure!(lat, 0.0)"
  ],
  "callbacks": [
    "(msm::MultiscaleMap, k::Int) -> begin
      if k % 20 == 0
        const nj = size(msm.u)[2];

        x = linspace(-0.5, 0.5, nj);
        y1 = vec(msm.u[3,:,1]);
        y2 = vec(msm.u[75,:,1]);
        y3 = vec(msm.u[148,:,1]);

        clf();
        plot(x,y1,x,y2,x,y3);
        legend([\"3\",\"75\",\"148\"]);
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
