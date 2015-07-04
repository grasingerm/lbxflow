using Base.Test
const __test_constitutive_root__ = dirname(@__FILE__);
require(abspath(joinpath(__test_constitutive_root__, "..", "inc", "col", 
                         "constitutive.jl")));
require(abspath(joinpath(__test_constitutive_root__, "..", "inc", "col", 
                         "equilibrium.jl")));
require(abspath(joinpath(__test_constitutive_root__, "..", "inc",
                         "multiscale.jl")));
require(abspath(joinpath(__test_constitutive_root__, "..", "inc", "sim",
                         "simtypes.jl")));

if length(ARGS) == 1
  const skip = int(ARGS[1]);
else
  const skip = 0;
end

f(mu_p, tau_y, m, gamma) = mu_p + tau_y / gamma * (1 - exp(-m*gamma));

# our consitutive model should be able to accurately
# characterize yield stress = 0
if skip < 1
  println("Testing @mu_papanstasiou macro...");
  for i=1:10000
    mu_p = rand(0:0.00001:10);
    tau_y = rand(0:0.00001:10);
    m = rand(4:12);
    gamma = rand(0.000001:0.000001:5);

    #println("mu_p = $mu_p, tau_y = 0.0, m = $m, gamma = $gamma");
    @test_approx_eq mu_p @mu_papanstasiou(mu_p, 0.0, m, gamma)

    #println("mu_p = $mu_p, tau_y = $tau_y, m = $m, gamma = $gamma");
    @test_approx_eq f(mu_p, tau_y, m, gamma) @mu_papanstasiou(mu_p, tau_y, m, gamma)
  end
  println("Tests passed.");
end

if skip < 2
  println("Testing srt constitutive relationships...");
  for i=1:1000
    if i % 10 == 0; println("test $i"); end
    mu_p = rand(0:0.00001:10);
    m = rand(4:12);
    gamma = rand(0.000001:0.000001:5);
    ni = rand(1:100);
    nj = rand(1:100);
    rho = rand(1:0.1:7.5);

    lat = LatticeD2Q9(1.0, 1.0, ni, nj, rho);
    lat.f += rand(lat.n, ni, nj); # add some random noise
    msm = MultiscaleMap(mu_p, lat, rho);
    sim = Sim(lat, msm);
    map_to_macro!(lat, msm);

    ccs = init_constit_srt_const(mu_p);
    lcs = init_constit_srt_const_local();
    bcse = init_constit_srt_bingham_explicit(mu_p, 0.0, m, 1e-9,
                                             rand(0.5:0.1:1.5));
    bcsi = init_constit_srt_bingham_implicit(mu_p, 0.0, m, 1e-9, rand(3:50), 
                                             10.0 ^ (rand(-8:-2)), rand(0.5:0.1:1.5));

    fneq = Array(Float64, lat.n);
    for j=1:nj, i=1:ni
      for k = 1:lat.n
        fneq[k] = lat.f[k,i,j] - feq_incomp(lat, msm.rho[i,j], msm.u[:,i,j], k);
      end
      try
        @test_approx_eq mu_p ccs(sim, fneq, i, j)
        @test_approx_eq mu_p lcs(sim, fneq, i, j)
        @test_approx_eq mu_p bcse(sim, fneq, i, j)
        @test_approx_eq mu_p bcsi(sim, fneq, i, j)
      catch e
        println("sim: ", sim, "feq: ", feq, "i: ", i, "j: ", j);
        println("mu_p: ", mu_p, "m: ", m);
        showerror(STDERR, e);
        println();
        Base.show_backtrace(STDERR, catch_backtrace()); # display callstack
        rethrow(e);
      end
    end
  end
  println("Tests passed");
end

if skip < 3
  println("Testing bingham srt with bounds...");
  for i=1:10000
    if i % 10 == 0; println("test $i"); end
    mu_p = rand(0:0.00001:10);
    m = rand(4:12);
    gamma = rand(0.000001:0.000001:5);
    ni = rand(1:100);
    nj = rand(1:100);
    rho = rand(1:0.1:7.5);
    rt_min = rand(0.5:0.1:2.5);
    rt_max = rand(rt_min+0.1:0.1:10.0);
    tau_y = rand(1.0e-5:1.0e-5:1.0);

    lat = LatticeD2Q9(1.0, 1.0, ni, nj, rho);
    lat.f += rand(lat.n, ni, nj); # add some random noise
    msm = MultiscaleMap(mu_p, lat, rho);
    sim = Sim(lat, msm);
    map_to_macro!(lat, msm);

    bcse = init_constit_srt_bingham_explicit(mu_p, tau_y, rt_min, rt_max, m, 1e-9,
                                             rand(0.5:0.1:1.5));
    bcsi = init_constit_srt_bingham_implicit(mu_p, tau_y, rt_min, rt_max, m, 1e-9,
                                             rand(3:50), 10.0 ^ (rand(-8:-2)),
                                             rand(0.5:0.1:1.5));

    fneq = Array(Float64, lat.n);
    for j=1:nj, i=1:ni
      for k = 1:lat.n
        fneq[k] = lat.f[k,i,j] - feq_incomp(lat, msm.rho[i,j], msm.u[:,i,j], k);
      end
      try
        @test (bcse(sim, fneq, i, j) + 1e-8 >= @nu(1.0/rt_min, lat.cssq, lat.dt) &&
               bcse(sim, fneq, i, j) - 1e-8 <= @nu(1.0/rt_max, lat.cssq, lat.dt));
        @test (bcsi(sim, fneq, i, j) + 1e-8 >= @nu(1.0/rt_min, lat.cssq, lat.dt) &&
               bcsi(sim, fneq, i, j) - 1e-8 <= @nu(1.0/rt_max, lat.cssq, lat.dt));
      catch e
        println("sim: ", sim, "fneq: ", fneq, "i: ", i, "j: ", j);
        println("mu_p: ", mu_p, " tau_y: ", tau_y, " rt_min: ", rt_min, " rt_max: ", rt_max, " m: ", m);
        println(@nu(1.0/rt_min, lat.cssq, lat.dt), " <= ", bcse(sim, fneq, i, j), " <= ",  @nu(1.0/rt_max, lat.cssq, lat.dt));
        println(@nu(1.0/rt_min, lat.cssq, lat.dt), " <= ", bcsi(sim, fneq, i, j), " <= ",  @nu(1.0/rt_max, lat.cssq, lat.dt));
        showerror(STDERR, e);
        println();
        Base.show_backtrace(STDERR, catch_backtrace()); # display callstack
        rethrow(e);
      end
    end
  end
  println("Tests passed");
end
