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

f(mu_p, tau_y, m, gamma) = mu_p + tau_y / gamma * (1 - exp(-m*gamma));

# our consitutive model should be able to accurately
# characterize yield stress = 0
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

  ccs = init_constit_const(mu_p);
  lcs = init_constit_const_local();
  bcse = init_constit_srt_bingham_explicit(mu_p, 0.0, m, 1e-9,
                                           rand(0.5:0.1:1.5));
  bcsi = init_constit_srt_bingham_implicit(mu_p, 0.0, m, 1e-9, rand(3:50), 
                                           10.0 ^ (rand(-8:-2)), rand(0.5:0.1:1.5));

  fneq = Array(Float64, lat.n);
  for j=1:nj, i=1:ni
    for k = 1:lat.n
      fneq[k] = lat.f[k,i,j] - feq_incomp(lat, msm.rho[i,j], msm.u[:,i,j], k);
      @test_approx_eq mu_p ccs(sim, fneq, i, j)
      @test_approx_eq mu_p lcs(sim, fneq, i, j)
      @test_approx_eq mu_p bcse(sim, fneq, i, j)
      @test_approx_eq mu_p bcsi(sim, fneq, i, j)
    end
  end
end
println("Tests passed");
