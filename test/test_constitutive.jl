using Base.Test
const __test_constitutive_root__ = dirname(@__FILE__);
require(abspath(joinpath(__test_constitutive_root__, "..", "inc", "col", 
                         "constitutive.jl")));
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
