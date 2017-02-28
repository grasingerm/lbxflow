include("../inc/LBXFlow.jl");
using LBXFlow;
using Roots;

h = 64;
H = 2 * h;

for (re, bn) in zip([6.05; 4.42; 3.04; 1.92], [0.68; 1.85; 4.04; 8.52])
  u0 = 0.2 * re / H;
  tau_y = bn * 0.2 * u0 / H;
  #pgrad = fzero(pgrad -> 1 / (2 * 0.2) * -pgrad * (h^2 - (tau_y / pgrad)^2) - tau_y / 0.2 * (h - (tau_y / pgrad)) - u0, -2.0, 0.0);
  pgrad = fzero(pgrad -> (maximum(LBXFlow.analytical_poise_bingham(0.2, tau_y, pgrad, H)) - u0), -1e-5; order=8);
  println("re = $re, bn = $bn, u = $u0 => pgrad = $pgrad, tau_y = $tau_y");
  umax = maximum(LBXFlow.analytical_poise_bingham(0.2, tau_y, pgrad, H));
  println("   Re => $(LBXFlow.reynolds(umax, H, 0.2))");
  println("   Bn => $(LBXFlow.bingham_number(umax, H, 0.2, tau_y))");
end

for (n, re) in zip([0.5; 0.75; 1.25; 1.5], [0.0007; 0.9125; 423.2; 2213])
  u0 = (re * 0.2 / H^n)^(1 / (2 - n));
  #pgrad = fzero(pgrad -> n / (n + 1) * (-1 / 0.2 * pgrad)^(1/n) * h^( (n+1)/n ) - u0, -5.0, 0.0);
  pgrad = try
    fzero(pgrad -> maximum(LBXFlow.analytical_poise_power_law(0.2, n, pgrad, H)) - u0, -1e-5; order=8);
  catch
    fzero(pgrad -> maximum(LBXFlow.analytical_poise_power_law(0.2, n, pgrad, H)) - u0, -1.0, 0.0);
  end
  umax = maximum(LBXFlow.analytical_poise_power_law(0.2, n, pgrad, H));
  println("re = $re, u = $u0 => pgrad = $pgrad");
  println("    Re => $(umax^(2.0 - n) * H^n / 0.2)");
end
