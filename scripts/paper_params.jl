using Roots;

h = 64;

for (re, bn) in zip([6.05; 4.42; 3.04; 1.92], [0.68; 1.85; 4.04; 8.52])
  u0 = 0.2 * re / h;
  tau_y = bn * 0.2 * u0 / h;
  pgrad = fzero(pgrad -> 1 / (2 * 0.2) * -pgrad * (h^2 - (tau_y / pgrad)^2) - tau_y / 0.2 * (h - (tau_y / pgrad)) - u0, -2.0, 0.0);
  println("re = $re, bn = $bn, u = $u0 => pgrad = $pgrad, tau_y = $tau_y");
  #println("   pgrad = $pgrad => u = $(1 / (2 * 0.2) * -pgrad * (h^2 - (tau_y / pgrad)^2) - tau_y / 0.2 * (h - (tau_y / pgrad))) ?= $u0");
end

for (n, re) in zip([0.5; 0.75; 1.25; 1.5], [0.0007; 0.9125; 423.2; 2213])
  u0 = (re * 0.2 / h^n)^(1 / (2 - n));
  pgrad = fzero(pgrad -> n / (n + 1) * (-1 / 0.2 * pgrad)^(1/n) * h^( (n+1)/n ) - u0, -5.0, 0.0);
  println("re = $re, u = $u0 => pgrad = $pgrad");
end
