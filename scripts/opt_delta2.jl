using Optim

function f(δ)
  r = open("example_sims/for_paper/poise/ehren_bgk-8.yaml");
  s = readall(r);
  close(r);
  w = open("example_sims/for_paper/poise/ehren_bgk-8.yaml", "w");
  write(w, replace(s, r"ds_threshold.*\)", "ds_threshold=$(δ))"));
  close(w);
  run(`julia lbxflow.jl -f example_sims/for_paper/poise/ehren_bgk-8.yaml`);
  e = open("data/for_paper/ehren/bgk-8/opt/rerrors.txt");
  l = readline(e);
  c = parse(strip(split(l, '=')[2]));
  close(e);
  return c;
end

result = optimize(f, 1e-9, 2.0, Optim.Brent());
println("Minimizer: $(Optim.minimizer(result))");
println("Minimum: $(Optim.minimum(result))");
