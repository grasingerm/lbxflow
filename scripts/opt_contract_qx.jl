using Optim

const input_fpath = "publication_sims/for_paper/poise/for-opt_median_bgk-12.yaml"
const indent = "                  "; 

function f(x)
  r = open(input_fpath);
  lines = readlines(r);
  close(r);
  w = open(input_fpath, "w");
  for line in lines
    if startswith(line, indent * "scale=")
      write(w, indent * "scale=(args...) -> scale_root_median(args...) / $(x[1]),\n");
    elseif startswith(line, indent * "stds=")
      write(w, indent * "stds=$(x[2]),\n");
    elseif startswith(line, indent * "ds_threshold=")
      write(w, indent * "ds_threshold=$(x[3])\n");
    elseif startswith(line, indent * "diss=")
      write(w, indent * "diss=(args...)->contract_qx!(args...; weight_a=$(x[4])),\n");
    else
      write(w, line);
    end
  end
  close(w);
  run(`julia lbxflow.jl -f $input_fpath`);
  e = open("data/for_paper/median/bgk-12/opt/rerrors.txt");
  l = readline(e);
  c = parse(strip(split(l, '=')[2]));
  close(e);
  return c;
end

result = optimize(f, [16.0; 2.7; 1e-5; 0.5], Optim.NelderMead(),
                  OptimizationOptions(x_tol=1e-12, f_tol=1e-12,
                  iterations=1000, store_trace=true,
                  show_trace=true, extended_trace=true, show_every=10)
                  );
