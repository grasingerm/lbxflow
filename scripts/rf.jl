top_d = ARGS[1];

fltr_bad_val(val::Number) = if (val < 0)
                              0.0
                            elseif (val > 1.0)
                              1.0
                            else
                              val
                            end;
fltr_bad_val(val::AbstractString) = 0;

function _r(args...)
  paths = readdir(joinpath(args...));
  for path in paths
    if isdir(path)
      _r(args..., path);
    elseif path == "ff.csv"
      report(args..., path);
    end
  end
end

function report(args...)
  ffdata = readdlm(joinpath(args...), ',');
  map!(fltr_bad_val, ffdata);

  ni, nj = size(ffdata);
  fill_max = 0;
  ma_max = 0;
  for j=1:nj
    ts = sub(ffdata, j, :);
    fill_sum = sum(ts);
    if fill_sum > fill_max
      fill_max = fill_sum;
    end
  end
  println(join(args), " fill_max: ", fill_max);
end

_r(top_d);
