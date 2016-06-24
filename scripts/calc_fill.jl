fltr_bad_val(val::Number) = if (val < 0)
                              0.0
                            elseif (val > 1.0)
                              1.0
                            else
                              val
                            end;
fltr_bad_val(val::AbstractString) = 0;

d = readdlm(ARGS[1], ',');
map!(fltr_bad_val, d);

println(norm(d, Inf));
