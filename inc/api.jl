import JSON

#! Load simulation definitions from inputfile
function load_sim_definitions(inputfile::String)
  j = JSON.parse(readall(inputfile));

  # parse and evaluate collision function to function pointer
  j["col_f"] = eval(parse(j["col_f"]));

  # parse and evaluate boundary conditions to function pointers
  new_bcs = Array(Function, size(j["bcs"]));
  for (i, bc) in enumerate(j["bcs"])
    new_bcs[i] = eval(parse(bc));
  end
  j["bcs"] = new_bcs;

  # if callbacks exist, parse and evaluate them
  if in("callbacks", keys(j))
    for (i, cb) in enumerate(j["callbacks"])
      j["callbacks"][i] = eval(parse(cb));
    end
  end

  return j;
end
