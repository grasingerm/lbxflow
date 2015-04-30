import JSON

#! Load simulation definitions from inputfile
function load_sim_definitions(inputfile::String)
  self = JSON.parse(readall(inputfile));

  # parse and evaluate setup code
  if in("preamble", keys(self))
    preamble = self["preamble"];
    println("evalling: $preamble")
    eval(parse(preamble));
  end

  # parse and evaluate collision function to function pointer
  self["col_f"] = eval(parse(self["col_f"]));

  # parse and evaluate boundary conditions to function pointers
  new_bcs = Array(Function, size(self["bcs"]));
  for (i, bc) in enumerate(self["bcs"])
    println("evalling: $bc");
    new_bcs[i] = eval(parse(bc));
  end
  self["bcs"] = new_bcs;

  # if callbacks exist, parse and evaluate them
  if in("callbacks", keys(self))
    new_callbacks = Array(Function, size(self["callbacks"]));
    for (i, cb) in enumerate(self["callbacks"])
      println("evalling: $cb");
      new_callbacks[i] = eval(parse(cb));
    end
    self["callbacks"] = new_callbacks;
  end

  if in("test_for_term", keys(self))
    self["test_for_term"] = eval(parse(self["test_for_term"]));
  end

  return self;
end
