@everywhere const main = begin;
  if isfile("lbxflow.jl")
    abspath("lbxflow.jl");
  elseif isfile(joinpath("..", "lbxflow.jl"))
    abspath(joinpath("..", "lbxflow.jl"));
  else
    error("Cannot find lbxflow.jl. Please run script from test directory");
  end
end

@everywhere const __run_tests_root__  = dirname(@__FILE__); 

const test_input_files    = filter(s -> endswith(s, ".yaml"), readdir(__run_tests_root__));

function did_pass(file::AbstractString)
  full_path_to_input_file = joinpath(__run_tests_root__, test_input_file);
  println("Parsing and running input file: ", full_path_to_input_file);
  passed = true;
  try
    run(`julia --color=yes $main -vf $full_path_to_input_file`);
  catch e
    passed = false;
  end
  return passed;
end

println(test_input_files);

test_results              = pmap(did_pass, test_input_files);
println("Test results");
println("============");
for (test, result) in zip(test_input_files, test_results)
  @printf("%60s: ", test);
  if result
    print_with_color(:green, "Passed.");
  else
    print_with_color(:red, "FAILED.");
  end
  print("\n");
end
