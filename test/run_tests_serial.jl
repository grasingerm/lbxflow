if isfile("lbxflow.jl")
  const main = abspath("lbxflow.jl");
elseif isfile(joinpath("..", "lbxflow.jl"))
  const main = abspath(joinpath("..", "lbxflow.jl"));
else
  error("Cannot find lbxflow.jl. Please run script from test directory");
end

const __run_tests_root__  = dirname(@__FILE__); 

const test_input_files    = filter(s -> endswith(s, ".yaml"), readdir(__run_tests_root__));
test_results              = fill(true, length(test_input_files));
for (i, test_input_file) in enumerate(test_input_files)
  full_path_to_input_file = joinpath(__run_tests_root__, test_input_file);
  try
    run(`julia $main -vf $full_path_to_input_file`);
  catch e
    test_results[i] = false; 
  end
end

println("Test results");
println("============");
for (test, result) in zip(test_input_files, results)
  @printf("%40s: ", test);
  if result
    print_with_color(:green, "Passed.");
  else
    print_with_color(:red, "FAILED.");
  end
  print("\n");
end
