include("test_multiscale.jl");
include("test_constitutive.jl");
include("test_tracking.jl");

if isfile("lbxflow.jl")
  const main = abspath("lbxflow.jl");
elseif isfile(joinpath("..", "lbxflow.jl"))
  const main = abspath(joinpath("..", "lbxflow.jl"));
else
  error("Cannot find lbxflow.jl. Please run script from test directory");
end

const __run_tests_root__ = dirname(@__FILE__); 

const test_input_files = filter(s -> endswith(s, ".yaml"), readdir(__run_tests_root__));
for test_input_file in test_input_files
  full_path_to_input_file = joinpath(__run_tests_root__, test_input_file);
  run(`julia $main -vf $full_path_to_input_file`);
end
