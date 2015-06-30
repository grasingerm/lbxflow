const __run_tests_root__ = dirname(@__FILE__);
require(abspath(joinpath(__run_tests_root__, "test_tracking.jl")));

const main = abspath(joinpath(__run_tests_root__, "..", "lbxflow.jl"));
const test_input_files = filter(s -> endswith(s, ".yaml"), readdir(__run_tests_root__));
for test_input_file in test_input_files
  full_path_to_input_file = joinpath(__run_tests_root__, test_input_file);
  run(`julia $main -vf $full_path_to_input_file`);
end
