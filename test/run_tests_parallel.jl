using ArgParse;

s = ArgParseSettings();
@add_arg_table s begin
  "--nprocs", "-p"
    help = "Number of processors to use"
    arg_type = Int
    default = 1
  "--skip-serial-tests"
    help = "Skip serial unit tests"
    action = :store_true
end
pargs = parse_args(s);

const nprocs = pargs["nprocs"];
if (nprocs < 1); error("Number of processes must be greater than 0"); end

if !pargs["skip-serial-tests"]
  __run_tests_root__ = dirname(@__FILE__);
  require(abspath(joinpath(__run_tests_root__, "test_multiscale.jl")));
  require(abspath(joinpath(__run_tests_root__, "test_constitutive.jl")));
  require(abspath(joinpath(__run_tests_root__, "test_tracking.jl")));
end

addprocs(nprocs-1);
@everywhere __run_tests_root__ = dirname(@__FILE__);
@everywhere const main = abspath(joinpath(__run_tests_root__, "..", "lbxflow.jl"));
@everywhere const test_input_files = filter(s -> endswith(s, ".yaml"), readdir(__run_tests_root__));
@everywhere run_test(test_input_file) = begin
  full_path_to_input_file = joinpath(__run_tests_root__, test_input_file);
  try
    run(`julia $main -vf $full_path_to_input_file`);
  catch e
    return false;
  end
  return true;
end;

passes = pmap(run_test, test_input_files);

# print results
for (file, passed) in zip(test_input_files, passes)
  if passed
    print_with_color(:green, file, ": passed.\n");
  else
    print_with_color(:red, file, ": failed.\n");
  end
end
