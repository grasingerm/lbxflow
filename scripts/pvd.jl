using HDF5, JLD, ProfileView

if length(ARGS) != 1
  println("Usage: pvd.jl path/to/profile.jld");
  exit(1);
end

const profpath = ARGS[1];

if isfile(profpath)
  d = jldopen(profpath) do file
    read(file);
  end
  ProfileView.view(d["bt"], lidict=d["lidict"]);
  println("Press enter to continue...");
  readline(STDIN);
else
  error(profpath, " not found!");
end
