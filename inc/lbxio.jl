__lbxio_root__ = dirname(@__FILE__);
require(abspath(joinpath(__lbxio_root__, "multiscale.jl")));

#! Create a callback function for writing to a delimited file
#!
#! \param pre Prefix of filename
#! \param stepout Number of steps in between writing
#! \param A Function for extracting a 2D array from the multiscale map
#! \param delim Delimiter to separate values with
function write_datafile_callback (pre::String, stepout::Int, A::Function, 
  dir=".", delim=',')

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % stepout == 0
      writedlm(joinpath(dir, pre*"_step-$k.dsv"), A(msm), delim);
    end
  end;

end
