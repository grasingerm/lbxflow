# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

@everywhere import HDF5, JLD

function __create_anon_counter()
  name = parse("__" * join(rand('a':'z', 100)) * "__");
  eval(quote
        global $name = 0.0;
        return (x::Real) -> begin; global $name += x; end;
       end);
end

macro create_callback(step, code)
  return quote
    (sim::AbstractSim, k::Real) -> begin
      if k % $step < sim.Δt
        $code
      end
    end
  end
end

include(joinpath("io", "readwrite.jl"));
include(joinpath("io", "animate.jl"));
