# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

include("lattice.jl");

#! Calculate the p norm distance from equilibrium
macro dist_from_eq_norm(f_neq::Vector{Float64}, p::Int)
  return :(norm(f_neq, p));
end
