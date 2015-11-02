# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

const __forcing_root__ = dirname(@__FILE__);
require(abspath(joinpath(__forcing_root__, "..", "lattice.jl")));
require(abspath(joinpath(__forcing_root__, "..", "multiscale.jl")));
require(abspath(joinpath(__forcing_root__, "..", "sim", "simtypes.jl")));

#! Initialize a sukop forcing function
#!
#! \param F Body force vector
#! \return (momentum_function, forcing_function)
function init_sukop_Fk(F::Vector{Float64})
  return (
    (lat::Lattice, u::Vector{Float64}) -> begin
      return u;
    end,
    (lat::Lattice, omega::AbstractFloat, u::Vector{Float64}, k::Int) -> begin
      return lat.w[k] * lat.dt / lat.cssq * dot(F, lat.c[:,k]);
    end
    );
end

#! Initialize a guo forcing function
#!
#! \param F Body force vector
#! \return (momentum_function, forcing_function)
function init_guo_Fk(F::Vector{Float64})
  return (
    (lat::Lattice, u::Vector{Float64}) -> begin
      return u + lat.dt / 2.0 * F;
    end,
    (lat::Lattice, omega::AbstractFloat, u::Vector{Float64}, k::Int) -> begin
      return ((1 - 0.5 * omega) * lat.w[k] * dot(((lat.c[:,k] - u) / lat.cssq +
              dot(lat.c[:,k], u) / (lat.cssq * lat.cssq) * lat.c[:,k]), F));
    end
    );
end

#! Initialize a korner forcing function
#!
#! \param F Body force vector
#! \return (momentum_function, forcing_function)
function init_korner_Fk(F::Vector{Float64})
  return (
    (lat::Lattice, u::Vector{Float64}) -> begin
      return u;
    end,
    (lat::Lattice, omega::AbstractFloat, u::Vector{Float64}, k::Int) -> begin
      return lat.w[k] * lat.dt / lat.cssq * dot(lat.c[:,k], F);
    end
    );
end
