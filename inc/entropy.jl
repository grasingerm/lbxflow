# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

include("lattice.jl");

#_EPS_POS_F = -1e-2;
_EPS_POS_F = -2*eps();

#! Calculate the classical Boltzmann entropy
#!
#! \param lat Lattice
#! \param i   ith index at which to calculate the Boltzmann entropy
#! \param j   jth index at which to calculate the Boltzmann entropy
#! \return    Entropy
function entropy_lat_boltzmann(lat::Lattice, i::Int, j::Int)
  ent = 0.0;
  for k = 1:lat.n
    if lat.f[k, i, j] > 0.0
      ent -= lat.f[k, i, j] * log(lat.f[k, i, j] / lat.w[k]);
    elseif lat.f[k, i, j] < _EPS_POS_F
      error("f[$k, $i, $j] = $(f[$k, $i, $j]), f[$k, $i, $j] < 0.0 in " * "
            $(@__FILE__)");
    end
  end
  return ent
end

#! Calculate the classical Boltzmann entropy
#!
#! \param   lat Lattice
#! \return      Map of Boltzmann entropy over the domain 
function entropy_lat_boltzmann(lat::Lattice)
  ni, nj = size(lat.f, 2), size(lat.f, 3);
  ent = Array{Float64,2}(lat.ni, lat.nj);
  for i = 1:ni, j = 1:nj
    ent[i, j] = lat_boltzmann_entropy(lat, i, j);
  end

  return ent;
end

#! Calculate the classical Boltzmann entropy
#!
#! \param lat Lattice
#! \param f   Particle distribution vector
#! \return    Entropy
function entropy_lat_boltzmann(lat::Lattice, f::AbstractVector{Float64})
  ent = 0.0;
  for k = 1:lat.n
    if f[k] > 0.0
      ent -= f[k] * log(f[k] / lat.w[k]);
    elseif f[k] < _EPS_POS_F
      error("f[$k] = $(f[k]), f[$k] < 0.0 in $(@__FILE__)");
    end
  end
  return ent
end

#! Calculate relative non-equilibrium entropy density
#!
#! \param   f       Particle distributions
#! \param   f_eq    Equilibrium distributions
#! \param   f_neq   Non-equilibrium distributions
#! \return          Relative non-equilibrium entropy density
function entropy_noneq_density(f::AbstractVector{Float64}, 
                               f_eq::AbstractVector{Float64},
                               f_neq::AbstractVector{Float64})
  nk = length(f);
  ds = 0.0;
  for k = 1:nk
    ds += f[k] * log(f[k] / f_eq[k]) - f_neq;
  end

  return ds;
end

#! Calculate relative non-equilibrium entropy density
#!
#! \param   f     Particle distributions
#! \param   f_eq  Equilibrium distributions
#! \return        Relative non-equilibrium entropy density
function entropy_noneq_density(f::AbstractVector{Float64}, 
                               f_eq::AbstractVector{Float64})
  return entropy_noneq_density(f, f_eq, f - f_eq);
end

#! Calculate the quadratic entropy
#!
#! \param   f       Particle distributions
#! \param   f_eq    Equilibrium distributions
#! \param   f_neq   Non-equilibrium distributions
#! \return          Quadratic entropy
function entropy_quadratic(f::AbstractVector{Float64}, 
                           f_eq::AbstractVector{Float64},
                           f_neq::AbstractVector{Float64})
  nk = length(f);
  ds = 0.0;
  for k = 1:nk
    ds += f_neq[k]^2 / (2 * f_eq[k]);
  end
 
  return ds;
end

#! Calculate the quadratic entropy
#!
#! \param   f       Particle distributions
#! \param   f_eq    Equilibrium distributions
#! \param   f_neq   Non-equilibrium distributions
#! \return          Quadratic entropy
function entropy_quadratic(f::AbstractVector{Float64}, 
                           f_eq::AbstractVector{Float64})
  return entropy_quadratic(f, f_eq, f - f_eq);
end
