# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Equilibrium frequency distribution for incompressible Newtonian flow
#!
#! \param lat D2Q9 lattice
#! \param rho Macroscopic density at lattice site
#! \param u Macroscopic flow at lattice site
#! \param k Lattice velocity vector index
#! \return Equilibrium frequency
function feq_incomp(lat::LatticeD2Q9, rho::AbstractFloat, u::Vector{Float64},
                    k::Int)
  const cssq = 1/3;
  const ckdotu = dot(lat.c[:,k], u);

  return rho * lat.w[k] * (1.0 + ckdotu/(cssq) + 0.5*(ckdotu*ckdotu)/(cssq*cssq)
                           - 0.5 * dot(u, u) / (cssq));
end

#! Equilibrium frequency distribution for incompressible Newtonian flow
#!
#! \param lat D2Q9 lattice
#! \param rho Density at lattice site
#! \param rho_0 Nominal or average density
#! \param u Macroscopic flow at lattice site
#! \param k Lattice velocity vector index
#! \return Equilibrium frequency
function feq_incomp_HL(lat::LatticeD2Q9, rho::AbstractFloat,
                       rho_0::AbstractFloat, u::Vector{Float64}, k::Int)
  const cssq = 1/3;
  const ckdotu = dot(lat.c[:,k], u);

  return (lat.w[k] * (rho + rho_0 * (ckdotu/(cssq)
              + 0.5*(ckdotu*ckdotu)/(cssq*cssq)
              - 0.5 * dot(u, u) / (cssq))));
end

#! Binds rho_0 to an HL equilibrium function
function init_feq_incomp_HL(rho_0::AbstractFloat)
  return ((lat::Lattice, rho::AbstractFloat, u::Vector{Float64},
          k::Int) -> feq_incomp_HL(lat, rho, rho_0, u, k));
end

#! Equilibrium frequency distribution for incompressible Newtonian flow
#!
#! \param lat D2Q4 lattice
#! \param rho Macroscopic density at lattice site
#! \param u Macroscopic flow at lattice site
#! \return Equilibrium frequency
function feq_incomp(lat::LatticeD2Q4, rho::AbstractFloat, u::Vector{Float64},
                    k::Int)
  const cssq = 1/2;
  const ckdotu = dot(lat.c[:,k], u);

  return rho * lat.w[k] * (1 + 3 * ckdotu);
end

#! Equilibrium frequency distribution for incompressible Newtonian flow
#!
#! \param lat D2Q4 lattice
#! \param rho Density at lattice site
#! \param rho_0 Nominal or average density
#! \param u Macroscopic flow at lattice site
#! \param k Lattice velocity vector index
#! \return Equilibrium frequency
function feq_incomp_HL(lat::LatticeD2Q4, rho::AbstractFloat,
                       rho_0::AbstractFloat, u::Vector{Float64}, k::Int)
  const cssq = 1/2;
  const ckdotu = dot(lat.c[:,k], u);

  return lat.w[k] * (rho + rho_0 * 3 * ckdotu);
end

#! Equilibrium frequency distribution that maximizes entropy
#! (Gorban and Packwood, Physica A 2014)
#!
#! \param lat D2Q9 lattice
#! \param rho Macroscopic density at lattice site
#! \param u Macroscopic flow at lattice site
#! \param k Lattice velocity vector index
#! \return Equilibrium frequency distribution that maximizes entropy
function feq_incomp_max_entropy(lat::LatticeD2Q9, rho::AbstractFloat, 
                                u::Vector{Float64}, k::Int)
  const nj  =   length(u);
  prod      =   1.0;

  for j=1:nj
    x1      =   sqrt(1 + 3.0*u[j]^2);
    prod    *=  (2 - x1) * ( (2*u[j] + x1) / (1 - u[j]) )^lat.c[j,k]; 
  end

  return lat.w[k] * rho * prod;
end
