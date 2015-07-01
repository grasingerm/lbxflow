const __equilibrium_root__ = dirname(@__FILE__);
require(abspath(joinpath(__equilibrium_root__, "..", "lattice.jl")));
require(abspath(joinpath(__equilibrium_root__, "..", "multiscale.jl")));

#! Equilibrium frequency distribution for incompressible Newtonian flow
#!
#! \param lat D2Q9 lattice
#! \param rho Macroscopic density at lattice site
#! \param u Macroscopic flow at lattice site
#! \return Equilibrium frequency
function feq_incomp(lat::LatticeD2Q9, rho::FloatingPoint, u::Vector{Float64},
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
#! \return Equilibrium frequency
function feq_incomp_HL(lat::LatticeD2Q9, rho::FloatingPoint,
                       rho_0::FloatingPoint, u::Vector{Float64})
  const cssq = 1/3;
  const ckdotu = dot(lat.c[:,k], u);

  return w * (rho + rho_0 * (ckdotu/(cssq) + 0.5*(ckdotu*ckdotu)/(cssq*cssq)
              - 0.5 * dot(u, u) / (cssq)));
end

#! Equilibrium frequency distribution for incompressible Newtonian flow
#!
#! \param lat D2Q4 lattice
#! \param rho Macroscopic density at lattice site
#! \param u Macroscopic flow at lattice site
#! \return Equilibrium frequency
function feq_incomp(lat::LatticeD2Q4, rho::FloatingPoint, u::Vector{Float64},
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
#! \return Equilibrium frequency
function feq_incomp_HL(lat::LatticeD2Q4, rho::FloatingPoint,
                       rho_0::FloatingPoint, u::Vector{Float64})
  const cssq = 1/2;
  const ckdotu = dot(lat.c[:,k], u);

  return lat.w[k] * (rho + rho_0 * 3 * ckdotu);
end