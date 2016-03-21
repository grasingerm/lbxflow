# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Calculate strain rate from strain rate tensor
macro strain_rate(D)
  return :(sqrt(2.0 * ddot($D, $D)));
end

# TODO: create a type `ConstitutiveModel`
#! Calculate apparent viscosity of a Bingham plastic using Papanstasiou's model
macro mu_papanstasiou(mu_p, tau_y, m, gamma)
  return :($mu_p + $tau_y / $gamma * (1 - exp(-$m * abs($gamma))));
end

#! Calculate relaxation time
macro relax_t(nu, cssq, dt)
  return :($nu / ($cssq*$dt) + 0.5);
end

#! Calculate collision frequency
macro omega(nu, cssq, dt)
	return :(1.0 / @relax_t($nu, $cssq, $dt));
end

#! Viscosity from collision frequency
macro nu(omega, cssq, dt)
  return :($cssq*$dt * (1.0/$omega - 0.5));
end

#! Macroscopic pressure from rho
macro p(rho, cssq)
  return :($rho * $cssq);
end

#! Macroscopic rho from pressure
macro rho(p, cssq)
  return :($p / $cssq);
end

#! Multiscale map for resolving macroscopic parameters
immutable MultiscaleMap
  rho_0::AbstractFloat;
  omega::Matrix{Float64};
  u::Array{Float64,3};
  rho::Matrix{Float64};

  function MultiscaleMap(nu::AbstractFloat, lat::Lattice, rho::AbstractFloat = 1.0)
    const ni, nj, = (size(lat.f, 2), size(lat.f, 3));

    new(rho, fill(@omega(nu, lat.cssq, lat.dt), (ni, nj)),
        zeros(Float64, (2, ni, nj)), fill(rho, (ni, nj)));
  end
  
  MultiscaleMap(rho_0::AbstractFloat, omega::Matrix{Float64},
                u::Array{Float64,3}, rho::Matrix{Float64}) =
    new(rho_0, omega, u, rho);

  MultiscaleMap(msm::MultiscaleMap) =
    new(msm.rho_0, copy(msm.omega), copy(msm.u), copy(msm.rho));
end

#! Map particle distribution frequencies to macroscopic variables
function map_to_macro!(lat::Lattice, msm::MultiscaleMap)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  const nk = lat.n;

  for j=1:nj, i=1:ni
    msm.rho[i,j] = 0;
    msm.u[1,i,j] = 0;
    msm.u[2,i,j] = 0;

    for k=1:nk
      msm.rho[i,j] += lat.f[k,i,j];
      msm.u[1,i,j] += lat.f[k,i,j] * lat.c[1,k];
      msm.u[2,i,j] += lat.f[k,i,j] * lat.c[2,k];
    end

    for a=1:2
      msm.u[a,i,j] = msm.u[a,i,j] / msm.rho[i,j];
    end
  end

end

#! Reynold's number
#!
#! \param u Magnitude of macroscopic flow
#! \param L Characteristic length of flow
#! \param nu Kinematic viscosity
#! \return Reyond's number
function reynolds(u::Real, L::Real, nu::Real)
  return u * L / nu;
end

#! Bingham number
#!
#! \param u     Magnitude of macroscopic flow
#! \param L     Characteristic length of flow
#! \param mu_p  Plastic viscosity
#! \param tau_y Yield stress
#! \return      Bingham number of the flow
function bingham_number(u::Real, L::Real, mu_p::Real, tau_y::Real)
  return tau_y * L / (mu_p * u);
end

#! Calculate the magnitude of velocity at each lattice node
#!
#! \param msm Multiscale map
#! \return Velocity magnitudes
function u_mag(msm::MultiscaleMap)
  const ni, nj = size(msm.u, 2), size(msm.u, 3);
  u_mag_res = Array(Float64, (ni, nj));

  for j = 1:nj, i = 1:ni
    u = msm.u[1,i,j];
    v = msm.u[2,i,j];
    u_mag_res[i,j] = sqrt(u*u + v*v);
  end

  return u_mag_res;
end

#! Calculate local strain rate tensor
#!
#! \param lat Lattice
#! \param rho Local density
#! \param f_neq Non-equilibrium distribution (f_neq1, f_neq2, ..., f_neq9)
#! \param omega Collision frequency
#! \return D Strain rate matrix
function strain_rate_tensor(lat::Lattice, rho::Number, fneq::Vector{Float64},
                            omega::Number)

  D = zeros(Float64, (2, 2)); #!< Heuristic, 2D so 2x2
  const ni = length(fneq);

  for alpha=1:2, beta=1:2
    sum = 0;
    for i=1:ni
      sum += lat.c[alpha,i] * lat.c[beta,i] * fneq[i];
    end
    D[alpha, beta] = -omega / (2.0 * rho * lat.cssq * lat.dt) * sum;
  end

  return D;
end


#! Calculate local strain rate tensor
#!
#! \param lat Lattice
#! \param rho Local density
#! \param f_neq Non-equilibrium distribution (f_neq1, f_neq2, ..., f_neq9)
#! \param M Transformation matrix to map f from velocity space to momentum space
#! \param S (Sparse) diagonal relaxation matrix
#! \return D Strain rate matrix
function strain_rate_tensor(lat::Lattice, rho::Number, fneq::Vector{Float64},
                            M::Matrix{Float64}, iM::Matrix{Float64},
                            S::SparseMatrixCSC)

  D = zeros(Float64, (2, 2)); #!< Heuristic, 2D so 2x2
  const ni = length(fneq);
  const MiSM = iM * S * M;

  for alpha=1:2, beta=1:2
    sum = 0;
    for i=1:ni
      MiSMsum = 0;
      for j=1:ni
        MiSMsum += MiSM[i,j] * fneq[j];
      end
      sum += lat.c[alpha,i] * lat.c[beta,i] * MiSMsum;
    end
    D[alpha, beta] = -1.0 / (2.0 * rho * lat.cssq * lat.dt) * sum;
  end

  return D;
end

#! Calculate divergence of the strain rate tensor
function div_strain_rate(D::Matrix{Float64}, c::Matrix{Float64})
  const ni, = size(c);
  divD = zeros(ni);

  for beta = 1:2
    for i = 1:ni, alpha = 1:2
      divD[alpha, beta] += c[i, alpha] * D[alpha, beta];
    end
  end

  return divD / 6.0;
end
