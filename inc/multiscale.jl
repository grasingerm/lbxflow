const __multiscale_root__ = dirname(@__FILE__);
require(abspath(joinpath(__multiscale_root__, "lattice.jl")));
require(abspath(joinpath(__multiscale_root__, "numerics.jl")));

# TODO: storing dx and dt in both the lattice and multiscale maps seems redundant
# TODO: probably not necessary storing dx or dt AT ALL
type MultiscaleMap
  dx::FloatingPoint;
  dt::FloatingPoint;
  rho_0::FloatingPoint;
  omega::Array{Float64,2};
  u::Array{Float64,3};
  rho::Array{Float64,2};

  MultiscaleMap(dx::FloatingPoint, dt::FloatingPoint, rho_0::FloatingPoint,
    omega::Array{Float64,2}, u::Array{Float64,3}, rho::Array{Float64,2}) =
    new(dx, dt, rho_0, omega, u, rho);

  MultiscaleMap(nu::FloatingPoint, dx::FloatingPoint, dt::FloatingPoint,
    ni::Unsigned, nj::Unsigned) =
    new(dx, dt, 1.0, fill(1.0/(3 * nu * dt / (dx*dx) + 0.5), (ni, nj)),
      zeros(Float64, (ni, nj, 2)), zeros(Float64, (ni, nj)));

  function MultiscaleMap(nu::FloatingPoint, lat::Lattice, rho::FloatingPoint)
    const dx = lat.dx;
    const dt = lat.dt;
    const ni, nj = size(lat.f);

    new(dx, dt, rho, fill(1.0/(3 * nu * dt / (dx*dx) + 0.5), (ni, nj)),
      zeros(Float64, (ni, nj, 2)), fill(rho, (ni, nj)));
  end

  function MultiscaleMap(nu::FloatingPoint, lat::Lattice)
    const dx = lat.dx;
    const dt = lat.dt;
    const ni, nj = size(lat.f);

    new(dx, dt, 1.0, fill(1.0/(3 * nu * dt / (dx*dx) + 0.5), (ni, nj)),
      zeros(Float64, (ni, nj, 2)), zeros(Float64, (ni, nj)));
  end

  function MultiscaleMap(msm::MultiscaleMap)
    new(msm.dx, msm.dt, msm.rho_0, copy(msm.omega), copy(msm.u), copy(msm.rho));
  end
end

#! Map particle distribution frequencies to macroscopic variables
function map_to_macro!(lat::Lattice, msm::MultiscaleMap)
  const ni, nj = size(lat.f);
  const nk = length(lat.w);

  for i=1:ni, j=1:nj
    msm.rho[i,j] = 0;

    for k=1:nk
      msm.rho[i,j] += lat.f[i,j,k];
    end

    msm.u[i,j,1] = 0;
    msm.u[i,j,2] = 0;

    for k=1:nk
      msm.u[i,j,1] += lat.f[i,j,k] * lat.c[k,1];
      msm.u[i,j,2] += lat.f[i,j,k] * lat.c[k,2];
    end

    for a=1:2
      msm.u[i,j,a] = msm.u[i,j,a] / msm.rho[i,j];
    end
  end

end

#! Reynold's number
#!
#! \param u Magnitude of macroscopic flow
#! \param L Characteristic length of flow
#! \param nu Kinematic viscosity
#! \return Reyond's number
function reynolds(u::Number, L::Number, nu::Number)
  return u * L / nu;
end

#! Calculate the magnitude of velocity at each lattice node
#!
#! \param msm Multiscale map
#! \return Velocity magnitudes
function u_mag(msm::MultiscaleMap)
  const ni, nj = size(msm.u);
  u_mag_res = Array(Float64, (ni, nj));

  for i = 1:ni, j = 1:nj
    u = msm.u[i,j,1];
    v = msm.u[i,j,2];
    u_mag_res[i,j] = sqrt(u*u + v*v);
  end

  return u_mag_res;
end

#! Calculate local strain rate tensor
#!
#! \param rho Local density
#! \param f_neq Non-equilibrium distribution (f_neq1, f_neq2, ..., f_neq9)
#! \param c Lattice vectors
#! \param c_sq Lattice speed of sound squared
#! \param dt Change in time
#! \param M Transformation matrix to map f from velocity space to momentum space
#! \param S (Sparse) diagonal relaxation matrix
#! \return D Strain rate matrix
function strain_rate_tensor(rho::Float64, f_neq::Array{Float64, 1},
  c::Array{Float64, 2}, c_ssq::Float64, dt::Float64, M::Array{Float64,2},
  S::SparseMatrixCSC)

  D = zeros(Float64, (2, 2)); #!< Heuristic, 2D so 2x2
  const ni = length(f_neq);

  const MiSM = inv(M) * S * M;

  for alpha=1:2, beta=1:2
    sum = 0;
    for i=1:ni
      MiSMsum = 0;
      for j=1:ni
        MiSMsum += MiSM[i, j] * f_neq[j];
      end
      sum += c[i, alpha] * c[i, beta] * MiSMsum;
    end
    D[alpha, beta] = -1.0 / (2.0 * rho * c_ssq * dt) * sum;
  end

  return D;
end

#! Calculate local strain rate tensor
#!
#! \param rho Local density
#! \param f_neq Non-equilibrium distribution (f_neq1, f_neq2, ..., f_neq9)
#! \param c Lattice vectors
#! \param c_sq Lattice speed of sound squared
#! \param dt Change in time
#! \param M Transformation matrix to map f from velocity space to momentum space
#! \param S (Sparse) diagonal relaxation matrix
#! \return D Strain rate tensor
function strain_rate_tensora(rho::Float64, f_neq::Array{Float64, 1},
  c::Array{Float64, 2}, c_ssq::Float64, dt::Float64, M::Array{Float64,2},
  S::SparseMatrixCSC)

  D = zeros(Float64, (2, 2)); #!< Heuristic, 2D so 2x2
  const ni = length(f_neq);
  const MiSM = inv(M) * S * M;
  const mult = -1.0 / (2.0 * rho * c_ssq * dt);

  for a = 1:2, b = 1:2
    outsum = 0;
    for i = 1:ni
      insum = 0;
      for j = 1:ni
        insum += MiSM[i,j] * f_neq[j];
      end
      outsum += c[i,a] * c[i,b] * insum;
    end
    D[a,b] = mult * outsum;
  end

  return D;
end

#! Calculate divergence of the strain rate tensor
function div_strain_rate(D::Array{Float64, 2}, c::Array{Float64, 2})
  const ni, = size(c);
  divD = zeros(ni);

  for beta = 1:2
    for i = 1:ni, alpha = 1:2
      divD[alpha, beta] += c[i, alpha] * D[alpha, beta];
    end
  end

  return divD / 6.0;
end

#! Calculate strain rate from strain rate tensor
macro strain_rate(D)
  return :(sqrt(2.0 * ddot($D, $D)));
end

#! Calculate apparent viscosity of a Bingham plastic using Papanstasiou's model
macro papanstasiou_mu(mu_p, tau_o, m, gamma)
  return :($mu_p + $tau_o / $gamma * (1 - exp(-$m * abs($gamma))));
end

#! Calculate relaxation time
macro relax_t(mu, rho, dx, dt)
  return :(3 * $mu / $rho * (($dt * $dt) / ($dx * $dx)) + 0.5 * $dt);
end

#! Calculate collision frequency
macro omega(mu, rho, dx, dt)
	return :(1.0 / (3.0 * $mu / $rho * ($dt)/($dx * $dx) + 0.5));
end

#! Viscosity from collision frequency
macro mu(omega, rho, c_ssq)
  return :((1.0/$omega - 0.5) * $rho * $c_ssq);
end

#! Macroscopic pressure from rho
macro p(rho, c_ssq)
  return :($rho * $c_ssq);
end

#! Macroscopic rho from pressure
macro rho(p, c_ssq)
  return :($p / $c_ssq);
end

function dumpsf(w::IOStream, msm::MultiscaleMap)
  ni, nj = size(msm.rho);

  write(w, "msm:\n");
  write(w, "  ni: $ni\n");
  write(w, "  nj: $nj\n");
  write(w, "  dx: ", msm.dx, "\n");
  write(w, "  dt: ", msm.dt, "\n");
  write(w, "  rho_0: ", msm.rho_0, "\n");

  write(w, "  omega: ");
  for omegaij in msm.omega
    write(w, omegaij);
  end
  write(w, "\n");

  write(w, "  u: ");
  for uij in msm.u
    write(w, uij);
  end
  write(w, "\n");

  write(w, "  rho: ");
  for rhoij in msm.rho
    write(w, rhoij);
  end
  write(w, "\n\n");
end
