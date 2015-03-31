const __collision_root__ = dirname(@__FILE__);
require(abspath(joinpath(__collision_root__, "debug.jl")));
require(abspath(joinpath(__collision_root__, "lattice.jl")));
require(abspath(joinpath(__collision_root__, "multiscale.jl")));
require(abspath(joinpath(__collision_root__, "numerics.jl")));

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
function srt_col_f! (lat::Lattice, msm::MultiscaleMap)

  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  for i=1:ni, j=1:nj, k=1:9
    rhoij = msm.rho[i,j];
    uij = vec(msm.u[i,j,:]);
    f_eq = incomp_f_eq(rhoij, lat.w[k], c_ssq, vec(lat.c[k,:]), uij);
    lat.f[i,j,k] = msm.omega[i,j] * f_eq + (1.0 - msm.omega[i,j]) * lat.f[i,j,k];
  end

end

#! Equilibrium frequency distribution for incompressible Newtonian flow
#!
#! \param rho Density at lattice site
#! \param w Weight for lattice direction
#! \param c_ssq Lattice speed of sound squared
#! \param c_k Vector for lattice direction
#! \param u Macroscopic flow at lattice site
#! \return Equilibrium frequency
function incomp_f_eq(rho::Float64, w::Float64, c_ssq::Float64,
  c_k::Array{Float64, 1}, u::Array{Float64, 1})

  const ckdotu = dot(c_k, u);

  return rho * w * (1.0 + ckdotu/(c_ssq) + 0.5*(ckdotu*ckdotu)/(c_ssq*c_ssq)
    - 0.5 * dot(u, u) / (c_ssq));
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
#! \param M Transmation matrix to map f from velocity space to momentum space
#! \param S (Sparse) diagonal relaxation matrix
function mrt_col_f! (lat::Lattice, msm::MultiscaleMap, M::Array{Float64,2},
  S::SparseMatrixCSC)

  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  for i=1:ni, j=1:nj
    rhoij = msm.rho[i, j];
    uij = vec(msm.u[i, j, :]);

    for k=1:9
      f_eq[k] = incomp_f_eq(rhoij, lat.w[k], c_ssq, vec(lat.c[k,:]), uij);
    end

    f = vec(lat.f[i,j,:]);
    m = M * f;
    m_eq = M * f_eq;

    lat.f[i,j,:] = f - iM * S * (m - m_eq); # perform collision
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
#! \param M Transmation matrix to map f from velocity space to momentum space
#! \param S Function for creating a sparse diagonal relaxation matrix
function mrt_col_f! (lat::Lattice, msm::MultiscaleMap, M::Array{Float64,2},
  S::Function)

  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  for i=1:ni, j=1:nj
    rhoij = msm.rho[i, j];
    uij = vec(msm.u[i, j, :]);

    for k=1:9
      f_eq[k] = incomp_f_eq(rhoij, lat.w[k], c_ssq, vec(lat.c[k,:]), uij);
    end

    f = vec(lat.f[i,j,:]);
    m = M * f;
    m_eq = M * f_eq;

    lat.f[i,j,:] = f - iM * S(lat, msm) * (m - m_eq); # perform collision
  end
end

#! Initializes the default multiple relaxation time transformation matrix
macro DEFAULT_MRT_M()
  return :( [1.0    1.0    1.0    1.0    1.0    1.0    1.0    1.0    1.0;
            -4.0   -1.0   -1.0   -1.0   -1.0    2.0    2.0    2.0    2.0;
             4.0   -2.0   -2.0   -2.0   -2.0    1.0    1.0    1.0    1.0;
             0.0    1.0    0.0   -1.0    0.0    1.0   -1.0   -1.0    1.0;
             0.0   -2.0    0.0    2.0    0.0    1.0   -1.0   -1.0    1.0;
             0.0    0.0    1.0    0.0   -1.0    1.0    1.0   -1.0   -1.0;
             0.0    0.0   -2.0    0.0    2.0    1.0    1.0   -1.0   -1.0;
             0.0    1.0   -1.0    1.0   -1.0    0.0    0.0    0.0    0.0;
             0.0    0.0    0.0    0.0    0.0    1.0   -1.0    1.0   -1.0]
          );
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
#! \param S (Sparse) diagonal relaxation matrix
function mrt_col_f! (lat::Lattice, msm::MultiscaleMap, S::SparseMatrixCSC)
  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  for i=1:ni, j=1:nj
		rhoij = msm.rho[i, j];
		uij = vec(msm.u[i, j, :]);

    for k=1:9
      f_eq[k] = incomp_f_eq(rhoij, lat.w[k], c_ssq, vec(lat.c[k,:]), uij);
    end

    f = vec(lat.f[i,j,:]);
    m = M * f;
    m_eq = M * f_eq;

    lat.f[i,j,:] = f - iM * S * (m - m_eq); # perform collision
  end
end

#! Vikhansky relaxation coefficient for s77, s88
#!
#! \param mu Dynamic viscosity
#! \param rho Local density
#! \param c_ssq Lattice speed of sound squared
#! \param dt Change in time
#! \return Vikhansky relaxation coefficient for s77 and s88
macro viks_8(mu, rho, c_ssq, dt)
  return :(1.0/($mu / ($rho * $c_ssq * $dt) + 0.5));
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \param mu_p Plastic viscosity
#! \param tau_y Yield stress
#! \param m Stress growth exponent
#! \param max_iters Maximum iterations for determining apparent viscosity
#! \param tol Tolerance for apparent viscosity convergence
function mrt_bingham_col_f! (lat::Lattice, msm::MultiscaleMap, S::Function,
  mu_p::Number, tau_y::Number, m::Number, max_iters::Int, tol::FloatingPoint)
  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  for i=1:ni, j=1:nj
    rhoij = msm.rho[i, j];
    uij = vec(msm.u[i, j, :]);

    for k=1:9
      f_eq[k] = incomp_f_eq(rhoij, lat.w[k], c_ssq, vec(lat.c[k,:]), uij);
    end

    # initialize density, viscosity, and relaxation matrix at node i,j
		muij = @mu(msm.omega[i,j], rhoij, c_ssq);
    Sij = S(muij, rhoij, c_ssq, lat.dt);

    f = vec(lat.f[i,j,:]);
    mij = M * f;
    mij_eq = M * f_eq;
    f_neq = f - f_eq;
    muo = muij;

    #=
    println("rhoij = ", rhoij);
    println("uij = ", uij);
    println("muij = ", muij);
    println("Sij = ", Sij);
    println("f = ", f);
    println("mij = ", mij);
    println("mij_eq = ", mij_eq);
    println("f_neq = ", f_neq);
    println("muo = ", muo);
    =#

    # iteratively determine mu
    iters = 0;
    mu_prev = muo;

    while true
      iters += 1;

      D = strain_rate_tensor(msm.rho[i,j], f_neq, lat.c, c_ssq, lat.dt, M, Sij);
      gamma = @strain_rate(D);

      # update relaxation matrix
      muij = mu_p + tau_y / gamma * (1 - exp(-m * abs(gamma)));
      s_8 = @viks_8(muij, rhoij, c_ssq, lat.dt);
      Sij[8,8] = s_8;
      Sij[9,9] = s_8;

      # check for convergence
      if abs(mu_prev - muij) / muo <= tol || iters > max_iters
        break;
      end

      #=println("$iters, ", abs(mu_prev - muij) / muo);
      println("D = ", D);
      println("gamma = ", gamma);
      println("muij = ", muij);
      println("tau = ", @relax_t(muij, rhoij, lat.dx, lat.dt));
      println("Sij = ", Sij);
      println("Enter to continue...");
      readline(STDIN);=#

      mu_prev = muij;
    end

    #=@mdebug(
      @relax_t(muij, rhoij, lat.dx, lat.dt) > 0.5 &&
      @relax_t(muij, rhoij, lat.dx, lat.dt) <= 8.0,
      "Warning: relaxation time should be between 0.5 and 8.0"
    );=#

    lat.f[i,j,:] = f - iM * Sij * (mij - mij_eq); # perform collision
    # update collision frequency matrix
    msm.omega[i,j] = @omega(muij, rhoij, lat.dx, lat.dt);
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \param mu_p Plastic viscosity
#! \param tau_y Yield stress
#! \param m Stress growth exponent
#! \param max_iters Maximum iterations for determining apparent viscosity
#! \param tol Tolerance for apparent viscosity convergence
#! \param tau_min Minimum relaxation time
#! \param tau_max Maximum relaxation time
function mrt_bingham_col_f! (lat::Lattice, msm::MultiscaleMap, S::Function,
  mu_p::Number, tau_y::Number, m::Number, max_iters::Int, tol::FloatingPoint,
  tau_min::Number, tau_max::Number)

  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  for i=1:ni, j=1:nj
    rhoij = msm.rho[i, j];
    uij = vec(msm.u[i, j, :]);

    for k=1:9
      f_eq[k] = incomp_f_eq(rhoij, lat.w[k], c_ssq, vec(lat.c[k,:]), uij);
    end

    # initialize density, viscosity, and relaxation matrix at node i,j
    muij = @mu(msm.omega[i,j], rhoij, c_ssq);
    Sij = S(muij, rhoij, c_ssq, lat.dt);

    f = vec(lat.f[i,j,:]);
    mij = M * f;
    mij_eq = M * f_eq;
    f_neq = f - f_eq;
    muo = muij;

    # iteratively determine mu
    iters = 0;
    mu_prev = muo;

    while true
      iters += 1;

      D = strain_rate_tensor(msm.rho[i,j], f_neq, lat.c, c_ssq, lat.dt, M, Sij);
      gamma = @strain_rate(D);

      # update relaxation matrix
      muij = mu_p + tau_y / gamma * (1 - exp(-m * abs(gamma)));
      s_8 = @viks_8(muij, rhoij, c_ssq, lat.dt);
      Sij[8,8] = s_8;
      Sij[9,9] = s_8;

      # check for convergence
      if abs(mu_prev - muij) / muo <= tol || iters > max_iters
        break;
      end

      mu_prev = muij;
    end

    # correct relaxation time based on stability limits
    tau = @relax_t(muij, rhoij, lat.dx, lat.dt);
    if tau > tau_max || tau < tau_min
      if tau > tau_max
        @mdebug("Relaxation time, $tau > $tau_max, maximum allowed. Setting
          relaxation time to $tau_max for numerical stability.")
        tau = tau_max;
      elseif tau < tau_min
        @mdebug("Relaxation time, $tau < $tau_min, minimum allowed. Setting
          relaxation time to $tau_min for numerical stability.")
        tau = tau_min;
      end

      muij = @mu((lat.dx / tau), rhoij, c_ssq);

      # update relaxation matrix
      s_8 = @viks_8(muij, rhoij, c_ssq, lat.dt);
      Sij[8,8] = s_8;
      Sij[9,9] = s_8;
    end

    lat.f[i,j,:] = f - iM * Sij * (mij - mij_eq); # perform collision

    # update collision frequency matrix
    msm.omega[i,j] = @omega(muij, rhoij, lat.dx, lat.dt);
  end
end

#! Vikhansky relaxation matrix
#!
#! \param mu Dynamic viscosity
#! \param rho Local density
#! \param c_ssq Lattice speed of sound squared
#! \param dt Change in time
#! \return Vikhansky relaxation matix
function vikhansky_relax_matrix(mu::Number, rho::Number, c_ssq::Number,
	dt::Number)

	const s_8 = @viks_8(mu, rho, c_ssq, dt);
	return spdiagm([0.0; 1.1; 1.1; 0.0; 1.1; 0.0; 1.1; s_8; s_8]);

end

