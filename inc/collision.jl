const __collision_root__ = dirname(@__FILE__);
require(abspath(joinpath(__collision_root__, "debug.jl")));
require(abspath(joinpath(__collision_root__, "lattice.jl")));
require(abspath(joinpath(__collision_root__, "multiscale.jl")));
require(abspath(joinpath(__collision_root__, "numerics.jl")));

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
#! \param rho Density at lattice site
#! \param w Weight for lattice direction
#! \param c_ssq Lattice speed of sound squared
#! \param c_k Vector for lattice direction
#! \param u Macroscopic flow at lattice site
#! \return Equilibrium frequency
function feq_incomp_HL(lat::LatticeD2Q9, rho::FloatingPoint,
                       rho_0::FloatingPoint, u::Vector{Float64})
  const cssq = 1/3;
  const ckdotu = dot(lat.c[:,k], u);

  return w * (rho + rho_0 * (ckdotu/(cssq) + 0.5*(ckdotu*ckdotu)/(cssq*cssq)
              - 0.5 * dot(u, u) / (cssq)));
end

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
function srt_col_f! (lat::Lattice, msm::MultiscaleMap, bounds::Array{Int64,2})

  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  const nbounds, = size(bounds);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];
    for i = i_min:i_max, j = j_min:j_max, k = 1:9
      rhoij = msm.rho[i,j];
      uij = vec(msm.u[i,j,:]);
      f_eq = incomp_f_eq(rhoij, lat.w[k], c_ssq, vec(lat.c[k,:]), uij);
      lat.f[i,j,k] = (msm.omega[i,j] * f_eq + (1.0 - msm.omega[i,j])
                        * lat.f[i,j,k]);
    end
  end

end

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
#! \param f Body force vector
function srt_col_f! (lat::Lattice, msm::MultiscaleMap, f::Array{Float64,1},
  bounds::Array{Int64,2})

  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  const nbounds, = size(bounds);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];

    for i = i_min:i_max, j = j_min:j_max
      rhoij = msm.rho[i,j];
      uij = vec(msm.u[i,j,:]);
      omegaij = msm.omega[i,j];

      for k=1:9
        ck = vec(lat.c[k,:]);
        wk = lat.w[k];

        f_eq = incomp_f_eq(rhoij, wk, c_ssq, ck, uij);

        # body force is incorporated with w*dt/c_ssq * dot(f,c)
        lat.f[i,j,k] = (omegaij * f_eq
                        + (1.0 - omegaij) * lat.f[i,j,k]
                        + wk * lat.dt / c_ssq * dot(f, ck));
      end
    end
  end

end

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
#! \param f Body force vector
function srt_guo_col_f! (lat::Lattice, msm::MultiscaleMap, f::Array{Float64,1},
  bounds::Array{Int64,2})

  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);
  const nbounds, = size(bounds);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];
    for i = i_min:i_max, j = j_min:j_max
      rhoij = msm.rho[i,j];

      uij = vec(msm.u[i,j,:]) + lat.dt / 2.0 * f;
      omegaij = msm.omega[i,j];

      for k=1:9
        ck = vec(lat.c[k,:]);
        wk = lat.w[k];

        f_eq = incomp_f_eq(rhoij, wk, c_ssq, ck, uij);

        fdlk = (1 - 0.5 * omegaij) * wk * dot(((ck - uij) / c_ssq +
                  dot(ck, uij) / (c_ssq * c_ssq) * ck), f);

        lat.f[i,j,k] = (omegaij * f_eq
                        + (1.0 - omegaij) * lat.f[i,j,k] + fdlk);
      end
    end
  end

end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
#! \param M Transmation matrix to map f from velocity space to momentum space
#! \param S (Sparse) diagonal relaxation matrix
function mrt_col_f! (lat::Lattice, msm::MultiscaleMap, M::Array{Float64,2},
  S::SparseMatrixCSC, bounds::Array{Int64,2})

  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);

  const nbounds, = size(bounds);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];
    for i = i_min:i_max, j = j_min:j_max
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
function mrt_col_f! (lat::Lattice, msm::MultiscaleMap, S::SparseMatrixCSC,
  bounds::Array{Int64,2})

  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  const nbounds, = size(bounds);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];
    for i = i_min:i_max, j = j_min:j_max
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
end

# TODO: this actually comes from Fallah 2012 and a few other studies

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
#! \param gamma_min Minimum strain rate to use in apparent viscosity calculation
#! \param bounds 2D array, each row is i_min, i_max, j_min, j_max
#! \param relax Relaxation factor for updating apparent viscosity
function mrt_bingham_explicit_col_f! (lat::Lattice, msm::MultiscaleMap,
  S::Function, mu_p::Number, tau_y::Number, m::Number, gamma_min::FloatingPoint,
  bounds::Array{Int64,2}, relax::Number = 1.0)

  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  const nbounds, = size(bounds);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];
    for i = i_min:i_max, j = j_min:j_max
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

      D = strain_rate_tensor(msm.rho[i,j], f_neq, lat.c, c_ssq, lat.dt, M, Sij);
      gamma = @strain_rate(D);

      # update relaxation matrix
      gamma = gamma < gamma_min ? gamma_min : gamma;
      muij = (
               (1 - relax) * muij
                + relax * (mu_p + tau_y / gamma * (1 - exp(-m * gamma)))
             );
      s_8 = @viks_8(muij, rhoij, c_ssq, lat.dt);
      Sij[8,8] = s_8;
      Sij[9,9] = s_8;

      lat.f[i,j,:] = f - iM * Sij * (mij - mij_eq); # perform collision

      # update collision frequency matrix
      msm.omega[i,j] = @omega(muij, rhoij, lat.dx, lat.dt);
    end
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
#! \param gamma_min Minimum strain rate to use in apparent viscosity calculation
#! \param f Body force vector
#! \param bounds 2D array, each row is i_min, i_max, j_min, j_max
#! \param relax Relaxation factor for updating apparent viscosity
function mrt_bingham_explicit_col_f! (lat::Lattice, msm::MultiscaleMap,
  S::Function, mu_p::Number, tau_y::Number, m::Number, gamma_min::FloatingPoint,
  f::Array{Float64,1}, bounds::Array{Int64,2}, relax::Number = 1.0)

  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  const nbounds, = size(bounds);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];
    for i = i_min:i_max, j = j_min:j_max
      rhoij = msm.rho[i, j];
      uij = vec(msm.u[i,j,:]) + lat.dt / 2.0 * f;

      for k=1:9
        f_eq[k] = incomp_f_eq(rhoij, lat.w[k], c_ssq, vec(lat.c[k,:]), uij);
      end

      # initialize density, viscosity, and relaxation matrix at node i,j
      muij = @mu(msm.omega[i,j], rhoij, c_ssq);
      Sij = S(muij, rhoij, c_ssq, lat.dt);

      fij = vec(lat.f[i,j,:]);
      mij = M * fij;
      mij_eq = M * f_eq;
      f_neq = fij - f_eq;
      muo = muij;

      D = strain_rate_tensor(msm.rho[i,j], f_neq, lat.c, c_ssq, lat.dt, M, Sij);
      gamma = @strain_rate(D);

      # update relaxation matrix
      gamma = gamma < gamma_min ? gamma_min : gamma;
      muij = (
               (1 - relax) * muij
                + relax * (mu_p + tau_y / gamma * (1 - exp(-m * gamma)))
             );
               
      s_8 = @viks_8(muij, rhoij, c_ssq, lat.dt);
      Sij[8,8] = s_8;
      Sij[9,9] = s_8;

      omegaij = @omega(muij, rhoij, lat.dx, lat.dt);

      fdl = Array(Float64, 9);
      for k=1:9
        ck = vec(lat.c[k,:]);
        fdl[k] = (1 - 0.5 * omegaij) * lat.w[k] * dot(((ck - uij) / c_ssq +
                  dot(ck, uij) / (c_ssq * c_ssq) * ck), f);
      end

      lat.f[i,j,:] = fij - iM * Sij * (mij - mij_eq) + fdl; # perform collision

      # update collision frequency matrix
      msm.omega[i,j] = omegaij;
    end
  end
end


#TODO: reconsider order of parameters... come up with a convenction
# phys models, const/material params, Functions, knobs

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
#! \param gamma_min Minimum strain rate to use in apparent viscosity calculation
#! \param bounds 2D array, each row is i_min, i_max, j_min, j_max
#! \param relax Relaxation factor for updating apparent viscosity
function mrt_bingham_implicit_col_f! (lat::Lattice, msm::MultiscaleMap,
  S::Function, mu_p::Number, tau_y::Number, m::Number, max_iters::Int,
  tol::FloatingPoint, gamma_min::FloatingPoint, bounds::Array{Int64,2},
  relax::Number = 1.0)

  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  const nbounds, = size(bounds);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];
    for i = i_min:i_max, j = j_min:j_max
      rhoij = msm.rho[i, j];
      uij = vec(msm.u[i, j, :]);

      for k=1:9
        f_eq[k] = incomp_f_eq(rhoij, lat.w[k], c_ssq, vec(lat.c[k,:]), uij);
      end

      # initialize density, viscosity, and relaxation matrix at node i,j
      muij = @mu(msm.omega[i,j], rhoij, c_ssq);
      Sij = S(muij, rhoij, c_ssq, lat.dt);

      fij = vec(lat.f[i,j,:]);
      mij = M * fij;
      mij_eq = M * f_eq;
      f_neq = fij - f_eq;
      muo = muij;

      # iteratively determine mu
      iters = 0;
      mu_prev = muo;

      while true
        iters += 1;

        D = strain_rate_tensor(msm.rho[i,j], f_neq, lat.c, c_ssq, lat.dt, M, Sij);
        gamma = @strain_rate(D);

        # update relaxation matrix
        gamma = (gamma < gamma_min) ? gamma_min : gamma;

        muij = (
                  (1 - relax) * muij
                  + relax * (mu_p + tau_y / gamma * (1 - exp(-m * gamma)))
               );

        s_8 = @viks_8(muij, rhoij, c_ssq, lat.dt);
        Sij[8,8] = s_8;
        Sij[9,9] = s_8;

        # check for convergence
        if abs(mu_prev - muij) / muo <= tol
          break;
        end

        if iters > max_iters
          break;
        end

        mu_prev = muij;
      end

      lat.f[i,j,:] = fij - iM * Sij * (mij - mij_eq); # perform collision

      # update collision frequency matrix
      msm.omega[i,j] = @omega(muij, rhoij, lat.dx, lat.dt);
    end
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
#! \param gamma_min Minimum strain rate to use in apparent viscosity calculation
#! \param f Body force
#! \param bounds 2D array, each row is i_min, i_max, j_min, j_max
#! \param relax Relaxation factor
function mrt_bingham_implicit_col_f! (lat::Lattice, msm::MultiscaleMap,
  S::Function, mu_p::Number, tau_y::Number, m::Number, max_iters::Int,
  tol::FloatingPoint, gamma_min::FloatingPoint, f::Array{Float64,1},
  bounds::Array{Int64,2}, relax::Number = 1.0)

  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  const nbounds, = size(bounds);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];
    for i = i_min:i_max, j = j_min:j_max
      rhoij = msm.rho[i, j];
      uij = vec(msm.u[i,j,:]) + lat.dt / 2.0 * f;

      for k=1:9
        f_eq[k] = incomp_f_eq(rhoij, lat.w[k], c_ssq, vec(lat.c[k,:]), uij);
      end

      # initialize density, viscosity, and relaxation matrix at node i,j
      muij = @mu(msm.omega[i,j], rhoij, c_ssq);
      Sij = S(muij, rhoij, c_ssq, lat.dt);

      fij = vec(lat.f[i,j,:]);
      mij = M * fij;
      mij_eq = M * f_eq;
      f_neq = fij - f_eq;
      muo = muij;

      # iteratively determine mu
      iters = 0;
      mu_prev = muo;

      while true
        iters += 1;

        D = strain_rate_tensor(msm.rho[i,j], f_neq, lat.c, c_ssq, lat.dt, M, Sij);
        gamma = @strain_rate(D);

        # update relaxation matrix
        gamma = (gamma < gamma_min) ? gamma_min : gamma;

        muij = (
                  (1 - relax) * muij
                  + relax * (mu_p + tau_y / gamma * (1 - exp(-m * gamma)))
               );

        s_8 = @viks_8(muij, rhoij, c_ssq, lat.dt);
        Sij[8,8] = s_8;
        Sij[9,9] = s_8;

        # check for convergence
        if abs(mu_prev - muij) / muo <= tol
          break;
        end

        if iters > max_iters
          break;
        end

        mu_prev = muij;
      end

      const omegaij = @omega(muij, rhoij, lat.dx, lat.dt);

      fdl = Array(Float64, 9);
      for k=1:9
        ck = vec(lat.c[k,:]);
        fdl[k] = (1 - 0.5 * omegaij) * lat.w[k] * dot(((ck - uij) / c_ssq +
                  dot(ck, uij) / (c_ssq * c_ssq) * ck), f);
      end

      lat.f[i,j,:] = fij - iM * Sij * (mij - mij_eq) + fdl; # perform collision

      # update collision frequency matrix
      msm.omega[i,j] = omegaij;
    end
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

#! Chen relaxation matrix
#!
#! \param omega Collision frequency
#! \return Chen relaxation matix
function chen_relax_matrix(omega::Number)
  return spdiagm([0.0; 1.1; 1.0; 0.0; 1.2; 0.0; 1.2; omega; omega]);
end

# =============================================================================
# ========================== not fully implemented functions! =================
# =============================================================================

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
#! \param gamma_min Minimum strain rate to use in apparent viscosity calculation
#! \param f Body force vector
function mrt_bingham_farnf_col_f! (lat::Lattice, msm::MultiscaleMap, S::Function,
  mu_p::Number, tau_y::Number, m::Number, max_iters::Int, tol::FloatingPoint,
  gamma_min::FloatingPoint, f::Array{Float64,1}, bounds::Array{Int64,2})

  error("Function not yet implemented!!!!");

  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  Gij = Array(Float64, 9);

  const nbounds, = size(bounds);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];
    for i = i_min:i_max, j = j_min:j_max
      rhoij = msm.rho[i, j];
      uij = vec(msm.u[i, j, :]);
      Fij = msm.F[i,j];

      for k=1:9
        wk = lat.w[k];
        ck = vec(lat.c[k,:]);

        f_eq[k] = incomp_f_eq_HL(rhoij, msm.rho_0, wk, c_ssq, ck, uij);
      end

      # initialize density, viscosity, and relaxation matrix at node i,j
      muij = @mu(msm.omega[i,j], rhoij, c_ssq);
      Sij = S(muij, rhoij, c_ssq, lat.dt);
      D_ext = zeros(2, 2);

      fij = vec(lat.f[i,j,:]);
      mij = M * fij;
      mij_eq = M * f_eq;
      f_neq = fij - f_eq;
      muo = muij;

      # iteratively determine mu
      iters = 0;
      mu_prev = muo;

      while true
        iters += 1;

        D = (strain_rate_tensor(msm.rho_0, f_neq, lat.c, c_ssq, lat.dt, M, Sij)
              + D_ext);
        divD = div_strain_rate(D, lat.c);
        gamma = @strain_rate(D);

        # update relaxation matrix
        gamma = gamma < gamma_min ? gamma_min : gamma;

        muij = mu_p + tau_y / gamma * (1 - exp(-m * gamma));

        s_8 = @viks_8(muij, rhoij, c_ssq, lat.dt);
        Sij[8,8] = s_8;
        Sij[9,9] = s_8;

        D_ext = -1.0 / (2 * c_ssq) * (Fij * uij' + uij * Fij');

        # check for convergence
        if abs(mu_prev - muij) / muo <= tol || iters > max_iters
          break;
        end

        mu_prev = muij;
      end

      lat.f[i,j,:] = fij - iM * Sij * (mij - mij_eq) + Gij; # perform collision

      # update collision frequency matrix
      msm.omega[i,j] = @omega(muij, rhoij, lat.dx, lat.dt);
    end
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
#! \param gamma_min Minimum strain rate to use in apparent viscosity calculation
#! \param f Body force vector
function mrt_bingham_farnf_col_f! (lat::Lattice, msm::MultiscaleMap, S::Function,
  mu_p::Number, tau_y::Number, m::Number, max_iters::Int, tol::FloatingPoint,
  gamma_min::FloatingPoint, f::Array{Float64,1}, bounds::Array{Int64,2})

  error("Function not yet implemented");

  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  Gij = Array(Float64, 9);

  const nbounds, = size(bounds);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];
    for i = i_min:i_max, j = j_min:j_max
      rhoij = msm.rho[i, j];
      uij = vec(msm.u[i, j, :]);
      Fij = msm.F[i,j];

      for k=1:9
        wk = lat.w[k];
        ck = vec(lat.c[k,:]);

        f_eq[k] = incomp_f_eq_HL(rhoij, msm.rho_0, wk, c_ssq, ck, uij);
      end

      # initialize density, viscosity, and relaxation matrix at node i,j
      muij = @mu(msm.omega[i,j], rhoij, c_ssq);
      Sij = S(muij, rhoij, c_ssq, lat.dt);
      D_ext = zeros(2, 2);

      fij = vec(lat.f[i,j,:]);
      mij = M * fij;
      mij_eq = M * f_eq;
      f_neq = fij - f_eq;
      muo = muij;

      # iteratively determine mu
      iters = 0;
      mu_prev = muo;

      while true
        iters += 1;

        D = (strain_rate_tensor(msm.rho_0, f_neq, lat.c, c_ssq, lat.dt, M, Sij)
              + D_ext);
        divD = div_strain_rate(D, lat.c);
        gamma = @strain_rate(D);

        # update relaxation matrix
        gamma = gamma < gamma_min ? gamma_min : gamma;

        muij = mu_p + tau_y / gamma * (1 - exp(-m * gamma));

        s_8 = @viks_8(muij, rhoij, c_ssq, lat.dt);
        Sij[8,8] = s_8;
        Sij[9,9] = s_8;

        D_ext = -1.0 / (2 * c_ssq) * (Fij * uij' + uij * Fij');

        # check for convergence
        if abs(mu_prev - muij) / muo <= tol || iters > max_iters
          break;
        end

        mu_prev = muij;
      end

      lat.f[i,j,:] = fij - iM * Sij * (mij - mij_eq) + Gij; # perform collision

      # update collision frequency matrix
      msm.omega[i,j] = @omega(muij, rhoij, lat.dx, lat.dt);
    end
  end
end

# =============================================================================
# =========================== DEPRECATED FUNCTIONS ============================
# =============================================================================

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
#! \param bounds 2D array, each row is i_min, i_max, j_min, j_max
#! \param relax Relaxation factor for updating apparent viscosity
function mrt_bingham_col_f! (lat::Lattice, msm::MultiscaleMap, S::Function,
  mu_p::Number, tau_y::Number, m::Number, max_iters::Int, tol::FloatingPoint,
  bounds::Array{Int64,2}, relax::Number = 1.0)

  global LBX_VERSION;
  if LBX_VERSION > 0.1
    error(string("`function mrt_bingham_col_f! (lat::Lattice,",
      " msm::MultiscaleMap, S::Function, mu_p::Number, tau_y::Number,",
      " m::Number, max_iters::Int, tol::FloatingPoint, bounds::Array{Int64,2},",
      " relax::Number = 1.0) is deprecated."));
  end

  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  const nbounds, = size(bounds);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[r,:];
    for i = i_min:i_max, j = j_min:j_max
      rhoij = msm.rho[i, j];
      uij = vec(msm.u[i,j,:]) + lat.dt / 2.0 * f;

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
        muij = ( 
                 (1 - relax) * muij +
                 relax * (mu_p + tau_y / gamma * (1 - exp(-m * abs(gamma))))
               );

        s_8 = @viks_8(muij, rhoij, c_ssq, lat.dt);
        Sij[8,8] = s_8;
        Sij[9,9] = s_8;

        # check for convergence
        if abs(mu_prev - muij) / muo <= tol || iters > max_iters
          break;
        end

        mu_prev = muij;
      end

      lat.f[i,j,:] = f - iM * Sij * (mij - mij_eq); # perform collision
      # update collision frequency matrix
      msm.omega[i,j] = @omega(muij, rhoij, lat.dx, lat.dt);
    end
  end
end
