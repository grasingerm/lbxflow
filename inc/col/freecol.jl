# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

const __freecol_root__ = dirname(@__FILE__);
require(abspath(joinpath(__freecol_root__, "..", "debug.jl")));
require(abspath(joinpath(__freecol_root__, "equilibrium.jl")));
require(abspath(joinpath(__freecol_root__, "..", "lattice.jl")));
require(abspath(joinpath(__freecol_root__, "mrt_matrices.jl")));
require(abspath(joinpath(__freecol_root__, "..", "multiscale.jl")));
require(abspath(joinpath(__freecol_root__, "..", "numerics.jl")));
require(abspath(joinpath(__freecol_root__, "..", "sim", "simtypes.jl")));

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param sim Simulation object
#! \param bounds Boundaries that define active parts of the lattice
function col_srt!(sim::FreeSurfSim, bounds::Matrix{Int64})
  lat = sim.lat;
  msm = sim.msm;
  const ni, nj = size(msm.rho);
  const nbounds = size(bounds, 2);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max, k = 1:lat.n
      if sim.tracker.state[i,j] == GAS; continue; end
      rhoij = msm.rho[i,j];
      uij = msm.u[:,i,j];
      feq = feq_incomp(lat, rhoij, uij, k);
      lat.f[k,i,j] = (msm.omega[i,j] * feq + (1.0 - msm.omega[i,j])
                        * lat.f[k,i,j]);
    end
  end
end

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param sim Simulation object
#! \param F Body force vector
#! \param bounds Boundaries that define active parts of the lattice
function col_srt_sukop!(sim::FreeSurfSim, F::Vector{Float64}, bounds::Matrix{Int64})
  lat = sim.lat;
  msm = sim.msm;
  const ni, nj = size(msm.rho);
  const nbounds = size(bounds, 2);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max, k = 1:lat.n 
      if sim.tracker.state[i,j] == GAS; continue; end
      rhoij = msm.rho[i,j];
      uij = msm.u[:,i,j];
      omegaij = msm.omega[i,j];
      feq = feq_incomp(lat, rhoij, uij, k);
      lat.f[k,i,j] = (omegaij * feq
                      + (1.0 - omegaij) * lat.f[k,i,j]
                      + lat.w[k] * lat.dt / lat.cssq * dot(F, lat.c[:,k]));
    end # body force is incorporated with w*dt/c_ssq * dot(f,c)
  end        
end

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param sim Simulation object
#! \param F Body force vector
#! \param bounds Boundaries that define active parts of the lattice
function col_srt_korner!(sim::FreeSurfSim, F::Vector{Float64}, bounds::Matrix{Int64})
  lat = sim.lat;
  msm = sim.msm;
  const ni, nj = size(msm.rho);
  const nbounds = size(bounds, 2);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max, k = 1:lat.n 
      if sim.tracker.state[i,j] == GAS; continue; end
      rhoij = msm.rho[i,j];
      uij = msm.u[:,i,j];
      omegaij = msm.omega[i,j];
      feq = feq_incomp(lat, rhoij, uij, k);
      lat.f[k,i,j] = (omegaij * feq
                      + (1.0 - omegaij) * lat.f[k,i,j]
                      + lat.w[k] * lat.dt / lat.cssq * dot(F, lat.c[:,k]));
    end # body force is incorporated with w*dt/c_ssq * dot(f,c)
  end        
end

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param sim Simulation object
#! \param F Body force vector
#! \param bounds Boundaries that define active parts of the lattice
function col_srt_guo!(sim::FreeSurfSim, F::Vector{Float64}, bounds::Matrix{Int64})
  lat = sim.lat;
  msm = sim.msm;
  const ni, nj = size(msm.rho);
  const nbounds = size(bounds, 2);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max, k = 1:lat.n
      if sim.tracker.state[i,j] == GAS; continue; end
      wk = lat.w[k];
      ck = lat.c[:,k];
      rhoij = msm.rho[i,j];
      uij = msm.u[:,i,j] + lat.dt / 2.0 * F;
      omegaij = msm.omega[i,j];

      fdlk = (1 - 0.5 * omegaij) * wk * dot(((ck - uij) / lat.cssq +
                  dot(ck, uij) / (lat.cssq * lat.cssq) * ck), F);
      feq = feq_incomp(lat, rhoij, uij, k);

      lat.f[k,i,j] = omegaij * feq + (1.0 - omegaij) * lat.f[k,i,j] + fdlk;
    end
  end
 end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param sim Simulation object
#! \param M Transmation matrix to map f from velocity space to momentum space
#! \param S (Sparse) diagonal relaxation matrix
#! \param bounds Boundaries that define active parts of the lattice
function col_mrt!(sim::FreeSurfSim, M::Matrix{Float64}, S::SparseMatrixCSC,
                  bounds::Array{Int64,2})
  lat = sim.lat;
  msm = sim.msm;
  const iM = inv(M);
  const ni, nj = size(msm.rho);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_n))
  feq = Array(Float64, lat.n);

  const nbounds = size(bounds, 2);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max 
      if sim.tracker.state[i,j] == GAS; continue; end
      rhoij = msm.rho[i,j];
      uij = msm.u[:,i,j];

      for k=1:lat.n; feq[k] = feq_incomp(lat, rhoij, uij, k); end

      f = lat.f[:,i,j];
      m = M * f;
      meq = M * feq;

      lat.f[:,i,j] = f - iM * S * (m - meq); # perform collision
    end
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param sim Simulation object
#! \param S (Sparse) diagonal relaxation matrix
#! \param bounds Boundaries that define active parts of the lattice
function col_mrt!(sim::FreeSurfSim, S::SparseMatrixCSC, bounds::Matrix{Int64})
  lat = sim.lat;
  msm = sim.msm;
  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const ni, nj = size(msm.rho);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_n))
  feq = Array(Float64, lat.n);

  const nbounds = size(bounds, 2);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max 
      if sim.tracker.state[i,j] == GAS; continue; end
      rhoij = msm.rho[i,j];
      uij = msm.u[:,i,j];

      for k=1:lat.n; feq[k] = feq_incomp(lat, rhoij, uij, k); end

      f = lat.f[:,i,j];
      m = M * f;
      meq = M * feq;

      lat.f[:,i,j] = f - iM * S * (m - meq); # perform collision
    end
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param sim Simulation object
#! \param S (Sparse) diagonal relaxation matrix
#! \param bounds Boundaries that define active parts of the lattice
#! \param F Body force vector
function col_mrt!(sim::FreeSurfSim, S::SparseMatrixCSC, F::Vector{Float64},
                  bounds::Matrix{Int64})
  lat = sim.lat;
  msm = sim.msm;
  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const ni, nj = size(msm.rho);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_n))
  feq = Array(Float64, lat.n);
  fdl = Array(Float64, lat.n);

  const nbounds = size(bounds, 2);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max 
      if sim.tracker.state[i,j] == GAS; continue; end
      rhoij = msm.rho[i,j];
      uij = msm.u[:,i,j] + lat.dt / 2.0 * F;
      omegaij = msm.omega[i,j];

      for k=1:lat.n
        wk = lat.w[k];
        ck = lat.c[:,k];
        feq[k] = feq_incomp(lat, rhoij, uij, k);
        fdl[k] = (1 - 0.5 * omegaij) * wk * dot(((ck - uij) / lat.cssq +
                  dot(ck, uij) / (lat.cssq * lat.cssq) * ck), F);
      end

      f = lat.f[:,i,j];
      m = M * f;
      meq = M * feq;

      lat.f[:,i,j] = f - iM * S * (m - meq) + fdl; # perform collision
    end
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param sim Simulation object
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \param mu_p Plastic viscosity
#! \param tau_y Yield stress
#! \param m Stress growth exponent
#! \param gamma_min Minimum strain rate to use in apparent viscosity calculation
#! \param bounds 2D array, each row is i_min, i_max, j_min, j_max
#! \param relax Relaxation factor for updating apparent viscosity
function col_mrt_bingham_explicit!(sim::FreeSurfSim, S::Function, mu_p::Number,
                                   tau_y::Number, m::Number,
                                   gamma_min::AbstractFloat,
                                   bounds::Matrix{Int64}, relax::Number = 1.0)
  lat = sim.lat;
  msm = sim.msm;
  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const ni, nj = size(msm.rho);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  feq = Array(Float64, lat.n);
  const nbounds = size(bounds, 2);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max
      if sim.tracker.state[i,j] == GAS; continue; end
      rhoij = msm.rho[i,j];
      uij = msm.u[:,i,j];

      for k=1:lat.n; feq[k] = feq_incomp(lat, rhoij, uij, k); end

      # initialize density, viscosity, and relaxation matrix at node i,j
      muij = @nu(msm.omega[i,j], lat.cssq, lat.dt);
      Sij = S(muij, rhoij, lat.cssq, lat.dt);

      f = lat.f[:,i,j];
      mij = M * f;
      meq = M * feq;
      fneq = f - feq;

      D = strain_rate_tensor(lat, rhoij, fneq, M, iM, Sij);
      gamma = @strain_rate(D);

      # update relaxation matrix
      gamma = gamma < gamma_min ? gamma_min : gamma;
      muij = (
               (1 - relax) * muij
                + relax * @mu_papanstasiou(mu_p, tau_y, m, gamma);
             );
      Sij = S(muij, rhoij, lat.cssq, lat.dt);

      lat.f[:,i,j] = f - iM * Sij * (mij - meq); # perform collision

      # update collision frequency matrix
      msm.omega[i,j] = @omega(muij, lat.cssq, lat.dt);
    end
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param sim Simulation object
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \param mu_p Plastic viscosity
#! \param tau_y Yield stress
#! \param m Stress growth exponent
#! \param gamma_min Minimum strain rate to use in apparent viscosity calculation
#! \param f Body force vector
#! \param bounds 2D array, each row is i_min, i_max, j_min, j_max
#! \param relax Relaxation factor for updating apparent viscosity
function col_mrt_bingham_explicit!(sim::FreeSurfSim, S::Function, mu_p::Number,
                                   tau_y::Number, m::Number,
                                   gamma_min::AbstractFloat,
                                   F::Vector{Float64}, bounds::Matrix{Int64},
                                   relax::Number = 1.0)
  lat = sim.lat;
  msm = sim.msm;
  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const ni, nj = size(msm.rho);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  feq = Array(Float64, lat.n);
  const nbounds = size(bounds, 2);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max
      if sim.tracker.state[i,j] == GAS; continue; end
      rhoij = msm.rho[i,j];
      uij = msm.u[:,i,j] + lat.dt / 2.0 * F;

      for k=1:lat.n; feq[k] = feq_incomp(lat, rhoij, uij, k); end

      # initialize density, viscosity, and relaxation matrix at node i,j
      muij = @nu(msm.omega[i,j], lat.cssq, lat.dt);
      Sij = S(muij, rhoij, lat.cssq, lat.dt);

      f = lat.f[:,i,j];
      mij = M * f;
      meq = M * feq;
      fneq = f - feq;

      D = strain_rate_tensor(lat, rhoij, fneq, M, iM, Sij);
      gamma = @strain_rate(D);

      # update relaxation matrix
      gamma = gamma < gamma_min ? gamma_min : gamma;
      muij = (
               (1 - relax) * muij
                + relax * @mu_papanstasiou(mu_p, tau_y, m, gamma);
             );
      Sij = S(muij, rhoij, lat.cssq, lat.dt);
      omegaij = @omega(muij, lat.cssq, lat.dt);

      fdl = Array(Float64, lat.n);
      for k=1:lat.n
        ck = lat.c[:,k];
        fdl[k] = (1 - 0.5 * omegaij) * lat.w[k] * dot(((ck - uij) / lat.cssq +
                  dot(ck, uij) / (lat.cssq * lat.cssq) * ck), F);
      end

      lat.f[:,i,j] = f - iM * Sij * (mij - meq) + fdl; # perform collision

      # update collision frequency matrix
      msm.omega[i,j] = omegaij;
    end
  end
end

#TODO: reconsider order of parameters... come up with a convenction
# phys models, const/material params, Functions, knobs

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param sim Simulation object
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \param mu_p Plastic viscosity
#! \param tau_y Yield stress
#! \param m Stress growth exponent
#! \param max_iters Maximum iterations for determining apparent viscosity
#! \param tol Tolerance for apparent viscosity convergence
#! \param gamma_min Minimum strain rate to use in apparent viscosity calculation
#! \param bounds 2D array, each row is i_min, i_max, j_min, j_max
#! \param relax Relaxation factor for updating apparent viscosity
function col_mrt_bingham_implicit!(sim::FreeSurfSim, S::Function, mu_p::Number,
                                   tau_y::Number, m::Number, max_iters::Int,
                                   tol::AbstractFloat,
                                   gamma_min::AbstractFloat,
                                   bounds::Matrix{Int64}, relax::Number = 1.0)
  lat = sim.lat;
  msm = sim.msm;
  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const ni, nj = size(msm.rho);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  feq = Array(Float64, lat.n);
  const nbounds = size(bounds, 2);

  #! Stream
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max
      if sim.tracker.state[i,j] == GAS; continue; end
      rhoij = msm.rho[i,j];
      uij = msm.u[:,i,j];

      for k=1:lat.n; feq[k] = feq_incomp(lat, rhoij, uij, k); end

      # initialize density, viscosity, and relaxation matrix at node i,j
      muij = @nu(msm.omega[i,j], lat.cssq, lat.dt);
      Sij = S(muij, rhoij, lat.cssq, lat.dt);

      f = lat.f[:,i,j];
      mij = M * f;
      meq = M * feq;
      fneq = f - feq;
      muo = muij;

      # iteratively determine mu
      iters = 0;
      mu_prev = muo;

      while true
        iters += 1;

        D = strain_rate_tensor(lat, rhoij, fneq, M, iM, Sij);
        gamma = @strain_rate(D);

        # update relaxation matrix
        gamma = gamma < gamma_min ? gamma_min : gamma;
        muij = (
                 (1 - relax) * muij
                  + relax * @mu_papanstasiou(mu_p, tau_y, m, gamma);
               );
        Sij = S(muij, rhoij, lat.cssq, lat.dt);

        # check for convergence
        if abs(mu_prev - muij) / muo <= tol
          break;
        end

        if iters > max_iters
          break;
        end

        mu_prev = muij;
      end

      lat.f[:,i,j] = f - iM * Sij * (mij - meq); # perform collision

      # update collision frequency matrix
      msm.omega[i,j] = @omega(muij, lat.cssq, lat.dt);
    end
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param sim Simulation object
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
function col_mrt_bingham_implicit!(sim::FreeSurfSim, S::Function, mu_p::Number,
                                   tau_y::Number, m::Number, max_iters::Int,
                                   tol::AbstractFloat,
                                   gamma_min::AbstractFloat,
                                   F::Vector{Float64},
                                   bounds::Matrix{Int64}, relax::Number = 1.0)
  lat = sim.lat;
  msm = sim.msm;
  const M = @DEFAULT_MRT_M();
  const iM = inv(M);
  const ni, nj = size(msm.rho);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  feq = Array(Float64, lat.n);
  const nbounds = size(bounds, 2);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max
      if sim.tracker.state[i,j] == GAS; continue; end
      rhoij = msm.rho[i,j];
      uij = msm.u[:,i,j] + lat.dt / 2.0 * F;

      for k=1:lat.n; feq[k] = feq_incomp(lat, rhoij, uij, k); end

      # initialize density, viscosity, and relaxation matrix at node i,j
      muij = @nu(msm.omega[i,j], lat.cssq, lat.dt);
      Sij = S(muij, rhoij, lat.cssq, lat.dt);

      f = lat.f[:,i,j];
      mij = M * f;
      meq = M * feq;
      fneq = f - feq;
      muo = muij;

      # iteratively determine mu
      iters = 0;
      mu_prev = muo;

      while true
        iters += 1;

        D = strain_rate_tensor(lat, rhoij, fneq, M, iM, Sij);
        gamma = @strain_rate(D);

        # update relaxation matrix
        gamma = gamma < gamma_min ? gamma_min : gamma;
        muij = (
                 (1 - relax) * muij
                  + relax * @mu_papanstasiou(mu_p, tau_y, m, gamma);
               );
        Sij = S(muij, rhoij, lat.cssq, lat.dt);

        # check for convergence
        if abs(mu_prev - muij) / muo <= tol
          break;
        end

        if iters > max_iters
          break;
        end

        mu_prev = muij;
      end

      const omegaij = @omega(muij, lat.cssq, lat.dt);

      fdl = Array(Float64, lat.n);
      for k=1:lat.n
        ck = lat.c[:,k];
        fdl[k] = (1 - 0.5 * omegaij) * lat.w[k] * dot(((ck - uij) / lat.cssq +
                  dot(ck, uij) / (lat.cssq * lat.cssq) * ck), F);
      end

      lat.f[:,i,j] = f - iM * Sij * (mij - meq) + fdl; # perform collision

      # update collision frequency matrix
      msm.omega[i,j] = omegaij;
    end
  end
end
