# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

const __constitutive_root__ = dirname(@__FILE__);
require(abspath(joinpath(__constitutive_root__, "..", "lattice.jl")));
require(abspath(joinpath(__constitutive_root__, "..", "multiscale.jl")));
require(abspath(joinpath(__constitutive_root__, "..", "sim", "simtypes.jl")));

#! Initialize a constant constitutive relationship
#!
#! \param mu Dynamic viscosity
#! \return Constitutive relation function
function init_constit_srt_const(mu::FloatingPoint)
  return (sim::AbstractSim, fneq::Vector{Float64}, i::Int, j::Int) -> begin
    return mu;
  end
end

#! Initialize a constant constitutive relationship
#!
#! \param mu Dynamic viscosity
#! \return Constitutive relation function
function init_constit_srt_const_local()
  return (sim::AbstractSim, fneq::Vector{Float64}, i::Int, j::Int) -> begin
    return @nu(sim.msm.omega[i,j], sim.lat.cssq, sim.lat.dt);
  end
end

#! Initialize an explicit bingham constitutive relationship
#!
#! \param mu_p Plastic vicosity
#! \param tau_y Yield stress
#! \param m Papanstasiou exponent
#! \param gamma_min Minimum allowable strain rate
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_bingham_explicit(mu_p::FloatingPoint,
                                           tau_y::FloatingPoint,
                                           m::Number,
                                           gamma_min::FloatingPoint,
                                           relax::Number = 1.0)

  return (sim::AbstractSim, fneq::Vector{Float64}, i::Int, j::Int) -> begin
    const rhoij = sim.msm.rho[i,j];
    const omegaij = sim.msm.omega[i,j];

    D = strain_rate_tensor(sim.lat, rhoij, fneq, omegaij);
    gamma = @strain_rate(D);

    # update relaxation matrix
    gamma = gamma < gamma_min ? gamma_min : gamma;
    muij = @nu(omegaij, sim.lat.cssq, sim.lat.dt);
    muij = (
             (1 - relax) * muij
              + relax * @mu_papanstasiou(mu_p, tau_y, m, gamma);
           );
    return muij;
  end
end

#! Initialize an implicit bingham constitutive relationship
#!
#! \param mu_p Plastic vicosity
#! \param tau_y Yield stress
#! \param m Papanstasiou exponent
#! \param gamma_min Minimum allowable strain rate
#! \param max_iters Maximum iterations
#! \param tol Convergence tolerance
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_bingham_implicit(mu_p::FloatingPoint,
                                           tau_y::FloatingPoint,
                                           m::Number,
                                           gamma_min::FloatingPoint,
                                           max_iters::Int,
                                           tol::FloatingPoint,
                                           relax::Number = 1.0)

  return (sim::AbstractSim, fneq::Vector{Float64}, i::Int, j::Int) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const rhoij = sim.msm.rho[i,j];
    omegaij = msm.omega[i,j];

    # initialize density, viscosity, and relaxation matrix at node i,j
    const muo = @nu(omegaij, sim.lat.cssq, sim.lat.dt);

    # iteratively determine mu
    iters = 0;
    mu_prev = muo;
    muij = muo;

    while true
      iters += 1;

      D = strain_rate_tensor(sim.lat, rhoij, fneq, omegaij);
      gamma = @strain_rate(D);

      # update relaxation matrix
      gamma = gamma < gamma_min ? gamma_min : gamma;
      muij = (
               (1 - relax) * muij
                + relax * @mu_papanstasiou(mu_p, tau_y, m, gamma);
             );
      omegaij = @omega(muij, sim.lat.cssq, sim.lat.dt);

      # check for convergence
      if abs(mu_prev - muij) / muo <= tol
        break;
      end

      if iters > max_iters
        break;
      end

      mu_prev = muij;
    end
    return muij;
  end
end

#! Initialize an explicit bingham constitutive relationship
#!
#! \param mu_p Plastic vicosity
#! \param tau_y Yield stress
#! \param rt_min Minimum relaxation time
#! \param rt_max Maximum relaxation time
#! \param m Papanstasiou exponent
#! \param gamma_min Minimum allowable strain rate
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_bingham_explicit(mu_p::FloatingPoint,
                                           tau_y::FloatingPoint,
                                           rt_min::FloatingPoint,
                                           rt_max::FloatingPoint,
                                           m::Number,
                                           gamma_min::FloatingPoint,
                                           relax::Number = 1.0)

  const inner_f = init_constit_srt_bingham_explicit(mu_p, tau_y, m, gamma_min, 
                                                    relax);

  return (sim::AbstractSim, fneq::Vector{Float64}, i::Int, j::Int) -> begin
    const mu_min = @nu(1.0/rt_min, sim.lat.cssq, sim.lat.dt);
    const mu_max = @nu(1.0/rt_min, sim.lat.cssq, sim.lat.dt);
    mu = inner_f(sim, fneq, i, j);
    if mu < mu_min
      return mu_min;
    elseif mu > mu_max
      return mu_max;
    else 
      return mu;
    end
  end
end

#! Initialize an implicit bingham constitutive relationship
#!
#! \param mu_p Plastic vicosity
#! \param tau_y Yield stress
#! \param rt_min Minimum relaxation time
#! \param rt_max Maximum relaxation time
#! \param m Papanstasiou exponent
#! \param gamma_min Minimum allowable strain rate
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_bingham_implicit(mu_p::FloatingPoint,
                                           tau_y::FloatingPoint,
                                           rt_min::FloatingPoint,
                                           rt_max::FloatingPoint,
                                           m::Number,
                                           gamma_min::FloatingPoint,
                                           max_iters::Int,
                                           tol::FloatingPoint,
                                           relax::Number = 1.0)

  const inner_f = init_constit_srt_bingham_implicit(mu_p, tau_y, m, gamma_min, 
                                                    max_iters, tol, relax);

  return (sim::AbstractSim, fneq::Vector{Float64}, i::Int, j::Int) -> begin
    const mu_min = @nu(1.0/rt_min, sim.lat.cssq, sim.lat.dt);
    const mu_max = @nu(1.0/rt_min, sim.lat.cssq, sim.lat.dt);
    mu = inner_f(sim, fneq, i, j);
    if mu < mu_min
      return mu_min;
    elseif mu > mu_max
      return mu_max;
    else 
      return mu;
    end
  end
end

#! Initialize a constant constitutive relationship
#!
#! \param mu Dynamic viscosity
#! \return Constitutive relation function
function init_constit_mrt_const(mu::FloatingPoint)
  return (sim::AbstractSim, fneq::Vector{Float64}, S::Function, 
          M::Matrix{Float64}, iM::Matrix{Float64}, i::Int, j::Int) -> begin
    return mu;
  end
end

#! Initialize a constant constitutive relationship
#!
#! \param mu Dynamic viscosity
#! \return Constitutive relation function
function init_constit_mrt_const_local()
  return (sim::AbstractSim, fneq::Vector{Float64}, S::Function, 
          M::Matrix{Float64}, iM::Matrix{Float64}, i::Int, j::Int) -> begin
    return @nu(sim.msm.omega[i,j], sim.lat.cssq, sim.lat.dt);
  end
end

#! Initialize an explicit bingham constitutive relationship
#!
#! \param mu_p Plastic vicosity
#! \param tau_y Yield stress
#! \param m Papanstasiou exponent
#! \param gamma_min Minimum allowable strain rate
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_mrt_bingham_explicit(mu_p::FloatingPoint,
                                           tau_y::FloatingPoint,
                                           m::Number,
                                           gamma_min::FloatingPoint,
                                           relax::Number = 1.0)

  return (sim::AbstractSim, fneq::Vector{Float64}, S::Function, 
          M::Matrix{Float64}, iM::Matrix{Float64}, i::Int, j::Int) -> begin
    const rhoij = sim.msm.rho[i,j];

    # initialize density, viscosity, and relaxation matrix at node i,j
    muij = @nu(sim.msm.omega[i,j], sim.lat.cssq, sim.lat.dt);
    Sij = S(muij, rhoij, sim.lat.cssq, sim.lat.dt);
    
    D = strain_rate_tensor(sim.lat, rhoij, fneq, M, iM, Sij);
    gamma = @strain_rate(D);

    # update relaxation matrix
    gamma = gamma < gamma_min ? gamma_min : gamma;
    muij = (
             (1 - relax) * muij
              + relax * @mu_papanstasiou(mu_p, tau_y, m, gamma);
           );
    return muij;
  end
end

#! Initialize an implicit bingham constitutive relationship
#!
#! \param mu_p Plastic vicosity
#! \param tau_y Yield stress
#! \param m Papanstasiou exponent
#! \param gamma_min Minimum allowable strain rate
#! \param max_iters Maximum iterations
#! \param tol Convergence tolerance
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_mrt_bingham_implicit(mu_p::FloatingPoint,
                                           tau_y::FloatingPoint,
                                           m::Number,
                                           gamma_min::FloatingPoint,
                                           max_iters::Int,
                                           tol::FloatingPoint = 1e-6,
                                           relax::Number = 1.0)

  return (sim::AbstractSim, fneq::Vector{Float64}, S::Function, 
          M::Matrix{Float64}, iM::Matrix{Float64}, i::Int, j::Int) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const rhoij = sim.msm.rho[i,j];

    # initialize density, viscosity, and relaxation matrix at node i,j
    const muo = @nu(msm.omega[i,j], lat.cssq, lat.dt);
    muij = muo;
    Sij = S(muij, rhoij, lat.cssq, lat.dt);

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
    return muij;
  end
end
