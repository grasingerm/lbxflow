# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Two-phase constitutive equation that relaxs interface
#=type M2PhaseConstit
  constit_r::LBXFunction;
  constit_b::LBXFunction;

  M2PhaseConstit(constit_r::LBXFunction, constit_b::LBXFunction, δ::Real) = (
    new(constit_r, constit_b, δ));
end

function call(mc::M2PhaseConstit, sim::M2PhaseSim, args...)
  const ρ_r   =   sim.simr.msm.rho[i, j] 
  const ρ_b   =   sim.simb.msm.rho[i, j] 
  const ψ     =   (ρ_r - ρ_b) / (ρ_r + ρ_b);
end=#

# call definition for constant constitutive relationship
function call(cc::_ConstConstit, sim::AbstractSim, 
              fneq::AbstractArray{Float64, 1}, i::Int, j::Int)
  return cc.μ;
end

function call(cc::_ConstConstit, sim::AbstractSim, 
              fneq::AbstractArray{Float64, 1}, S::Function, 
              M::AbstractArray{Float64, 2}, iM::AbstractArray{Float64, 2}, 
              i::Int, j::Int)
  return cc.μ;
end


#! Initialize a constant constitutive relationship
#!
#! \param mu Dynamic viscosity
#! \return Constitutive relation function
function init_constit_srt_const(mu::AbstractFloat)
  return _ConstConstit(mu);
end

#! Initialize an explicit bingham constitutive relationship
#!
#! \param mu_p Plastic vicosity
#! \param tau_y Yield stress
#! \param m Papanstasiou exponent
#! \param gamma_min Minimum allowable strain rate
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_bingham_explicit(mu_p::AbstractFloat,
                                           tau_y::AbstractFloat,
                                           m::Number,
                                           gamma_min::AbstractFloat,
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
function init_constit_srt_bingham_implicit(mu_p::AbstractFloat,
                                           tau_y::AbstractFloat,
                                           m::Number,
                                           gamma_min::AbstractFloat,
                                           max_iters::Int,
                                           tol::AbstractFloat,
                                           relax::Number = 1.0)

  return (sim::AbstractSim, fneq::Vector{Float64}, i::Int, j::Int) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const rhoij = sim.msm.rho[i,j];
    omegaij = msm.omega[i,j];

    # initialize density, viscosity, and relaxation matrix at node i,j
    const muo = (omegaij != 0.0) ? @nu(omegaij, sim.lat.cssq, sim.lat.dt) : 1.0;
    @assert(!isnan(muo), "Initial guess is for mu is NaN");

    # iteratively determine mu
    iters = 0;
    mu_prev = muo;
    muij = muo;

    while true
      iters += 1;

      D = strain_rate_tensor(sim.lat, rhoij, fneq, omegaij);
      gamma = @strain_rate(D);
      @assert(!isnan(gamma), "GAMMA IS NAN");

      # update relaxation matrix
      gamma = gamma < gamma_min ? gamma_min : gamma;
      muij = (
               (1 - relax) * muij
                + relax * @mu_papanstasiou(mu_p, tau_y, m, gamma);
             );
      omegaij = @omega(muij, sim.lat.cssq, sim.lat.dt);
      @assert(!isnan(omegaij), "This shit should never be NaN");

      # check for convergence
      if abs(mu_prev - muij) / muo <= tol
        break;
      end

      if iters > max_iters
        #warn("Solution for constitutive equation did not converge");
        #@show i, j;
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
function init_constit_srt_bingham_explicit(mu_p::AbstractFloat,
                                           tau_y::AbstractFloat,
                                           rt_min::AbstractFloat,
                                           rt_max::AbstractFloat,
                                           m::Number,
                                           gamma_min::AbstractFloat,
                                           relax::Number = 1.0)

  const inner_f = init_constit_srt_bingham_explicit(mu_p, tau_y, m, gamma_min, 
                                                    relax);

  return (sim::AbstractSim, fneq::Vector{Float64}, i::Int, j::Int) -> begin
    const mu_min = sim.lat.cssq*sim.lat.dt * (rt_min - 0.5);
    const mu_max = sim.lat.cssq*sim.lat.dt * (rt_max - 0.5);
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
#! \param max_iters Maximum iterations
#! \param tol Convergence tolerance
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_bingham_implicit(mu_p::AbstractFloat,
                                           tau_y::AbstractFloat,
                                           rt_min::AbstractFloat,
                                           rt_max::AbstractFloat,
                                           m::Number,
                                           gamma_min::AbstractFloat,
                                           max_iters::Int,
                                           tol::AbstractFloat,
                                           relax::Number = 1.0)

  const inner_f = init_constit_srt_bingham_implicit(mu_p, tau_y, m, gamma_min, 
                                                    max_iters, tol, relax);

  return (sim::AbstractSim, fneq::Vector{Float64}, i::Int, j::Int) -> begin
    const mu_min = sim.lat.cssq*sim.lat.dt * (rt_min - 0.5);
    const mu_max = sim.lat.cssq*sim.lat.dt * (rt_max - 0.5);
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

#! Initialize an explicit herschel-bulkley constitutive relationship
#!
#! \param k Flow consistency index
#! \param n Power law index
#! \param tau_y Yield stress
#! \param m Papanstasiou exponent
#! \param gamma_min Minimum allowable strain rate
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_hb_explicit(k::AbstractFloat, n::Number,
                                      tau_y::AbstractFloat,
                                      m::Number,
                                      gamma_min::AbstractFloat,
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
              + relax * @mu_papanstasiou(k * gamma^(n-1), tau_y, m, gamma);
           );
    return muij;
  end
end

#! Initialize an implicit herschel-bulkley constitutive relationship
#!
#! \param k Flow consistency index
#! \param n Power law index
#! \param tau_y Yield stress
#! \param m Papanstasiou exponent
#! \param gamma_min Minimum allowable strain rate
#! \param max_iters Maximum iterations
#! \param tol Convergence tolerance
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_hb_implicit(k::AbstractFloat, n::Number,
                                      tau_y::AbstractFloat,
                                      m::Number,
                                      gamma_min::AbstractFloat,
                                      max_iters::Int,
                                      tol::AbstractFloat,
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
                + relax * @mu_papanstasiou(k * gamma^(n-1), tau_y, m, gamma);
             );
      omegaij = @omega(muij, sim.lat.cssq, sim.lat.dt);

      # check for convergence
      if abs(mu_prev - muij) / muo <= tol
        break;
      end

      if iters > max_iters
        #warn("Solution for constitutive equation did not converge");
        #@show i, j;
        break;
      end

      mu_prev = muij;
    end
    return muij;
  end
end

#! Initialize an explicit power law constitutive relationship
#!
#! \param k Flow consistency index
#! \param n Power law index
#! \param gamma_min Minimum allowable strain rate
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_power_law_explicit(k::AbstractFloat, n::Number,
                                             gamma_min::AbstractFloat,
                                             relax::Number = 1.0)

  return (sim::AbstractSim, fneq::Vector{Float64}, i::Int, j::Int) -> begin
    const rhoij = sim.msm.rho[i,j];
    const omegaij = sim.msm.omega[i,j];
    muij = @nu(omegaij, sim.lat.cssq, sim.lat.dt);

    D = strain_rate_tensor(sim.lat, rhoij, fneq, omegaij);
    gamma = @strain_rate(D);
    gamma = gamma < gamma_min ? gamma_min : gamma;

    # update relaxation matrix
    muij = (
             (1 - relax) * muij
              + relax * k * gamma^(n-1);
           );
    return muij;
  end
end

#! Initialize an implicit power law constitutive relationship
#!
#! \param k Flow consistency index
#! \param n Power law index
#! \param gamma_min Minimum allowable strain rate
#! \param max_iters Maximum iterations
#! \param tol Convergence tolerance
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_power_law_implicit(k::AbstractFloat, n::Number,
                                             gamma_min::AbstractFloat,
                                             max_iters::Int, tol::AbstractFloat,
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

      D = strain_rate_tensor(lat, rhoij, fneq, omegaij);
      gamma = @strain_rate(D);
      gamma = gamma < gamma_min ? gamma_min : gamma;

      # update relaxation matrix
      muij = (
               (1 - relax) * muij
                + relax * k * gamma^(n-1);
             );
      omegaij = @omega(muij, lat.cssq, lat.dt);

      # check for convergence
      if abs(mu_prev - muij) / muo <= tol
        break;
      end

      if iters > max_iters
        #warn("Solution for constitutive equation did not converge");
        #@show i, j;
        break;
      end

      mu_prev = muij;
    end
    return muij;
  end
end

#! Initialize an implicit power law constitutive relationship
#!
#! \param k Flow consistency index
#! \param n Power law index
#! \param rt_min Minimum relaxation time
#! \param rt_max Maximum relaxation time
#! \param gamma_min Minimum allowable strain rate
#! \param max_iters Maximum iterations
#! \param tol Convergence tolerance
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_power_law_implicit(k::AbstractFloat,
                                             n::AbstractFloat,
                                             rt_min::AbstractFloat,
                                             rt_max::AbstractFloat,
                                             gamma_min::AbstractFloat,
                                             max_iters::Int,
                                             tol::AbstractFloat,
                                             relax::Number = 1.0)

  const inner_f = init_constit_srt_power_law_implicit(k, n, gamma_min, 
                                                      max_iters, tol, relax);

  return (sim::AbstractSim, fneq::Vector{Float64}, i::Int, j::Int) -> begin
    const mu_min = @nu(1/rt_min, sim.lat.cssq, sim.lat.dt);
    const mu_max = @nu(1/rt_max, sim.lat.cssq, sim.lat.dt);
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

#! Initialize an explicit casson constitutive relationship
#!
#! \param k Flow consistency index
#! \param n Power law index
#! \param gamma_min Minimum allowable strain rate
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_casson_explicit(k::AbstractFloat, n::Number,
                                          gamma_min::AbstractFloat,
                                          relax::Number = 1.0)

  return (sim::AbstractSim, fneq::Vector{Float64}, i::Int, j::Int) -> begin
    const rhoij = sim.msm.rho[i,j];
    const omegaij = sim.msm.omega[i,j];
    muij = @nu(omegaij, sim.lat.cssq, sim.lat.dt);

    D = strain_rate_tensor(sim.lat, rhoij, fneq, omegaij);
    gamma = @strain_rate(D);
    gamma = gamma < gamma_min ? gamma_min : gamma;

    # update relaxation matrix
    muij = (
             (1 - relax) * muij
              + relax * k * gamma^(n-1);
           );
    return muij;
  end
end

#! Initialize an explicit power law constitutive relationship
#!
#! \param k Flow consistency index
#! \param n Power law index
#! \param gamma_min Minimum allowable strain rate
#! \param max_iters Maximum iterations
#! \param tol Convergence tolerance
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_srt_power_law_implicit(k::AbstractFloat, n::Number,
                                             gamma_min::AbstractFloat,
                                             max_iters::Int, tol::AbstractFloat,
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

      D = strain_rate_tensor(lat, rhoij, fneq, omegaij);
      gamma = @strain_rate(D);
      gamma = gamma < gamma_min ? gamma_min : gamma;

      # update relaxation matrix
      muij = (
               (1 - relax) * muij
                + relax * k * gamma^(n-1);
             );
      omegaij = @omega(muij, lat.cssq, lat.dt);

      # check for convergence
      if abs(mu_prev - muij) / muo <= tol
        break;
      end

      if iters > max_iters
        #warn("Solution for constitutive equation did not converge");
        #@show i, j;
        break;
      end

      mu_prev = muij;
    end
    return muij;
  end
end

#! Initialize a constant constitutive relationship
#!
#! \param mu Dynamic viscosity
#! \return Constitutive relation function
function init_constit_mrt_const(mu::AbstractFloat)
  return _ConstConstit(mu);
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
function init_constit_mrt_bingham_explicit(mu_p::AbstractFloat,
                                           tau_y::AbstractFloat,
                                           m::Number,
                                           gamma_min::AbstractFloat,
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
function init_constit_mrt_bingham_implicit(mu_p::AbstractFloat,
                                           tau_y::AbstractFloat,
                                           m::Number,
                                           gamma_min::AbstractFloat,
                                           max_iters::Int,
                                           tol::AbstractFloat = 1e-6,
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
        #warn("Solution for constitutive equation did not converge");
        #@show i, j;
        break;
      end

      mu_prev = muij;
    end
    return muij;
  end
end

#! Initialize an explicit herschel-bulkley constitutive relationship
#!
#! \param k Flow consistency index
#! \param n Power law index
#! \param tau_y Yield stress
#! \param m Papanstasiou exponent
#! \param gamma_min Minimum allowable strain rate
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_mrt_hb_explicit(k::AbstractFloat, n::Number,
                                      tau_y::AbstractFloat,
                                      m::Number,
                                      gamma_min::AbstractFloat,
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
              + relax * @mu_papanstasiou(k * gamma^(n-1), tau_y, m, gamma);
           );
    return muij;
  end
end

#! Initialize an implicit herschel-bulkley constitutive relationship
#!
#! \param k Flow consistency index
#! \param n Power law index
#! \param tau_y Yield stress
#! \param m Papanstasiou exponent
#! \param gamma_min Minimum allowable strain rate
#! \param max_iters Maximum iterations
#! \param tol Convergence tolerance
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_mrt_hb_implicit(k::AbstractFloat, n::Number,
                                      tau_y::AbstractFloat,
                                      m::Number,
                                      gamma_min::AbstractFloat,
                                      max_iters::Int,
                                      tol::AbstractFloat = 1e-6,
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
                + relax * @mu_papanstasiou(k * gamma^(n-1), tau_y, m, gamma);
             );
      Sij = S(muij, rhoij, lat.cssq, lat.dt);

      # check for convergence
      if abs(mu_prev - muij) / muo <= tol
        break;
      end

      if iters > max_iters
        #warn("Solution for constitutive equation did not converge");
        #@show i, j;
        break;
      end

      mu_prev = muij;
    end
    return muij;
  end
end

#! Initialize an explicit power law constitutive relationship
#!
#! \param k Flow consistency index
#! \param n Power law index
#! \param gamma_min Minimum allowable strain rate
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_mrt_power_law_explicit(k::AbstractFloat, n::Number,
                                             gamma_min::AbstractFloat,
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
              + relax * k * gamma^(n-1);
           );
    return muij;
  end
end

#! Initialize an implicit power law constitutive relationship
#!
#! \param k Flow consistency index
#! \param n Power law index
#! \param gamma_min Minimum allowable strain rate
#! \param max_iters Maximum iterations
#! \param tol Convergence tolerance
#! \param relax Relaxation coefficient
#! \return Constitutive relation function
function init_constit_mrt_power_law_implicit(k::AbstractFloat, n::Number,
                                             gamma_min::AbstractFloat,
                                             max_iters::Int, 
                                             tol::AbstractFloat = 1e-6,
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
                + relax * k * gamma^(n-1);
             );
      Sij = S(muij, rhoij, lat.cssq, lat.dt);

      # check for convergence
      if abs(mu_prev - muij) / muo <= tol
        break;
      end

      if iters > max_iters
        #warn("Solution for constitutive equation did not converge");
        #@show i, j;
        break;
      end

      mu_prev = muij;
    end
    return muij;
  end
end

