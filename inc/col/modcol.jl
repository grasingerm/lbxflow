# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

include("constitutive.jl");
include("forcing.jl");
include("equilibrium.jl");
include(joinpath("..", "lattice.jl"));
include("mrt_matrices.jl");
include(joinpath("..", "multiscale.jl"));
include(joinpath("..", "numerics.jl"));
include(joinpath("..", "sim", "simtypes.jl"));

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param bounds Boundaries that define active parts of the lattice
#! \param constit_relation_f Constitutive relationship
#! \return collision_function!(sim, bounds)
function init_col_srt(constit_relation_f::Function)
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        rhoij = msm.rho[i,j];
        uij = msm.u[:,i,j];
        feq = Array(Float64, lat.n); 
        fneq = Array(Float64, lat.n); 
        for k = 1:lat.n 
          feq[k] = feq_incomp(lat, rhoij, uij, k);
          fneq[k] = lat.f[k,i,j] - feq[k];
        end
        const mu = constit_relation_f(sim, fneq, i, j);
        const omega = @omega(mu, lat.cssq, lat.dt);
        for k = 1:lat.n
          lat.f[k,i,j] = (omega * feq[k] + (1.0 - omega) * lat.f[k,i,j]);
        end
        msm.omega[i,j] = omega;
      end
    end
  end
end

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param constit_relation_f Constitutive relationship
#! \param forcing_kf Forcing functions
#! \return collision_function!(sim, bounds)
function init_col_srt(constit_relation_f::Function,
                      forcing_kf::Tuple{Function, Function})
  const uf, colf = forcing_kf;
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        rhoij = msm.rho[i,j];
        uij = uf(lat, msm.u[:,i,j]);
        feq = Array(Float64, lat.n); 
        fneq = Array(Float64, lat.n); 
        for k = 1:lat.n 
          feq[k] = feq_incomp(lat, rhoij, uij, k);
          fneq[k] = lat.f[k,i,j] - feq[k];
        end
        const mu = constit_relation_f(sim, fneq, i, j);
        const omega = @omega(mu, lat.cssq, lat.dt);
        for k = 1:lat.n
          lat.f[k,i,j] = (omega * feq[k] + (1.0 - omega) * lat.f[k,i,j]
                          + colf(lat, omega, uij, k));
        end
        msm.omega[i,j] = omega;
      end
    end
  end
end

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param constit_relation_f Constitutive relationship
#! \param feq_f Equilibrium particle distribution function
#! \return collision_function!(sim, bounds)
function init_col_srt(constit_relation_f::Function, feq_f::Function)
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        rhoij = msm.rho[i,j];
        uij = msm.u[:,i,j];
        feq = Array(Float64, lat.n); 
        fneq = Array(Float64, lat.n); 
        for k = 1:lat.n 
          feq[k] = feq_f(lat, rhoij, uij, k);
          fneq[k] = lat.f[k,i,j] - feq[k];
        end
        const mu = constit_relation_f(sim, fneq, i, j);
        const omega = @omega(mu, lat.cssq, lat.dt);
        for k = 1:lat.n
          lat.f[k,i,j] = (omega * feq[k] + (1.0 - omega) * lat.f[k,i,j]);
        end
        msm.omega[i,j] = omega;
      end
    end
  end
end

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param constit_relation_f Constitutive relationship
#! \param forcing_kf Forcing functions
#! \param feq_f Equilibrium particle distribution function
#! \return collision_function!(sim, bounds)
function init_col_srt(constit_relation_f::Function,
                      forcing_kf::Tuple{Function, Function}, feq_f::Function)
  const uf, colf = forcing_kf;
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        rhoij = msm.rho[i,j];
        uij = uf(lat, msm.u[:,i,j]);
        feq = Array(Float64, lat.n); 
        fneq = Array(Float64, lat.n); 
        for k = 1:lat.n 
          feq[k] = feq_f(lat, rhoij, uij, k);
          fneq[k] = lat.f[k,i,j] - feq[k];
        end
        const mu = constit_relation_f(sim, fneq, i, j);
        const omega = @omega(mu, lat.cssq, lat.dt);
        for k = 1:lat.n
          lat.f[k,i,j] = (omega * feq[k] + (1.0 - omega) * lat.f[k,i,j]
                          + colf(lat, omega, uij, k));
        end
        msm.omega[i,j] = omega;
      end
    end
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param constit_relation_f Constitutive relationship
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \return collision_function!(sim, bounds)
function init_col_mrt(constit_relation_f::Function, S::Function)
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const M = @DEFAULT_MRT_M();
    const iM = @DEFAULT_MRT_IM();
    const ni, nj = size(msm.rho);

    # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
    feq = Array(Float64, lat.n);
    const nbounds = size(bounds, 2);

    #! Stream
    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        rhoij = msm.rho[i,j];
        uij = msm.u[:,i,j];

        for k=1:lat.n; feq[k] = feq_incomp(lat, rhoij, uij, k); end

        fij = lat.f[:,i,j];
        fneq = fij - feq;

        muij = constit_relation_f(sim, fneq, S, M, iM, i, j);
        Sij = S(muij, rhoij, lat.cssq, lat.dt);

        lat.f[:,i,j] = fij - iM * Sij * M * fneq; # perform collision

        # update collision frequency matrix
        msm.omega[i,j] = @omega(muij, lat.cssq, lat.dt);
      end
    end
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param constit_relation_f Constitutive relationship
#! \param forcing_kf Forcing functions
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \return collision_function!(sim, bounds)
function init_col_mrt(constit_relation_f::Function,
                      forcing_kf::Tuple{Function, Function}, S::Function)
  const uf, colf = forcing_kf;
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const M = @DEFAULT_MRT_M();
    const iM = @DEFAULT_MRT_IM();
    const ni, nj = size(msm.rho);

    # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
    feq = Array(Float64, lat.n);
    const nbounds = size(bounds, 2);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        rhoij = msm.rho[i,j];
        uij = uf(lat, msm.u[:,i,j]);

        for k=1:lat.n; feq[k] = feq_incomp(lat, rhoij, uij, k); end

        fij = lat.f[:,i,j];
        fneq = fij - feq;

        muij = constit_relation_f(sim, fneq, S, M, iM, i, j);
        omegaij = @omega(muij, lat.cssq, lat.dt);
        Sij = S(muij, rhoij, lat.cssq, lat.dt);

        fdl = Array(Float64, lat.n);
        for k = 1:lat.n
          fdl[k] = colf(lat, omegaij, uij, k);
        end
        lat.f[:,i,j] = fij - iM * Sij * M * fneq + fdl; # perform collision

        # update collision frequency matrix
        msm.omega[i,j] = omegaij;
      end
    end
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param constit_relation_f Constitutive relationship
#! \param feq_f Equilibrium particle distribution function
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \return collision_function!(sim, bounds)
function init_col_mrt(constit_relation_f::Function, feq_f::Function,
                      S::Function)
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const M = @DEFAULT_MRT_M();
    const iM = @DEFAULT_MRT_IM();
    const ni, nj = size(msm.rho);

    # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
    feq = Array(Float64, lat.n);
    const nbounds = size(bounds, 2);

    #! Stream
    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        rhoij = msm.rho[i,j];
        uij = msm.u[:,i,j];

        for k=1:lat.n; feq[k] = feq_f(lat, rhoij, uij, k); end

        fij = lat.f[:,i,j];
        fneq = fij - feq;

        muij = constit_relation_f(sim, fneq, S, M, iM, i, j);
        Sij = S(muij, rhoij, lat.cssq, lat.dt);

        lat.f[:,i,j] = fij - iM * Sij * M * fneq; # perform collision

        # update collision frequency matrix
        msm.omega[i,j] = @omega(muij, lat.cssq, lat.dt);
      end
    end
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param constit_relation_f Constitutive relationship
#! \param forcing_kf Forcing functions
#! \param feq_f Equilibrium particle distribution function
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \return collision_function!(sim, bounds)
function init_col_mrt(constit_relation_f::Function,
                      forcing_kf::Tuple{Function, Function}, feq_f::Function,
                      S::Function)
  const uf, colf = forcing_kf;
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const M = @DEFAULT_MRT_M();
    const iM = @DEFAULT_MRT_IM();
    const ni, nj = size(msm.rho);

    # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
    feq = Array(Float64, lat.n);
    const nbounds = size(bounds, 2);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      for j = j_min:j_max, i = i_min:i_max
        rhoij = msm.rho[i,j];
        uij = uf(lat, msm.u[:,i,j]);

        for k=1:lat.n; feq[k] = feq_f(lat, rhoij, uij, k); end

        fij = lat.f[:,i,j];
        fneq = fij - feq;

        muij = constit_relation_f(sim, fneq, S, M, iM, i, j);
        omegaij = @omega(muij, lat.cssq, lat.dt);
        Sij = S(muij, rhoij, lat.cssq, lat.dt);

        fdl = Array(Float64, lat.n);
        for k = 1:lat.n
          fdl[k] = colf(lat, omegaij, uij, k);
        end
        lat.f[:,i,j] = fij - iM * Sij * M * fneq + fdl; # perform collision

        # update collision frequency matrix
        msm.omega[i,j] = omegaij;
      end
    end
  end
end
