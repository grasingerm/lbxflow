# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param sim Simulation object
#! \param bounds Boundaries that define active parts of the lattice
#! \param constit_relation_f Constitutive relationship
function init_pcol_srt(constit_relation_f::Function)
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      lat.f[:,i_min:i_max,j_min:j_max], msm.omega[i_min:i_max, j_min, j_max] = (
        pmap((i,j) -> begin
          rhoij = msm.rho[i,j];
          uij = msm.u[:,i,j];
          f = Array(Float64, lat.n);
          feq = Array(Float64, lat.n); 
          fneq = Array(Float64, lat.n); 
          for k = 1:lat.n 
            feq[k] = feq_incomp(lat, rhoij, uij, k);
            fneq[k] = lat.f[k,i,j] - feq[k];
          end
          const mu = constit_relation_f(sim, fneq, i, j);
          const omega = @omega(mu, lat.cssq, lat.dt);
          for k = 1:lat.n
            f[k] = (omega * feq[k] + (1.0 - omega) * lat.f[k,i,j]);
          end
          return f, omega;
        end, [(i,j) for i=i_min:i_max, j=j_min:j_max])
      );
    end
  end
end

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param sim Simulation object
#! \param bounds Boundaries that define active parts of the lattice
#! \param constit_relation_f Constitutive relationship
#! \param forcing_kf Forcing functions
function init_pcol_srt(constit_relation_f::Function,
                       forcing_kf::Tuple{Function, Function})
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
    lat = sim.lat;
    msm = sim.msm;
    const ni, nj = size(msm.rho);
    const nbounds = size(bounds, 2);

    for r = 1:nbounds
      i_min, i_max, j_min, j_max = bounds[:,r];
      lat.f[:,i_min:i_max,j_min:j_max], msm.omega[i_min:i_max, j_min, j_max] = (
        pmap((i,j) -> begin
          rhoij = msm.rho[i,j];
          uij = uf(lat, msm.u[:,i,j]);
          f = Array(Float64, lat.n);
          feq = Array(Float64, lat.n); 
          fneq = Array(Float64, lat.n); 
          for k = 1:lat.n 
            feq[k] = feq_incomp(lat, rhoij, uij, k);
            fneq[k] = lat.f[k,i,j] - feq[k];
          end
          const mu = constit_relation_f(sim, fneq, i, j);
          const omega = @omega(mu, lat.cssq, lat.dt);
          for k = 1:lat.n
            f[k] = ((omega * feq[k] + (1.0 - omega) * lat.f[k,i,j])
                    + colf(lat, omega, uij, k));
          end
          return f, omega;
        end, [(i,j) for i=i_min:i_max, j=j_min:j_max])
      );
    end
  end
end

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param sim Simulation object
#! \param bounds Boundaries that define active parts of the lattice
#! \param constit_relation_f Constitutive relationship
#! \param feq_f Equilibrium particle distribution function
function init_pcol_srt(constit_relation_f::Function, feq_f::Function)
  error("Not yet implemented.");
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
#! \param sim Simulation object
#! \param bounds Boundaries that define active parts of the lattice
#! \param constit_relation_f Constitutive relationship
#! \param forcing_kf Forcing functions
#! \param feq_f Equilibrium particle distribution function
function init_pcol_srt(constit_relation_f::Function,
                       forcing_kf::Tuple{Function, Function}, feq_f::Function)
  error("Not yet implemented.");
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
#! \param sim Simulation object
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \param bounds 2D array, each row is i_min, i_max, j_min, j_max
#! \param constit_relation_f Constitutive relationship
function init_pcol_mrt(constit_relation_f::Function)
  error("Not yet implemented.");
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
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
        rhoij = msm.rho[i,j];
        uij = msm.u[:,i,j];

        for k=1:lat.n; feq[k] = feq_incomp(lat, rhoij, uij, k); end

        f = lat.f[:,i,j];
        mij = M * f;
        meq = M * feq;
        fneq = f - feq;

        muij = constit_relation_f(sim, S, M, iM, f, feq, fneq, mij, meq, i, j);
        Sij = S(mu, rhoij, lat.cssq, lat.dt);

        lat.f[:,i,j] = f - iM * Sij * (mij - meq); # perform collision

        # update collision frequency matrix
        msm.omega[i,j] = @omega(muij, lat.cssq, lat.dt);
      end
    end
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param sim Simulation object
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \param bounds 2D array, each row is i_min, i_max, j_min, j_max
#! \param constit_relation_f Constitutive relationship
#! \param forcing_kf Forcing functions
function init_pcol_mrt(constit_relation_f::Function,
                       forcing_kf::Tuple{Function, Function})
  error("Not yet implemented.");
  const uf, colf = forcing_kf;
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
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
        rhoij = msm.rho[i,j];
        uij = uf(msm.u[:,i,j]);

        for k=1:lat.n; feq[k] = feq_incomp(lat, rhoij, uij, k); end

        f = lat.f[:,i,j];
        mij = M * f;
        meq = M * feq;
        fneq = f - feq;

        muij = constit_relation_f(sim, S, M, iM, f, feq, fneq, mij, meq, i, j);
        Sij = S(mu, rhoij, lat.cssq, lat.dt);

        fdl = Array(Float64, lat.n);
        for k = 1:lat.n
          fdl[k] = colf(lat, omega, uij, k);
        end
        lat.f[:,i,j] = f - iM * Sij * (mij - meq) + fdl; # perform collision

        # update collision frequency matrix
        msm.omega[i,j] = @omega(muij, lat.cssq, lat.dt);
      end
    end
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param sim Simulation object
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \param bounds 2D array, each row is i_min, i_max, j_min, j_max
#! \param constit_relation_f Constitutive relationship
#! \param feq_f Equilibrium particle distribution function
function init_pcol_mrt(constit_relation_f::Function, feq_f::Function)
  error("Not yet implemented.");
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
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
        rhoij = msm.rho[i,j];
        uij = msm.u[:,i,j];

        for k=1:lat.n; feq[k] = feq_f(lat, rhoij, uij, k); end

        f = lat.f[:,i,j];
        mij = M * f;
        meq = M * feq;
        fneq = f - feq;

        muij = constit_relation_f(sim, S, M, iM, f, feq, fneq, mij, meq, i, j);
        Sij = S(mu, rhoij, lat.cssq, lat.dt);

        lat.f[:,i,j] = f - iM * Sij * (mij - meq); # perform collision

        # update collision frequency matrix
        msm.omega[i,j] = @omega(muij, lat.cssq, lat.dt);
      end
    end
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param sim Simulation object
#! \param S Function that returns (sparse) diagonal relaxation matrix
#! \param bounds 2D array, each row is i_min, i_max, j_min, j_max
#! \param constit_relation_f Constitutive relationship
#! \param forcing_kf Forcing functions
#! \param feq_f Equilibrium particle distribution function
function init_pcol_mrt(constit_relation_f::Function,
                       forcing_kf::Tuple{Function, Function}, feq_f::Function)
  error("Not yet implemented.");
  const uf, colf = forcing_kf;
  return (sim::Sim, bounds::Matrix{Int64}) -> begin
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
        rhoij = msm.rho[i,j];
        uij = uf(msm.u[:,i,j]);

        for k=1:lat.n; feq[k] = feq_f(lat, rhoij, uij, k); end

        f = lat.f[:,i,j];
        mij = M * f;
        meq = M * feq;
        fneq = f - feq;

        muij = constit_relation_f(sim, S, M, iM, f, feq, fneq, mij, meq, i, j);
        Sij = S(mu, rhoij, lat.cssq, lat.dt);

        fdl = Array(Float64, lat.n);
        for k = 1:lat.n
          fdl[k] = colf(lat, omega, uij, k);
        end
        lat.f[:,i,j] = f - iM * Sij * (mij - meq) + fdl; # perform collision

        # update collision frequency matrix
        msm.omega[i,j] = @omega(muij, lat.cssq, lat.dt);
      end
    end
  end
end
