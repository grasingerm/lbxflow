# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Multiple relaxation time collision operator 
type MRT <: ColFunction
  feq_f::LBXFunction;
  constit_relation_f::LBXFunction;
  M::Matrix{Float64};
  iM::Matrix{Float64};
  S::LBXFunction;

  function MRT(constit_relation_f::LBXFunction; feq_f::LBXFunction=feq_incomp,
               S::LBXFunction=S_fallah)
    const M   = @DEFAULT_MRT_M;
    const iM  = @DEFAULT_MRT_IM;
    return new(feq_f, constit_relation_f, M, iM, S);
  end
end

#! Multiple relaxation time collision operator with forcing 
type MRT_F <: ColFunction
  feq_f::LBXFunction;
  constit_relation_f::LBXFunction;
  forcing_f::Force;
  M::Matrix{Float64};
  iM::Matrix{Float64};
  S::LBXFunction;

  function MRT_F(constit_relation_f::LBXFunction,
                 forcing_f::Force; 
                 feq_f::LBXFunction=feq_incomp, S::LBXFunction=S_fallah)
    const M   = @DEFAULT_MRT_M;
    const iM  = @DEFAULT_MRT_IM;
    return new(feq_f, constit_relation_f, forcing_f, M, iM, S);
  end
end

#! Perform collision over box regions of the domain
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function call(col_f::MRT, sim::AbstractSim, bounds::Matrix{Int64})
  lat           =   sim.lat;
  msm           =   sim.msm;
  const ni, nj  =   size(msm.rho);
  const nbounds =   size(bounds, 2);

  for r = 1:nbounds
    @inbounds i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max

      @inbounds rhoij =   msm.rho[i,j];
      @inbounds uij   =   msm.u[:,i,j];
      feq             =   Vector{Float64}(lat.n);
      fneq            =   Vector{Float64}(lat.n);

      @simd for k = 1:lat.n 
        @inbounds feq[k]          =   col_f.feq_f(lat, rhoij, uij, k);
        @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
      end

      const mu        =   col_f.constit_relation_f(sim, fneq, col_f.S, col_f.M, 
                                                   col_f.iM, i, j);
      const Sij       =   col_f.S(mu, rhoij, lat.cssq, lat.dt);

      @inbounds lat.f[:, i, j] -= col_f.iM * Sij * col_f.M * fneq;

      # update collision frequency matrix
      @inbounds msm.omega[i,j] = @omega(mu, lat.cssq, lat.dt);

    end
  end
end

#! Perform collision over box regions of the domain
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function call(col_f::MRT_F, sim::AbstractSim, bounds::Matrix{Int64})
  lat           =   sim.lat;
  msm           =   sim.msm;
  const ni, nj  =   size(msm.rho);
  const nbounds =   size(bounds, 2);

  for r = 1:nbounds
    @inbounds i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max

      @inbounds rhoij =   msm.rho[i,j];
      @inbounds uij   =   col_f.forcing_f[1](sim, i, j);
      feq             =   Vector{Float64}(lat.n);
      fneq            =   Vector{Float64}(lat.n);

      @simd for k = 1:lat.n 
        @inbounds feq[k]          =   col_f.feq_f(lat, rhoij, uij, k);
        @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
      end

      const mu    =   col_f.constit_relation_f(sim, fneq, col_f.S, col_f.M, 
                                               col_f.iM, i, j);
      const omega =   @omega(mu, lat.cssq, lat.dt);
      const Sij   =   col_f.S(mu, rhoij, lat.cssq, lat.dt);
      const F     =   map(k -> col_f.forcing_f[2](sim, omega, k, i, j), 1:lat.n);

      @inbounds lat.f[:,i,j]  -= col_f.iM * Sij * col_f.M * fneq - F;

      # update collision frequency matrix
      @inbounds msm.omega[i,j] = omega;

    end
  end
end

#! Perform collision over box regions of the domain
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function call(col_f::MRT, sim::FreeSurfSim, bounds::Matrix{Int64})
  lat           =   sim.lat;
  msm           =   sim.msm;
  const ni, nj  =   size(msm.rho);
  const nbounds =   size(bounds, 2);

  for r = 1:nbounds
    @inbounds i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max

      @inbounds if sim.tracker.state[i, j] != GAS

        @inbounds rhoij =   msm.rho[i,j];
        @inbounds uij   =   msm.u[:,i,j];
        feq             =   Vector{Float64}(lat.n);
        fneq            =   Vector{Float64}(lat.n);

        @simd for k = 1:lat.n 
          @inbounds feq[k]          =   col_f.feq_f(lat, rhoij, uij, k);
          @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
        end

        const mu    =   col_f.constit_relation_f(sim, fneq, col_f.S, col_f.M, 
                                                 col_f.iM, i, j);
        const Sij       =   col_f.S(mu, rhoij, lat.cssq, lat.dt);

        @inbounds lat.f[:,i,j]  -= col_f.iM * Sij * col_f.M * fneq;

        # update collision frequency matrix
        @inbounds msm.omega[i,j] = @omega(mu, lat.cssq, lat.dt);

      end

    end
  end
end

#! Perform collision over box regions of the domain
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function call(col_f::MRT_F, sim::FreeSurfSim, bounds::Matrix{Int64})
  lat           =   sim.lat;
  msm           =   sim.msm;
  const ni, nj  =   size(msm.rho);
  const nbounds =   size(bounds, 2);

  for r = 1:nbounds
    @inbounds i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max

      @inbounds if sim.tracker.state[i, j] != GAS

        @inbounds rhoij =   msm.rho[i,j];
        @inbounds uij   =   col_f.forcing_f[1](sim, i, j);
        feq             =   Vector{Float64}(lat.n);
        fneq            =   Vector{Float64}(lat.n);

        @simd for k = 1:lat.n 
          @inbounds feq[k]          =   col_f.feq_f(lat, rhoij, uij, k);
          @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
        end

        const mu    =   col_f.constit_relation_f(sim, fneq, col_f.S, col_f.M, 
                                                 col_f.iM, i, j);
        const omega =   @omega(mu, lat.cssq, lat.dt);
        const Sij   =   col_f.S(mu, rhoij, lat.cssq, lat.dt);
        const F     =   map(k -> col_f.forcing_f[2](sim, omega, k, i, j), 1:lat.n);

        @inbounds lat.f[:, i, j]  -= col_f.iM * Sij * col_f.M * fneq - F;

        # update collision frequency matrix
        @inbounds msm.omega[i,j] = omega;

      end

    end
  end
end

#! Perform collision over active regions of the domain
#!
#! \param   sim           Simulation
#! \param   active_cells  Active flags for domain
function call(col_f::MRT, sim::AbstractSim, active_cells::Matrix{Bool})
  lat           =   sim.lat;
  msm           =   sim.msm;
  const ni, nj  =   size(msm.rho);

  for j = 1:nj, i = 1:ni

    @inbounds if active_cells[i, j]

      @inbounds rhoij =   msm.rho[i,j];
      @inbounds uij   =   msm.u[:,i,j];
      feq             =   Vector{Float64}(lat.n);
      fneq            =   Vector{Float64}(lat.n);

      @simd for k = 1:lat.n 
        @inbounds feq[k]          =   col_f.feq_f(lat, rhoij, uij, k);
        @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
      end

      const mu        =   col_f.constit_relation_f(sim, fneq, col_f.S, col_f.M, 
                                                   col_f.iM, i, j);
      const Sij       =   col_f.S(mu, rhoij, lat.cssq, lat.dt);

      @inbounds lat.f[:, i, j] -= col_f.iM * Sij * col_f.M * fneq;

      # update collision frequency matrix
      @inbounds msm.omega[i,j] = @omega(mu, lat.cssq, lat.dt);

    end
  end
end

#! Perform collision over active regions of the domain
#!
#! \param   sim           Simulation
#! \param   active_cells  Active flags for domain
function call(col_f::MRT_F, sim::AbstractSim, active_cells::Matrix{Bool})
  lat           =   sim.lat;
  msm           =   sim.msm;
  const ni, nj  =   size(msm.rho);

  for j = 1:nj, i = 1:ni

    @inbounds if active_cells[i, j]

      @inbounds rhoij =   msm.rho[i,j];
      @inbounds uij   =   col_f.forcing_f[1](sim, i, j);
      feq             =   Vector{Float64}(lat.n);
      fneq            =   Vector{Float64}(lat.n);

      @simd for k = 1:lat.n 
        @inbounds feq[k]          =   col_f.feq_f(lat, rhoij, uij, k);
        @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
      end

      const mu    =   col_f.constit_relation_f(sim, fneq, col_f.S, col_f.M, 
                                               col_f.iM, i, j);
      const omega =   @omega(mu, lat.cssq, lat.dt);
      const Sij   =   col_f.S(mu, rhoij, lat.cssq, lat.dt);
      const F     =   map(k -> col_f.forcing_f[2](sim, omega, k, i, j), 1:lat.n);

      @inbounds lat.f[:,i,j]  -= col_f.iM * Sij * col_f.M * fneq - F;

      # update collision frequency matrix
      @inbounds msm.omega[i,j] = omega;

    end
  end
end

#! Perform collision active regions of the domain
#!
#! \param   sim           Simulation
#! \param   active_cells  Active flags for domain
function call(col_f::MRT, sim::FreeSurfSim, active_cells::Matrix{Bool})
  lat           =   sim.lat;
  msm           =   sim.msm;
  const ni, nj  =   size(msm.rho);

  for j = 1:nj, i = 1:ni

    @inbounds if active_cells[i, j] && sim.tracker.state[i, j] != GAS

      @inbounds rhoij =   msm.rho[i,j];
      @inbounds uij   =   msm.u[:,i,j];
      feq             =   Vector{Float64}(lat.n);
      fneq            =   Vector{Float64}(lat.n);

      @simd for k = 1:lat.n 
        @inbounds feq[k]          =   col_f.feq_f(lat, rhoij, uij, k);
        @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
      end

      const mu        =   col_f.constit_relation_f(sim, fneq, col_f.S, col_f.M, 
                                                   col_f.iM, i, j);
      const Sij       =   col_f.S(mu, rhoij, lat.cssq, lat.dt);

      @inbounds lat.f[:,i,j]  -= col_f.iM * Sij * col_f.M * fneq;

      # update collision frequency matrix
      @inbounds msm.omega[i,j] = @omega(mu, lat.cssq, lat.dt);

    end

  end
end

#! Perform collision active regions of the domain
#!
#! \param   sim           Simulation
#! \param   active_cells  Active flags for domain
function call(col_f::MRT_F, sim::FreeSurfSim, active_cells::Matrix{Bool})
  lat           =   sim.lat;
  msm           =   sim.msm;
  const ni, nj  =   size(msm.rho);

  for j = 1:nj, i = 1:ni

    @inbounds if active_cells[i, j] && sim.tracker.state[i, j] != GAS

      @inbounds rhoij =   msm.rho[i,j];
      @inbounds uij   =   col_f.forcing_f[1](sim, i, j);
      feq             =   Vector{Float64}(lat.n);
      fneq            =   Vector{Float64}(lat.n);

      @simd for k = 1:lat.n 
        @inbounds feq[k]          =   col_f.feq_f(lat, rhoij, uij, k);
        @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
      end

      const mu    =   col_f.constit_relation_f(sim, fneq, col_f.S, col_f.M, 
                                               col_f.iM, i, j);
      const omega =   @omega(mu, lat.cssq, lat.dt);
      const Sij   =   col_f.S(mu, rhoij, lat.cssq, lat.dt);
      const F     =   map(k -> col_f.forcing_f[2](sim, omega, k, i, j), 1:lat.n);

      @inbounds lat.f[:,i,j]  -= col_f.iM * Sij * col_f.M * fneq - F;

      # update collision frequency matrix
      @inbounds msm.omega[i,j] = omega;

    end
  end
end
