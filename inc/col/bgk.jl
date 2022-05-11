# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#TODO should we simplify and not differentiate between collision function with
#     forcing and without?

#! Bhatnagar-Gross-Krook collision operator 
struct BGK <: ColFunction
  feq_f::LBXFunction;
  constit_relation_f::LBXFunction;

  BGK(constit_relation_f::LBXFunction; feq_f::LBXFunction=feq_incomp) = (
                                                new(feq_f, constit_relation_f));
end

#! Bhatnagar-Gross-Krook collision operator with forcing 
struct BGK_F <: ColFunction
  feq_f::LBXFunction;
  constit_relation_f::LBXFunction;
  forcing_f::Force;

  BGK_F(constit_relation_f::LBXFunction,
        forcing_f::Force; 
        feq_f::LBXFunction=feq_incomp) = (
                                    new(feq_f, constit_relation_f, forcing_f));
end

#TODO consider wrapping a lot of code in kernal functions to help JIT n'at

#! Perform collision over box regions of the domain
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function (col_f::BGK)(sim::AbstractSim, bounds::AbstractMatrix{Int64})
  lat           =   sim.lat;
  msm           =   sim.msm;
  ni, nj  =   size(msm.rho);
  nbounds =   size(bounds, 2);

  for r = 1:nbounds
    @inbounds i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max


      @inbounds uij   =   view(msm.u, :, i, j);
      feq             =   zeros(lat.n);
      fneq            =   zeros(lat.n);

      for k = 1:lat.n 
        @inbounds feq[k]          =   col_f.feq_f(lat, msm, uij, i, j, k);
        @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
      end

      mu        =   col_f.constit_relation_f(sim, fneq, i, j);
      omegaij   =   omega(mu, lat.cssq, lat.dt);

      for k = 1:lat.n
        @inbounds lat.f[k,i,j] = (omegaij * feq[k] + (1.0 - omegaij) * 
                                  lat.f[k,i,j]);
      end

      @inbounds msm.omega[i,j]  =   omegaij;

    end
  end
end

#! Perform collision over box regions of the domain
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function (col_f::BGK_F)(sim::AbstractSim, bounds::AbstractMatrix{Int64})
  lat           =   sim.lat;
  msm           =   sim.msm;
  ni, nj  =   size(msm.rho);
  nbounds =   size(bounds, 2);

  for r = 1:nbounds
    @inbounds i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max

      @inbounds uij   =   col_f.forcing_f[1](sim, i, j);
      feq             =   zeros(lat.n);
      fneq            =   zeros(lat.n);

      for k = 1:lat.n 
        @inbounds feq[k]          =   col_f.feq_f(lat, msm, uij, i, j, k);
        @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
      end

      mu        =   col_f.constit_relation_f(sim, fneq, i, j);
      omegaij   =   omega(mu, lat.cssq, lat.dt);

      for k = 1:lat.n
        @inbounds lat.f[k,i,j]    =   ((omegaij * feq[k] + (1.0 - omegaij) * 
                                        lat.f[k,i,j]) +
                                       col_f.forcing_f[2](sim, omegaij, k, i, j));
      end

      @inbounds msm.omega[i,j]  =   omegaij;

    end
  end
end

#! Perform collision over box regions of the domain for free surface flow
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function (col_f::BGK)(sim::FreeSurfSim, bounds::AbstractMatrix{Int64})
  lat           =   sim.lat;
  msm           =   sim.msm;
  ni, nj  =   size(msm.rho);
  nbounds =   size(bounds, 2);

  for r = 1:nbounds
    @inbounds i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max

      @inbounds if sim.tracker.state[i, j] != GAS

        @inbounds uij   =   view(msm.u, :, i, j);
        feq             =   zeros(lat.n);
        fneq            =   zeros(lat.n);

        for k = 1:lat.n 
          @inbounds feq[k]          =   col_f.feq_f(lat, msm, uij, i, j, k);
          @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
        end

        mu        =   col_f.constit_relation_f(sim, fneq, i, j);
        omegaij   =   omega(mu, lat.cssq, lat.dt);

        for k = 1:lat.n
          @inbounds lat.f[k,i,j]    =   (omegaij * feq[k] + (1.0 - omegaij) * 
                                         lat.f[k,i,j]);
        end

        @inbounds msm.omega[i,j]  =   omegaij;

      end

    end
  end
end

#! Perform collision over box regions of the domain for free surface flow
#!
#! \param   sim     Simulation
#! \param   bounds  Each column of the matrix defines a box region
function (col_f::BGK_F)(sim::FreeSurfSim, bounds::AbstractMatrix{Int64})
  lat           =   sim.lat;
  msm           =   sim.msm;
  ni, nj  =   size(msm.rho);
  nbounds =   size(bounds, 2);

  for r = 1:nbounds
    @inbounds i_min, i_max, j_min, j_max = bounds[:,r];
    for j = j_min:j_max, i = i_min:i_max

      @inbounds if sim.tracker.state[i, j] != GAS

        @inbounds uij   =   col_f.forcing_f[1](sim, i, j);
        feq             =   zeros(lat.n);
        fneq            =   zeros(lat.n);

        for k = 1:lat.n 
          @inbounds feq[k]          =   col_f.feq_f(lat, msm, uij, i, j, k);
          @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
        end

        mu        =   col_f.constit_relation_f(sim, fneq, i, j);
        omegaij   =   omega(mu, lat.cssq, lat.dt);

        for k = 1:lat.n
          @inbounds lat.f[k,i,j]    =   ((omegaij * feq[k] + (1.0 - omegaij) * 
                                         lat.f[k,i,j]) +
                                         col_f.forcing_f[2](sim, omegaij, k, i, 
                                                            j));
        end

        @inbounds msm.omega[i,j]  =   omegaij;

      end

    end
  end
end

#! Perform collision over active parts of the domain
#!
#! \param   sim           Simulation
#! \param   active_cells  Active flags for domain
function (col_f::BGK)(sim::AbstractSim, active_cells::AbstractMatrix{Bool})
  lat           =   sim.lat;
  msm           =   sim.msm;
  ni, nj  =   size(msm.rho);

  for j = 1:nj, i = 1:ni
      
    @inbounds if active_cells[i, j]

      @inbounds uij   =   view(msm.u, :, i, j);
      feq             =   zeros(lat.n);
      fneq            =   zeros(lat.n);

      for k = 1:lat.n 
        @inbounds feq[k]          =   col_f.feq_f(lat, msm, uij, i, j, k);
        @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
      end

      mu        =   col_f.constit_relation_f(sim, fneq, i, j);
      omegaij   =   omega(mu, lat.cssq, lat.dt);

      for k = 1:lat.n
        @inbounds lat.f[k,i,j]    =   (omegaij * feq[k] + (1.0 - omegaij) * 
                                       lat.f[k,i,j]);
      end

      @inbounds msm.omega[i,j]  =   omegaij;

    end
  end
end

#! Perform collision over active parts of the domain
#!
#! \param   sim           Simulation
#! \param   active_cells  Active flags for domain
function (col_f::BGK_F)(sim::AbstractSim, active_cells::AbstractMatrix{Bool})
  lat           =   sim.lat;
  msm           =   sim.msm;
  ni, nj  =   size(msm.rho);

  for j = 1:nj, i = 1:ni
    
    @inbounds if active_cells[i, j]

      @inbounds uij   =   col_f.forcing_f[1](sim, i, j);
      feq             =   zeros(lat.n);
      fneq            =   zeros(lat.n);

      for k = 1:lat.n 
        @inbounds feq[k]          =   col_f.feq_f(lat, msm, uij, i, j, k);
        @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
      end

      mu        =   col_f.constit_relation_f(sim, fneq, i, j);
      omegaij   =   omega(mu, lat.cssq, lat.dt);

      for k = 1:lat.n
        @inbounds lat.f[k,i,j]    =   ((omegaij * feq[k] + (1.0 - omegaij) * 
                                       lat.f[k,i,j]) +
                                       col_f.forcing_f[2](sim, omegaij, k, i, j));
      end

      @inbounds msm.omega[i,j]  =   omegaij;

    end
  end
end

#! Perform collision over active parts of the domain for free surface flow
#!
#! \param   sim           Simulation
#! \param   active_cells  Active flags for domain
function (col_f::BGK)(sim::FreeSurfSim, active_cells::AbstractMatrix{Bool})
  lat           =   sim.lat;
  msm           =   sim.msm;
  ni, nj  =   size(msm.rho);

  for j = 1:nj, i = 1:ni

    @inbounds if active_cells[i, j] && sim.tracker.state[i, j] != GAS

      @inbounds uij   =   view(msm.u, :, i, j);
      feq             =   zeros(lat.n);
      fneq            =   zeros(lat.n);

      for k = 1:lat.n 
        @inbounds feq[k]          =   col_f.feq_f(lat, msm, uij, i, j, k);
        @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
      end

      mu        =   col_f.constit_relation_f(sim, fneq, i, j);
      omegaij     =   omega(mu, lat.cssq, lat.dt);

      for k = 1:lat.n
        @inbounds lat.f[k,i,j]    =   (omegaij * feq[k] + (1.0 - omegaij) * 
                                       lat.f[k,i,j]);
      end

      @inbounds msm.omega[i,j]  =   omegaij;

    end

  end
end

#! Perform collision over active parts of the domain for free surface flow
#!
#! \param   sim           Simulation
#! \param   active_cells  Active flags for domain
function (col_f::BGK_F)(sim::FreeSurfSim, active_cells::AbstractMatrix{Bool})
  lat           =   sim.lat;
  msm           =   sim.msm;
  ni, nj  =   size(msm.rho);

  for j = 1:nj, i = 1:ni

    @inbounds if active_cells[i, j] && sim.tracker.state[i, j] != GAS

      @inbounds uij   =   col_f.forcing_f[1](sim, i, j);
      feq             =   zeros(lat.n);
      fneq            =   zeros(lat.n);

      for k = 1:lat.n 
        @inbounds feq[k]          =   col_f.feq_f(lat, msm, uij, i, j, k);
        @inbounds fneq[k]         =   lat.f[k,i,j] - feq[k];
      end

      mu        =   col_f.constit_relation_f(sim, fneq, i, j);
      omegaij   =   omega(mu, lat.cssq, lat.dt);

      for k = 1:lat.n
        @inbounds lat.f[k,i,j]    =   ((omegaij * feq[k] + (1.0 - omegaij) * 
                                       lat.f[k,i,j]) +
                                       col_f.forcing_f[2](sim, omegaij, k, i, j));
      end

      @inbounds msm.omega[i,j]  =   omegaij;

    end

  end
end
