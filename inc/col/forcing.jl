# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

typealias Force Tuple{LBXFunction, LBXFunction};

#! Initialize a Sukop forcing function
#!
#! \param F Body force vector
#! \return (momentum_function, forcing_function)
function init_sukop_Fk(F::Vector{Float64})
  return (
    (sim::AbstractSim, i::Int, j::Int) -> begin
      return sim.msm.u[:, i, j];
    end,
    (sim::AbstractSim, omega::AbstractFloat, k::Int, i::Int, j::Int) -> begin
      return sim.lat.w[k] * sim.lat.dt / sim.lat.cssq * dot(F, sim.lat.c[:, k]);
    end
    );
end

# Function aliases
init_korner_Fk    =   init_sukop_Fk;
init_zhang_Fk     =   init_sukop_Fk;

#TODO implement Silva, Guo, Buick forcing functions...
#=
#! Initialize a Silva forcing function
#!
#! \param F Body force vector
#! \return (momentum_function, forcing_function)
function init_silva_Fk(F::Vector{Float64})
  return (
    (sim::AbstractSim, i::Int, j::Int) -> begin
      return sim.msm.u[:, i, j] + sim.lat.dt / 2.0 * F / sim.msm.rho[i, j];
    end,
    (sim::AbstractSim, omega::AbstractFloat, k::Int, i::Int, j::Int) -> begin
      return ((1 - 0.5 * omega) * sim.lat.w[k] * 
              dot(((sim.lat.c[:, k] - sim.msm.u[:, i, j]) / sim.lat.cssq +
              dot(sim.lat.c[:, k], sim.msm.u[:, i, j]) / (sim.lat.cssq * sim.lat.cssq) * sim.lat.c[:, k]), F));
    end
    );
end
=#

#! Sukop forcing function but with gravitation acceleration instead of force
function _sukop_gravity(g::Vector{Float64}, sim::AbstractSim, 
                        k::Int, i::Int, j::Int)
  return (sim.lat.w[k] * sim.lat.dt / sim.lat.cssq * sim.msm.rho[i, j] * 
          dot(g, sim.lat.c[:, k]));
end

#TODO should this just be M? or will M not be the actual mass for fluid cells?
#! Sukop forcing function but with gravitation acceleration instead of force
function _sukop_gravity(g::Vector{Float64}, sim::FreeSurfSim, 
                        k::Int, i::Int, j::Int)
  return (sim.lat.w[k] * sim.lat.dt / sim.lat.cssq * sim.msm.rho[i, j] * 
          sim.tracker.eps[i, j] * dot(g, sim.lat.c[:, k]));
end

#! Initialize a Sukop forcing function with gravitation acceleration
#!
#! \param g Gravitation acceleration
#! \return (momentum_function, forcing_function)
function init_sukop_gravity_Fk(g::Vector{Float64})
  return (
    (sim::AbstractSim, i::Int, j::Int) -> begin
      return sim.msm.u[:, i, j];
    end,
    (sim::AbstractSim, omega::AbstractFloat, k::Int, i::Int, j::Int) -> begin
      return _sukop_gravity(g, sim, k, i, j);
    end
    );
end
