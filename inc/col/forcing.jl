# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

abstract Force;

type ConstForce <: Force
  fs::Tuple{LBXFunction, LBXFunction};

  ConstForce(f1::LBXFunction, f2::LBXFunction) = new((f1, f2));
end

Base.getindex(f::ConstForce, idx) = f.fs[idx];

# Kernal functions for Sukop forcing
function _sukop_f1(sim::AbstractSim, i::Int, j::Int)
  return sim.msm.u[:, i, j];
end
function _sukop_f2(sim::AbstractSim, omega::Real, k::Int, i::Int, j::Int,
                   F::Vector{Float64})
  return sim.lat.w[k] * sim.lat.dt / sim.lat.cssq * dot(F, sim.lat.c[:, k]);
end

#! Initialize a Sukop forcing function
#!
#! \param F Body force vector
#! \return (momentum_function, forcing_function)
function init_sukop_Fk(F::Vector{Float64}, ftype::Symbol=:ConstForce)
  return eval(:($ftype(
    _sukop_f1,
    @anon (sim, omega, k, i, j) -> _sukop_f2(sim, omega, k, i, j, $F)
    )));
end

# Function aliases
init_korner_Fk    =   init_sukop_Fk;
init_zhang_Fk     =   init_sukop_Fk;

#! Sukop forcing function but with gravitation acceleration instead of force
function _sukop_gravity(sim::AbstractSim, k::Int, i::Int, j::Int, 
                        g::Vector{Float64})
  return (sim.lat.w[k] * sim.lat.dt / sim.lat.cssq * sim.msm.rho[i, j] * 
          dot(g, sim.lat.c[:, k]));
end

#TODO should this just be M? or will M not be the actual mass for fluid cells?
#! Sukop forcing function but with gravitation acceleration instead of force
function _sukop_gravity(sim::FreeSurfSim, k::Int, i::Int, j::Int, 
                        g::Vector{Float64})
  return (sim.lat.w[k] * sim.lat.dt / sim.lat.cssq * sim.msm.rho[i, j] * 
          sim.tracker.eps[i, j] * dot(g, sim.lat.c[:, k]));
end

#! Initialize a Sukop forcing function with gravitation acceleration
#!
#! \param g Gravitation acceleration
#! \return (momentum_function, forcing_function)
function init_sukop_gravity_Fk(g::Vector{Float64}, ftype::Symbol=:ConstForce)
  return eval(:($ftype(
    _sukop_f1, 
    @anon (sim, omega, k, i, j) -> _sukop_gravity(sim, k, i, j, $g)
    )));
end


#TODO should this just be M? or will M not be the actual mass for fluid cells?
#! Velocity coupled gravitational acceleration
function _vel_coup_gravity(sim::AbstractSim, k::Int, i::Int, j::Int, 
                           g::Vector{Float64})
  const ck = sub(sim.lat.c ,:, k);
  const u  = sub(sim.msm.u, :, i, j);
  return (sim.lat.w[k] * sim.lat.dt / sim.lat.cssq * sim.msm.rho[i, j] * 
          dot(g, (ck - u) + dot(ck, u)*ck/sim.lat.cssq));
end

#! Velocity coupled gravitational acceleration
function _vel_coup_gravity(sim::FreeSurfSim, k::Int, i::Int, j::Int, 
                           g::Vector{Float64})
  const ck = sub(sim.lat.c ,:, k);
  const u  = sub(sim.msm.u, :, i, j);
  return (sim.tracker.eps[i, j] * sim.lat.w[k] * sim.lat.dt / sim.lat.cssq 
          * sim.msm.rho[i, j] * 
          dot(g, (ck - u) + dot(ck, u)*ck/sim.lat.cssq));
end

#! Initialize a velocity coupled body forcing function
#!
#! \param g Gravitation acceleration
#! \return (momentum_function, forcing_function)
function init_vel_coup_gravity_Fk(g::Vector{Float64}, ftype::Symbol=:ConstForce)
  return eval(:($ftype(
    _sukop_f1,
    @anon (sim, omega, k, i, j) -> _vel_coup_gravity(sim, k, i, j, $g)
    )));
end

_tp_vec(cssq) = [1/3; 1/3; 1/3; 1/3; 1/12; 1/12; 1/12; 1/12; 1 - 5/3 * cssq];


#! Initialize a Ginzburg and Steiner body forcing function
#!
#! \param g Gravitation acceleration
#! \return (momentum_function, forcing_function)
function init_gs_Fk(g::Vector{Float64}, ftype::Symbol=:ConstForce)
  return eval(:($ftype(
          (sim, i, j) -> if sim.msm.rho[i, j] != 0.0
                           sim.msm.u[:, i, j] + 0.5 * $g / sim.msm.rho[i, j]
                         else
                           sim.msm.u[:, i, j] + 0.5 * $g
                         end,
          (sim, omega, k, i, j) -> (_tp_vec(sim.lat.cssq)[k] * 
                                    dot(sub(sim.lat.c, :, k), $g))
    )));
end

type ScalableForce <: Force
  fs::Tuple{LBXFunction, LBXFunction};
  c::Real;

  ScalableForce(f1::LBXFunction, f2::LBXFunction) = new((f1, f2), 1.0);
  ScalableForce(f1::LBXFunction, f2::LBXFunction, c::Real) = new((f1, f2), c);
end

Base.getindex(f::ScalableForce, idx) = (args...) -> f.c * f.fs[idx](args...);
