# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Distributed;

#! Kernal function for accessing velocities
function _vel_acsr_kernal(sim::AbstractSim, c::Int, i::Int, j::Int)
  return sim.msm.u[c, i, j];
end

#! Kernal function for accessing velocities
function _vel_acsr_kernal(sim::AdaptiveTimeStepSim, c::Int, i::Int, j::Int)
  return sim.isim.msm.u[c, i, j] / sim.Δt;
end

#! Velocity profile accessor
#!
#! \param     c           Component index
#! \param     i           ith index
#! \param     j_range     Range of j values
#! \return                Accessor function
function vel_prof_acsr(c::Int, i::Int, j_range::UnitRange{Int})
  @assert(c <= 3, "Component should be less than or equal to 3"); 
  return (sim::AbstractSim) -> begin
    n = length(j_range);
    x = range(-0.5, 0.5; length=n);

    f = j -> _vel_acsr_kernal(sim, c, i, j);
    y = pmap(f, j_range);

    return x, y;
  end
end

#! Velocity profile accessor
#!
#! \param     c           Component index
#! \param     i_range     Range of i values
#! \param     j           jth index
#! \return                Accessor function
function vel_prof_acsr(c::Int, i_range::UnitRange{Int}, j::Int)
  @assert(c <= 3, "Component should be less than or equal to 3"); 
  return (sim::AbstractSim) -> begin
    n = length(i_range);
    x = range(-0.5, 0.5; length=n);

    f = i -> _vel_acsr_kernal(sim, c, i, j);
    y = pmap(f, i_range);

    return x, y;
  end
end

#! Normalized velocity profile accessor
#!
#! \param     c           Component index
#! \param     i           ith index
#! \param     j_range     Range of j values
#! \return                Accessor function
function vbar_prof_acsr(c::Int, i::Int, j_range::UnitRange{Int})
  @assert(c <= 3, "Component should be less than or equal to 3"); 
  return (sim::AbstractSim) -> begin
    n = length(j_range);
    x = range(-0.5, 0.5; length=n);

    f = j -> _vel_acsr_kernal(sim, c, i, j);
    y = pmap(f, j_range);
    y /= maximum(y);

    return x, y;
  end
end

#! Normalized velocity profile accessor
#!
#! \param     c           Component index
#! \param     i_range     Range of i values
#! \param     j           jth index
#! \return                Accessor function
function vbar_prof_acsr(c::Int, i_range::UnitRange{Int}, j::Int)
  @assert(c <= 3, "Component should be less than or equal to 3"); 
  return (sim::AbstractSim) -> begin
    n = length(i_range);
    x = range(-0.5, 0.5; length=n);

    f = i -> _vel_acsr_kernal(sim, c, i, j);
    y = pmap(f, i_range);
    y /= maximum(y);

    return x, y;
  end
end

#! Velocity magnitude accessor
#!
#! \param   sim   Simulation object
#! \return        Velocity magnitude over domain
function vel_mag_acsr(sim::AbstractSim)
  return transpose(u_mag(sim.msm));
end

#! Velocity magnitude accessor
#!
#! \param   sim   Simulation object
#! \return        Velocity magnitude over domain
function vel_mag_acsr(sim::M2PhaseSim)
  ni, nj = size(sim.isim.msm);
  u_mag = zeros(ni, nj);
  for j=1:nj, i=1:ni
    ux = sim.simr.msm.u[1, i, j] + sim.simb.msm.u[1, i, j];
    uy = sim.simr.msm.u[2, i, j] + sim.simb.msm.u[2, i, j];
    u_mag[j, i] = sqrt(ux*ux + uy*uy);
  end
  return u_mag;
end

#! Velocity magnitude accessor
#!
#! \param   sim   Simulation object
#! \return        Velocity magnitude over domain
function vel_mag_acsr(sim::AdaptiveTimeStepSim)
  ni, nj = size(sim.isim.msm);
  u_mag = zeros(ni, nj);
  for j=1:nj, i=1:ni
    ux = sim.isim.msm.u[1, i, j] / sim.Δt;
    uy = sim.isim.msm.u[2, i, j] / sim.Δt;
    u_mag[j, i] = sqrt(ux*ux + uy*uy);
  end
  return u_mag;
end

#! Velocity field accessor
#!
#! \param   sim   Simulation object
#! \return        Velocity magnitude over domain
function vel_field_acsr(sim::AbstractSim)
  ni, nj = size(sim.msm.rho);
  return (transpose(reshape(view(sim.msm.u, 1, :, :), (ni, nj))), 
          transpose(reshape(view(sim.msm.u, 2, :, :), (ni, nj))));
end

#! Velocity field accessor
#!
#! \param   sim   Simulation object
#! \return        Velocity magnitude over domain
function vel_field_acsr(sim::AdaptiveTimeStepSim)
  ni, nj = size(sim.isim.msm.rho);
  return (transpose(map(u -> u / sim.Δt, reshape(view(sim.msm.u, 1, :, :)), 
                    (ni, nj))), 
          transpose(map(u -> u / sim.Δt, reshape(view(sim.msm.u, 2, :, :)), 
                    (ni, nj))));
end

#! Density field accessor
#!
#! \param   sim   Simulation object
#! \return        Fluid density over domain
function density_acsr(sim::AbstractSim)
  return transpose(sim.msm.rho);
end

#! Density field accessor
#!
#! \param   sim   Simulation object
#! \return        Fluid density over domain
function density_acsr(sim::M2PhaseSim)
  return transpose(sim.simr.msm.rho) + transpose(sim.simb.msm.rho);
end

#! Density field accessor
#!
#! \param   sim   Simulation object
#! \return        Fluid density over domain
function density_acsr(sim::AdaptiveTimeStepSim)
  return transpose(map(ρ -> ρ / sim.Δt, sim.isim.msm.rho));
end

#! Pressure field accessor
#!
#! \param   sim   Simulation object
#! \return        Fluid pressure over domain
function pressure_acsr(sim::AbstractSim)
  f = rho -> sim.lat.cssq * rho;
  return transpose(map(f, sim.msm.rho));
end

#! Pressure field accessor
#!
#! \param   sim   Simulation object
#! \return        Fluid pressure over domain
function pressure_acsr(sim::M2PhaseSim)
  f = rho -> sim.simr.lat.cssq * rho;
  return transpose(map(f, sim.simr.msm.rho)) + transpose(map(f, sim.simb.msm.rho));
end

#! Mass field accessor
#!
#! \param   sim   Simulation object
#! \return        Fluid mass over domain
function mass_acsr(sim::FreeSurfSim)
  return transpose(sim.tracker.M);
end

#! Mass field accessor
#!
#! \param   sim   Simulation object
#! \return        Fluid mass over domain
function mass_acsr(sim::AdaptiveTimeStepSim)
  return transpose(map(m -> m / sim.Δt, sim.tracker.M));
end

#! Fluid fraction accessor
#!
#! \param   sim   Simulation object
#! \return        Fluid fraction over domain
function ff_acsr(sim::FreeSurfSim)
  return transpose(sim.tracker.eps);
end

#! Streamline fields accessor
#!
#! \param   sim   Simulation object
#! \return        x, y, u, v for streamlines
function streamlines_acsr(sim::AbstractSim)
  ni, nj = size(sim.msm.rho);
  return (collect(range(0.0, 1.0; length=ni)), collect(range(0.0, 1.0; length=nj)),
          transpose(reshape(view(sim.msm.u, 1, :, :), (ni, nj))), 
          transpose(reshape(view(sim.msm.u, 2, :, :), (ni, nj))));
end

#! Streamline fields accessor
#!
#! \param   sim   Simulation object
#! \return        x, y, u, v for streamlines
function streamlines_acsr(sim::AdaptiveTimeStepSim)
  ni, nj = size(sim.isim.msm);
  u = zeros(ni, nj);
  v = zeros(ni, nj);
  for j=1:nj, i=1:ni
    u[j, i] = sim.isim.msm[1, i, j] / sim.Δt;
    v[j, i] = sim.isim.msm[2, i, j] / sim.Δt;
  end
  return (collect(range(0.0, 1.0; length=ni)), collect(range(0.0, 1.0; length=nj)),
          u, v);
end

#! Wrapper for two-phase simulations
#!
#! \param   inner_acsr    Inner accessor function to apply to color
#! \param   color         Color of fluid to access
#! \return                Annonymous wrapper
function m2phase_acsr(inner_acsr::LBXFunction, color::Symbol)
  if      color == :red
    return (sim::M2PhaseSim) -> inner_acsr(sim.simr); 
  elseif  color == :blue
    return (sim::M2PhaseSim) -> inner_acsr(sim.simb); 
  else
    error("Color $color is not understood. `red` or `blue` are valid colors.");
  end
end

#! Fluid fraction accessor for two-phase flow simulations
#!
#! \param   color   Color of fluid to access
#! \return          Accessor for fluid fractions
function fluid_frac_acsr(color::Symbol=:red)
  if      color == :red
    return (sim::M2PhaseSim) -> begin
      ni, nj  =  size(sim.simr.msm.rho);
      ff            =  zeros(ni, nj); 
      for j=1:nj, i=1:ni
        ρ_r, ρ_b      =  sim.simr.msm.rho[i, j], sim.simb.msm.rho[i, j]; 
        ff[i, j]      =  ρ_r / (ρ_r + ρ_b);
      end
      return transpose(ff);
    end
  elseif  color == :blue
    return (sim::M2PhaseSim) -> begin
      ni, nj  =  size(sim.simr.msm.rho);
      ff            =  zeros(ni, nj); 
      for j=1:nj, i=1:ni
        ρ_r, ρ_b      =  sim.simr.msm.rho[i, j], sim.simb.msm.rho[i, j]; 
        ff[i, j]      =  ρ_b / (ρ_r + ρ_b);
      end
      return transpose(ff);
    end 
  else
    error("Color $color is not understood. `red` or `blue` are valid colors.");
  end  
end

#! Initialize a default fluid fraction accessor for red fluids
red_fluid_frac_acsr = fluid_frac_acsr();
