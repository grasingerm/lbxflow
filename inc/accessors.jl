# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Velocity profile accessor
#!
#! \param     c           Component index
#! \param     i           ith index
#! \param     j_range     Range of j values
#! \return                Accessor function
function vel_prof_acsr(c::Int, i::Int, j_range::UnitRange{Int})
  @assert(c <= 3, "Component should be less than or equal to 3"); 
  return (sim::AbstractSim) -> begin
    const n = length(j_range);
    const x = linspace(-0.5, 0.5, n);

    f = @anon j -> sim.msm.u[c, i, j];
    const y = pmap(f, j_range);

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
    const n = length(i_range);
    const x = linspace(-0.5, 0.5, n);

    f = @anon i -> sim.msm.u[c, i, j];
    const y = pmap(f, j_range);

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
    const n = length(j_range);
    const x = linspace(-0.5, 0.5, n);

    f = @anon j -> sim.msm.u[c, i, j];
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
    const n = length(i_range);
    const x = linspace(-0.5, 0.5, n);

    f = @anon i -> sim.msm.u[c, i, j];
    y = pmap(f, j_range);
    y /= maximum(y);

    return x, y;
  end
end

#! Velocity magnitude accessor
#!
#! \param   sim   Simulation object
#! \return        Velocity magnitude over domain
function vel_mag_ascr(sim::AbstractSim)
  return transpose(u_mag(sim.msm));
end

#! Velocity field accessor
#!
#! \param   sim   Simulation object
#! \return        Velocity magnitude over domain
function vel_field_ascr(sim::AbstractSim)
  const ni, nj = size(sim.msm.rho);
  return (transpose(reshape(sub(sim.msm.u, 1, :, :), (ni, nj))), 
          transpose(reshape(sub(sim.msm.u, 2, :, :), (ni, nj))));
end

#! Density field accessor
#!
#! \param   sim   Simulation object
#! \return        Fluid density over domain
function density_ascr(sim::AbstractSim)
  return transpose(sim.msm.rho);
end

#! Pressure field accessor
#!
#! \param   sim   Simulation object
#! \return        Fluid pressure over domain
function pressure_ascr(sim::AbstractSim)
  f = @anon rho -> sim.lat.cssq * rho;
  return transpose(map(f, sim.msm.rho));
end

#! Mass field accessor
#!
#! \param   sim   Simulation object
#! \return        Fluid mass over domain
function mass_ascr(sim::FreeSurfSim)
  return transpose(sim.tracker.M);
end

#! Streamline fields accessor
#!
#! \param   sim   Simulation object
#! \return        x, y, u, v for streamlines
function streamlines_ascr(sim::AbstractSim)
  const ni, nj = size(sim.msm.rho);
  return (collect(linspace(0.0, 1.0, ni)), collect(linspace(0.0, 1.0, nj)),
          transpose(reshape(sub(sim.msm.u, 1, :, :), (ni, nj))), 
          transpose(reshape(sub(sim.msm.u, 2, :, :), (ni, nj))));
end
