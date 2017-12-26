# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Check to see if a pair of indices fall inside a bounds
function inbounds(i::Int, j::Int, sbounds::Matrix{Int64})
  const nbounds = size(sbounds, 2);
  for n=1:nbounds
    if !inbounds(i, j, sbounds[:, n]); return false; end
  end
  return true;
end

#! Check to see if a pair of indices fall inside a bounds
function inbounds(i::Int, j::Int, sbound::Vector{Int64})
  if i < sbound[1]; return false; end
  if i > sbound[2]; return false; end
  if j < sbound[3]; return false; end
  if j > sbound[4]; return false; end
  return true;
end

# Setting up types that will be used throughout module
#using FastAnonymous;
macro anon(expr)
  return expr;
end # empty macro

abstract ColFunction;
abstract FltrColFunction <: ColFunction;
type _ConstConstit
  μ::Real
  _ConstConstit(μ::Real) = new(μ);
end
#typealias LBXFunction Union{Function, FastAnonymous.AbstractClosure, ColFunction, _ConstConstit};
typealias LBXFunction Union{Function, ColFunction, _ConstConstit};

abstract AbstractSim;

#! Simulation object
#! lat Lattice
#! msm Multiscale map
immutable Sim <: AbstractSim
  lat::Lattice;
  msm::MultiscaleMap;
  Δt::Real;

  Sim(lat::Lattice, msm::MultiscaleMap) = new(lat, msm, lat.dt);
end

#! Type system for cell states
immutable Gas; end;         const GAS = Gas();
immutable Interface; end;   const INTERFACE = Interface();
immutable Fluid; end;       const FLUID = Fluid();

#! Type alias for cell states
typealias State Union{Gas, Interface, Fluid};

#! Mass and state tracker
immutable Tracker
  state::Matrix{State};
  M::Matrix{Float64};
  eps::Matrix{Float64};
  interfacels::Set{Tuple{Int64, Int64}};

  function Tracker(ni::Int, nj::Int,
                   state::State = GAS) 
    return new(convert(Matrix{Union{Gas, Interface, Fluid}}, 
                       fill(state, (ni, nj))),
               zeros(ni, nj), zeros(ni, nj),
               Set{Tuple{Int64, Int64}}());
  end

  function Tracker(msm::MultiscaleMap,
                   state::Matrix{State})

    const ni, nj  = size(state);
    lst           = Set{Tuple{Int64, Int64}}();
    M             = Array{Float64}(ni, nj);
    eps           = Array{Float64}(ni, nj);

    for j=1:nj, i=1:ni
      if state[i,j] == GAS
        M[i,j]    = 0.0;
        eps[i,j]  = 0.0; 
      elseif state[i,j] == INTERFACE
        M[i,j]    = 0.5 * msm.rho[i,j];
        eps[i,j]  = 0.5;
        push!(lst, (i,j));
      elseif state[i,j] == FLUID
        M[i,j]    = msm.rho[i,j];
        eps[i,j]  = 1.0;
      else
        error("state not understood or invalid");
      end
    end

    return new(copy(state), M, eps, lst);
  end
end

macro _safe_calc_eps(mass, rho)
  return quote
    if $mass != 0.0 && $rho != 0.0
      $mass / $rho;
    else
      0.0;
    end
  end
end

#! Free surface flow simulation object
#! lat Lattice
#! msm Multiscal map
#! tracker Mass tracker
immutable FreeSurfSim <: AbstractSim
  lat::Lattice
  msm::MultiscaleMap
  tracker::Tracker
  rho_g::AbstractFloat
  Δt::Real;

  function FreeSurfSim(lat::Lattice, msm::MultiscaleMap, tracker::Tracker)
    return new(lat, msm, tracker, 1.0, lat.dt);
  end

  FreeSurfSim(lat::Lattice, msm::MultiscaleMap, t::Tracker,
              rho_g::Real) = new(lat, msm, t, rho_g, lat.dt);

  function FreeSurfSim(lat::Lattice, msm::MultiscaleMap, rho_0::Real, 
                       rho_g::Real, fill_x::Real, fill_y::Real)
    const ni, nj    =     size(msm.rho);
    const fill_ni   =     convert(Int, round(fill_x * ni));
    const fill_nj   =     convert(Int, round(fill_y * nj));
    const nk        =     length(lat.w);

    t               =     Tracker(ni, nj, GAS);

    #TODO consider initializing f with distribution functions that satisfy
    #     conservation of linear momentum at the interface
    for j=1:nj, i=1:ni, k=1:nk
      lat.f[k, i, j]  =  rho_0 * lat.w[k];
    end

    map_to_macro!(lat, msm);

    for j=1:fill_nj, i=1:fill_ni
      t.M[i, j]     =   msm.rho[i, j];
      t.eps[i, j]   =   1.0;
      t.state[i, j] =   FLUID;
    end

    if fill_ni < ni
      for j=1:fill_nj
        const i       =   fill_ni + 1;
        t.M[i, j]     =   msm.rho[i, j] / 2.0;
        t.eps[i, j]   =   0.5;
        t.state[i, j] =   INTERFACE;
        push!(t.interfacels, (i, j));
      end
    end

    if fill_nj < nj
      for i=1:fill_ni
        const j       =   fill_nj + 1;
        t.M[i, j]     =   msm.rho[i, j] / 2.0;
        t.eps[i, j]   =   0.5;
        t.state[i, j] =   INTERFACE;
        push!(t.interfacels, (i, j));
      end
    end

    if fill_ni < ni && fill_nj < nj
      const i, j    =   fill_ni + 1, fill_nj + 1;
      t.M[i, j]     =   msm.rho[i, j] / 2.0;
      t.eps[i, j]   =   0.5;
      t.state[i, j] =   INTERFACE;
      push!(t.interfacels, (i, j));
    end

    return new(lat, msm, t, rho_g, lat.dt);
  end

end

#! Type for simulations with adaptive time stepping
type AdaptiveTimeStepSim <: AbstractSim
  lat::Lattice
  msm::MultiscaleMap
  isim::AbstractSim;
  ξ::Real;
  Δt::Real;
  incr::Bool;
  decr::Bool;
  relax::Real;

  function AdaptiveTimeStepSim(isim::AbstractSim, ξ::Real=4/5; incr::Bool=true, 
                               decr::Bool=true, relax=1.0)
    return new(isim.lat, isim.msm, isim, ξ, isim.lat.dt, incr, decr, relax);
  end
end

type M2PhaseSim <: AbstractSim
  simr::Sim;
  simb::Sim;
  Ar::Real;
  Ab::Real;
  αr::Real;
  αb::Real;
  β::Real;
  Δt::Real;

  function M2PhaseSim(nur::Real, nub::Real, rho_0r::Real, rho_0b::Real, 
                      ni::Int, nj::Int, Ar::Real, Ab::Real, αr::Real, αb::Real,
                      β::Real;
                      fill_r = (0.0, 1.0, 0.0, 1.0),
                      fill_b = (0.0, 1.0, 0.0, 1.0))
    @assert(β >= 0 && β <= 1, "β must be between 0 and 1. $(β) is invalid."); 

    latr = LatticeD2Q9(1.0, 1.0, ni, nj);
    latb = LatticeD2Q9(1.0, 1.0, ni, nj);

    _fill_to_range(percents, nx) = begin
      start = convert(Int, round(percents[1] * nx)) + 1;
      fin   = convert(Int, round(percents[2] * nx));
      return start:fin;
    end
    
    i_range_r = _fill_to_range((fill_r[1], fill_r[2]), ni);
    j_range_r = _fill_to_range((fill_r[3], fill_r[4]), nj);
    i_range_b = _fill_to_range((fill_b[1], fill_b[2]), ni);
    j_range_b = _fill_to_range((fill_b[3], fill_b[4]), nj);

    _fill_lat(latr, i_range_r, j_range_r, rho_0r);
    _fill_lat(latb, i_range_b, j_range_b, rho_0b);

    msmr = MultiscaleMap(nur, latr, rho_0r); 
    msmb = MultiscaleMap(nub, latb, rho_0b);

    map_to_macro!(latr, msmr);
    map_to_macro!(latb, msmb);

    return new(Sim(latr, msmr), Sim(latb, msmb), Ar, Ab, αr, αb, β, 1.0);
  end
end
