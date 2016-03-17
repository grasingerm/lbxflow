# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

using List

#using FastAnonymous;
#typealias LBXFunction Union{Function, FastAnonymous.AbstractClosure};
typealias LBXFunction Function;

abstract AbstractSim

#! Simulation object
#! lat Lattice
#! msm Multiscale map
immutable Sim <: AbstractSim
  lat::Lattice;
  msm::MultiscaleMap;

  Sim(lat::Lattice, msm::MultiscaleMap) = new(lat, msm);
end

#! Type system for cell states
immutable Gas; end;         const GAS = Gas();
immutable Interface; end;   const INTERFACE = Interface();
immutable Fluid; end;       const FLUID = Fluid();

#! Mass and state tracker
immutable Tracker
  state::Matrix{Union{Gas, Interface, Fluid}};
  M::Matrix{Float64};
  eps::Matrix{Float64};
  interfacels::DoublyLinkedList{Tuple{Int64, Int64}};

  function Tracker(state::Union{Gas, Interface, Fluid} = GAS) 
    const ni, nj = size(state);
    return new(fill(state, (ni, nj)), zeros(ni, nj), zeros(ni, nj),
               DoublyLinkedList{Tuple{Int64, Int64}}());
  end

  function Tracker(msm::MultiscaleMap,
                   state::Matrix{Union{Gas, Interface, Fluid}})

    const ni, nj  = size(state);
    lst           = DoublyLinkedList{Tuple{Int64, Int64}}();
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

#! Free surface flow simulation object
#! lat Lattice
#! msm Multiscal map
#! tracker Mass tracker
immutable FreeSurfSim <: AbstractSim
  lat::Lattice
  msm::MultiscaleMap
  tracker::Tracker
  rhog::AbstractFloat

  function FreeSurfSim(lat::Lattice, msm::MultiscaleMap, tracker::Tracker)
    return new(lat, msm, tracker, 1.0);
  end

  FreeSurfSim(lat::Lattice, msm::MultiscaleMap, t::Tracker,
              rhog::Real) = new(lat, msm, t, rhog);

  function FreeSurfSim(lat::Lattice, msm::MultiscaleMap, rho_0::Real, 
                       rho_g::Real, fill_x::Real, fill_y::Real)
    const ni, nj    =     size(msm.rho);
    const fill_ni   =     convert(Int, fill_x * ni);
    const fill_nj   =     convert(Int, fill_y * nj);
    const nk        =     length(lat.w);

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
      end
    end

    if fill_nj < nj
      for i=1:fill_ni
        const j       =   fill_nj + 1;
        t.M[i, j]     =   msm.rho[i, j] / 2.0;
        t.eps[i, j]   =   0.5;
        t.state[i, j] =   INTERFACE;
      end
    end

    if fill_ni < ni && fill_nj < nj
      const i, j    =   fill_ni + 1, fill_nj + 1;
      t.M[i, j]     =   msm.rho[i, j] / 2.0;
      t.eps[i, j]   =   0.5;
      t.state[i, j] =   INTERFACE;
    end
  end

end
