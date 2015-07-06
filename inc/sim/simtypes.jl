using List

abstract AbstractSim

#! Simulation object
#! lat Lattice
#! msm Multiscale map
immutable Sim <: AbstractSim
  lat::Lattice
  msm::MultiscaleMap

  Sim(lat::Lattice, msm::MultiscaleMap) = new(lat, msm);
end

#! Type system for cell states
immutable Gas; end;         const GAS = Gas();
immutable Interface; end;   const INTERFACE = Interface();
immutable Fluid; end;       const FLUID = Fluid();

#! Mass and state tracker
immutable Tracker
  state::Matrix{Union(Gas, Interface, Fluid)};
  M::Matrix{Float64};
  eps::Matrix{Float64};
  interfacels::DoublyLinkedList{(Int64, Int64)};

  function Tracker(state::Union(Gas, Interface, Fluid) = Gas()) 
    const ni, nj = size(state);
    return new(fill(state, (ni, nj)), zeros(ni, nj), zeros(ni, nj),
               DoublyLinkedList{(Int64, Int64)});
  end

  function Tracker(msm::MultiscaleMap,
                   state::Matrix{Union(Gas, Interface, Fluid)})
    const ni, nj = size(state);
    lst = DoublyLinkedList{(Int64,Int64)}();
    M = Array(Float64, (ni, nj));
    eps = Array(Float64, (ni, nj));
    for j=1:nj, i=1:ni
      if state[i,j] == GAS
        M[i,j] = 0.0;
        eps[i,j] = 0.0; 
      elseif state[i,j] == INTERFACE
        M[i,j] = 0.5 * msm.rho[i,j];
        eps[i,j] = 0.5;
        push!(lst, (i,j));
      elseif state[i,j] == FLUID
        M[i,j] = msm.rho[i,j];
        eps[i,j] = 1.0;
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
  rhog::FloatingPoint

  function FreeSurfSim(lat::Lattice, msm::MultiscaleMap, tracker::Tracker)
    return new(lat, msm, tracker, 1.0);
  end

  FreeSurfSim(lat::Lattice, msm::MultiscaleMap, t::Tracker,
              rhog::FloatingPoint) = new(lat, msm, t, rhog);
end
