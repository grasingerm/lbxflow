const __tracking_root__ = dirname(@__FILE__);
require(abspath(joinpath(__tracking_root__, "multiscale.jl")));

immutable Gas; end;         const GAS = Gas();
immutable Interface; end;   const INTERFACE = Interface();
immutable Fluid; end;       const FLUID = Fluid();

immutable Tracker
  state::Matrix{Union(Gas, Interface, Fluid)};
  M::Matrix{Float64};
  interface::Vector{(Int64, Int64)};

  Tracker(ni::Int, nj::Int, state::Union(Gas, Interface, Fluid) = Gas()) = 
    new(fill(state, (ni, nj)), zeros(ni, nj));

  Tracker(msm::MultiscaleMap, state::Union(Gas, Interface, Fluid) = Fluid()) = 
    new(fill(state, size(msm.rho)), copy(msm.rho));

  Tracker(t::Tracker) = new(copy(t.state), copy(t.M));
end

#! Update cell states of tracker
#!
#! \param t Mass tracker
#! \param msm Multiscale map of the domain
function update!(t::Tracker, msm::MultiscaleMap)
  const ni, nj = size(t.M);

  for 
    updatecell!(t, t.state[i,j], msm.rho[i,j], 

