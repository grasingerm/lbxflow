# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Equilibrium frequency distribution for incompressible Newtonian flow
#!
#! \param lat D2Q9 lattice
#! \param msm Multiscale map
#! \param u Force modified velocity
#! \param i Lattice site x-direction index
#! \param j Lattice site y-direction index
#! \param k Lattice velocity vector index
#! \return Equilibrium frequency
function feq_incomp(lat::LatticeD2Q9, msm::MultiscaleMap, 
                    u::AbstractVector{Float64}, i::Int, j::Int, k::Int)
  @inbounds const rho = msm.rho[i, j];
  const cssq = 1/3;
  @inbounds const ckdotu = dot(view(lat.c, :, k), u);

  return rho * lat.w[k] * (1.0 + ckdotu/(cssq) + 0.5*(ckdotu*ckdotu)/(cssq*cssq)
                           - 0.5 * dot(u, u) / (cssq));
end

#! Equilibrium frequency distribution for incompressible Newtonian flow
#!
#! \param lat D2Q9 lattice
#! \param msm Multiscale map
#! \param u Force modified velocity
#! \param rho_0 Nominal or average density
#! \param i Lattice site x-direction index
#! \param j Lattice site y-direction index
#! \param k Lattice velocity vector index
#! \return Equilibrium frequency
function feq_incomp_HL(lat::LatticeD2Q9, msm::MultiscaleMap, 
                       u::AbstractVector{Float64}, rho_0::Real,
                       i::Int, j::Int, k::Int)
  @inbounds const rho = msm.rho[i, j];
  const cssq = 1/3;
  @inbounds const ckdotu = dot(view(lat.c, :, k), u);

  return (lat.w[k] * (rho + rho_0 * (ckdotu/(cssq)
              + 0.5*(ckdotu*ckdotu)/(cssq*cssq)
              - 0.5 * dot(u, u) / (cssq))));
end

#! Binds rho_0 to an HL equilibrium function
function init_feq_incomp_HL(rho_0::AbstractFloat)
  return (@anon (lat, msm, u, i, j, k) -> feq_incomp_HL(lat, msm, u, rho_0, i, j, k));
end

#! Equilibrium frequency distribution for incompressible Newtonian flow
#!
#! \param lat D2Q4 lattice
#! \param msm Multiscale map
#! \param u Force modified velocity
#! \param i Lattice site x-direction index
#! \param j Lattice site y-direction index
#! \param k Lattice velocity vector index
#! \return Equilibrium frequency
function feq_incomp(lat::LatticeD2Q4, msm::MultiscaleMap, 
                    u::AbstractVector{Float64}, i::Int, j::Int, k::Int)
  @inbounds const rho = msm.rho[i, j];
  const cssq = 1/2;
  @inbounds const ckdotu = dot(view(lat.c, :, k), u);

  return rho * lat.w[k] * (1 + 3 * ckdotu);
end

#! Equilibrium frequency distribution for incompressible Newtonian flow
#!
#! \param lat D2Q4 lattice
#! \param msm Multiscale map
#! \param u Force modified velocity
#! \param rho_0 Nominal or average density
#! \param i Lattice site x-direction index
#! \param j Lattice site y-direction index
#! \param k Lattice velocity vector index
#! \return Equilibrium frequency
function feq_incomp_HL(lat::LatticeD2Q4, msm::MultiscaleMap,
                       u::AbstractVector{Float64}, rho_0::AbstractFloat, i::Int, 
                       j::Int, k::Int)
  @inbounds const rho = msm.rho[i, j];
  const cssq = 1/2;
  @inbounds const ckdotu = dot(view(lat.c, :, k), u);

  return lat.w[k] * (rho + rho_0 * 3 * ckdotu);
end

#! Equilibrium frequency distribution that maximizes entropy
#! (Gorban and Packwood, Physica A 2014)
#!
#! \param lat D2Q9 lattice
#! \param msm Multiscale map
#! \param u Force modified velocity
#! \param i Lattice site x-direction index
#! \param j Lattice site y-direction index
#! \param k Lattice velocity vector index
#! \return Equilibrium frequency distribution that maximizes entropy
function feq_incomp_max_entropy(lat::LatticeD2Q9, msm::MultiscaleMap,
                                u::AbstractVector{Float64}, i::Int, j::Int, 
                                k::Int)
  @inbounds const rho = msm.rho[i, j];
  const nj  =   length(u);
  prod      =   1.0;

  for j=1:nj
    x1      =   sqrt(1 + 3.0*u[j]^2);
    prod    *=  if u[j] != 1
                  (2 - x1) * ( (2*u[j] + x1) / (1 - u[j]) )^lat.c[j,k];
                else
                  1.0;
                end;
  end

  return lat.w[k] * rho * prod;
end

#! Equilibrium distribution for immisible two-phase flow
#!
#! \param lat D2Q9 lattice
#! \param msm Multiscale map
#! \param u Force modified velocity
#! \param i Lattice site x-direction index
#! \param j Lattice site y-direction index
#! \param k Lattice velocity vector index
#! \param α Free paramter related to surface tension
#! \return Equilibrium particle distribution
function feq_incomp_mphase_immis(lat::LatticeD2Q9, msm::MultiscaleMap,
                                 u::AbstractVector{Float64}, i::Int, j::Int, 
                                 k::Int, α::Real)
  @inbounds const rho = msm.rho[i, j];
  @inbounds const ckdotu = dot(view(lat.c, :, k), u);

  if k == 9
    return rho * (α - 2/3 * dot(u, u));
  elseif k > 0 && k < 5
    return (rho * ((1-α)/5 + lat.w[k]*
                   (3*ckdotu + 9/2*ckdotu*ckdotu - 3/2 * dot(u, u))));
  elseif k > 4 && k < 9
    return (rho * ((1-α)/20 + lat.w[k]*
                   (3*ckdotu + 9/2*ckdotu*ckdotu - 3/2 * dot(u, u))));
  else
    error("k must be between 1 and 9 for D2Q9 lattice");
    return 0.0;
  end
end

#! Binds α to a two-phase immisible equilibrium function
function init_feq_mphase_immis(α::Real)
  return (@anon (lat, msm, u, i, j, k) -> feq_incomp_mphase_immis(lat, msm, u, i, j, k, α));
end

#! Equilibrium frequency distribution for incompressible Newtonian flow with smoothing
#!
#! \param lat D2Q9 lattice
#! \param msm Multiscale map
#! \param u Force modified velocity
#! \param i Lattice site x-direction index
#! \param j Lattice site y-direction index
#! \param k Lattice velocity vector index
#! \param alpha Weight for averaging
#! \return Equilibrium frequency
function feq_incomp_smoothing(lat::LatticeD2Q9, msm::MultiscaleMap, 
                              u::AbstractVector{Float64}, i::Int, j::Int, k::Int;
                              alpha::Real = 1.0)
  @inbounds const rho = msm.rho[i, j];
  const cssq = 1/3;

  const ni, nj = size(msm.rho);
  nbr_idxs = Tuple{Int, Int}[];
  if i > 1
    push!(nbr_idxs, (i-1, j));
    if i < ni
      push!(nbr_idxs, (i+1, j));
    end
  else
    push!(nbr_idxs, (i+1, j));
  end
  if j > 1
    push!(nbr_idxs, (i, j-1));
    if j < nj
      push!(nbr_idxs, (i, j+1));
    end
  else
    push!(nbr_idxs, (i, j+1));
  end

  const nα = (1.0 - alpha) / length(nbr_idxs);
  u_avg = alpha * u;
  rho_avg = alpha * rho;
  for nbr_idx in nbr_idxs
    @inbounds u_avg += nα * view(msm.u, :, nbr_idx[1], nbr_idx[2]);
    @inbounds rho_avg += nα * msm.rho[nbr_idx[1], nbr_idx[2]];
  end
  @inbounds const ckdotu = dot(view(lat.c, :, k), u_avg);

  return rho_avg * lat.w[k] * (1.0 + ckdotu/(cssq) + 0.5*(ckdotu*ckdotu)/(cssq*cssq)
                               - 0.5 * dot(u, u) / (cssq));
end

#! Binds α to smoothing equilibrium function
function init_feq_incomp_smoothing(α::Real)
  return (@anon (lat, msm, u, i, j, k) -> feq_incomp_smoothing(lat, msm, u, i, j, k; alpha=α));
end
