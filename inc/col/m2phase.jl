# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Collision function wrapper for two-phase flow
type M2PhaseColFunction <: ColFunction
  col_fr!::ColFunction;
  col_fb!::ColFunction;
  m2phase_col_f!::LBXFunction;

  M2PhaseColFunction(col_fr!::ColFunction, col_fr!::ColFunction, 
                     mcol_f!::LBXFunction=m2phase_col_f) = new(col_fr!, col_fb!, 
                                                               mcol_f!);
end

#! Foward calls for collision functions
function call(col_f::M2PhaseColFunction, sim::M2PhaseSim, args...)
  col_f.col_fr!(sim.simr, args...);
  col_f.col_fb!(sim.simb, args...);
  col_f.m2phase_col_f!(sim, args...);
end

#! Calculate color gradient
function _color_grad(sim::M2PhaseSim, cbounds::Matrix{Int64}, i::Int, j::Int)
  F = Array{Float64, 2}(2);
  for k=1:sim.simr.lat.n-1
    const i_nbr   =   i + sim.rsim.lat.c[1, k];
    const j_nbr   =   j + sim.rsim.lat.c[2, k];

    if inbounds(i_nbr, j_nbr, cbounds)
      F += sub(sim.rsim.lat.c, :, k) * (sim.simr.msm.rho[i_nbr, j_nbr] - 
                                        sim.simb.msm.rho[i_nbr, j_nbr]);
    end
  end
  return F;
end

#! Calculate color gradient
function _color_grad(sim::M2PhaseSim, active_cells::Matrix{Bool}, i::Int, j::Int)
  F = Array{Float64, 2}(2);
  for k=1:sim.simr.lat.n-1
    const i_nbr   =   i + sim.rsim.lat.c[1, k];
    const j_nbr   =   j + sim.rsim.lat.c[2, k];

    if active_cells[i_nbr, j_nbr]
      F += sub(sim.rsim.lat.c, :, k) * (sim.simr.msm.rho[i_nbr, j_nbr] - 
                                        sim.simb.msm.rho[i_nbr, j_nbr]);
    end
  end
  return F;
end

const __B = Float64[2/27, 2/27, 2/27, 2/27, 5/108, 5/108, 5/108, 5/108, -4/27];

#! Two-phase flow collision function
function m2phase_col_f!(sim::M2PhaseSim, cbounds::Matrix{Int64})
  const nbounds = size(cbounds, 2);
  for r = 1:nbounds
    i_min, i_max, j_min, j_max = sub(cbounds, :, r);
    for j=j_min:j_max, i=i_min:i_max
      const F = _color_grad(sim, cbounds, i, j);
      for k=1:sim.simr.lat.n
        Fdotcsq   =   (dot(F, sub(sim.simr.lat.c, :, k)))^2;
        Fmag      =   norm(F, 2);
        sim.simr.lat.f[k, i, j] += (sim.simr.lat.w[k] * sim.simr.Ar/2 * Fmag * 
                                    (Fdotcsq/(Fmag*Fmag) - __B[k]));
        sim.simb.lat.f[k, i, j] += (sim.simb.lat.w[k] * sim.simb.Ab/2 * Fmag * 
                                    (Fdotcsq/(Fmag*Fmag) - __B[k]));
      end
    end
  end
end

#! Two-phase flow collision function
function m2phase_col_f!(sim::M2PhaseSim, active_cells::Matrix{Bool})
  const ni, nj = size(sim.simr.msm.rho);
  for j=1:nj, i=1:ni
    if active_cells[i, j]
      const F = _color_grad(sim, active_cells, i, j);
      for k=1:sim.simr.lat.n
        Fdotcsq   =   (dot(F, sub(sim.simr.lat.c, :, k)))^2;
        Fmag      =   norm(F, 2);
        sim.simr.lat.f[k, i, j] += (sim.simr.lat.w[k] * sim.simr.Ar/2 * Fmag * 
                                    (Fdotcsq/(Fmag*Fmag) - __B[k]));
        sim.simb.lat.f[k, i, j] += (sim.simb.lat.w[k] * sim.simb.Ab/2 * Fmag * 
                                    (Fdotcsq/(Fmag*Fmag) - __B[k]));
      end
    end
  end
end
