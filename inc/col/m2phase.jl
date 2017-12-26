# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Collision function wrapper for two-phase flow
type M2PhaseColFunction <: ColFunction
  col_fr!::ColFunction;
  col_fb!::ColFunction;
  m2phase_col_f!::LBXFunction;

  M2PhaseColFunction(col_fr!::ColFunction, col_fb!::ColFunction, 
                     mcol_f!::LBXFunction=m2phase_col_f!) = new(col_fr!, col_fb!, 
                                                               mcol_f!);
end

#! Foward calls for collision functions
function (col_f::M2PhaseColFunction)(sim::M2PhaseSim, args...)
  println("Red collision function");
  col_f.col_fr!(sim.simr, args...);
  println("Is there an NaN?", true in isnan(sim.simr.lat.f));
  readline(STDIN);
  println("Blue collision function");
  col_f.col_fb!(sim.simb, args...);
  println("Is there an NaN?", true in isnan(sim.simb.lat.f));
  readline(STDIN);
  println("2 phase collision function");
  col_f.m2phase_col_f!(sim, args...);
  println("Is there an NaN?", true in isnan(sim.simr.lat.f));
  readline(STDIN);
end

#! Calculate color gradient
function _color_grad(sim::M2PhaseSim, cbounds::Matrix{Int64}, i::Int, j::Int)
  F = zeros(2);
  for k=1:sim.simr.lat.n-1
    const i_nbr   =   i + sim.simr.lat.c[1, k];
    const j_nbr   =   j + sim.simr.lat.c[2, k];

    if inbounds(i_nbr, j_nbr, cbounds)
      F += sub(sim.simr.lat.c, :, k) * (sim.simr.msm.rho[i_nbr, j_nbr] - 
                                        sim.simb.msm.rho[i_nbr, j_nbr]);
    end
  end
  return F;
end

#! Calculate color gradient
function _color_grad(sim::M2PhaseSim, active_cells::Matrix{Bool}, i::Int, j::Int)
  F = zeros(2);
  for k=1:sim.simr.lat.n-1
    const i_nbr   =   i + sim.simr.lat.c[1, k];
    const j_nbr   =   j + sim.simr.lat.c[2, k];

    if active_cells[i_nbr, j_nbr]
      F += sub(sim.simr.lat.c, :, k) * (sim.simr.msm.rho[i_nbr, j_nbr] - 
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
      const F         =   _color_grad(sim, cbounds, i, j);
      const Fmag      =   norm(F, 2);
      if (Fmag*Fmag > eps())
        for k=1:sim.simr.lat.n
          Fdotcsq   =   (dot(F, sub(sim.simr.lat.c, :, k)))^2;
          sim.simr.lat.f[k, i, j] += if Fmag != 0.0
                                       (sim.Ar/2 * Fmag * 
                                         (sim.simr.lat.w[k] * 
                                          Fdotcsq/(Fmag*Fmag) - __B[k]));
                                     else
                                       0.0;
                                     end
          if isnan(sim.simr.lat.f[k, i, j])
            @show F, sub(sim.simr.lat.c, :, k)
            @show Fdotcsq, sim.Ar/2, Fmag, sim.simr.lat.w[k], __B[k];
            @show k, i, j, sim.simr.lat.f[k, i, j]
          end
          sim.simb.lat.f[k, i, j] += if Fmag != 0.0
                                       (sim.Ab/2 * Fmag * 
                                        (sim.simb.lat.w[k] * 
                                         Fdotcsq/(Fmag*Fmag) - __B[k]));
                                     else
                                       0.0;
                                     end
        end
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
      const Fmag      =   norm(F, 2);
      if (Fmag*Fmag > eps())
        for k=1:sim.simr.lat.n
          Fdotcsq   =   (dot(F, sub(sim.simr.lat.c, :, k)))^2;
          @show Fdotcsq, sim.Ar/2, Fmag, sim.simr.lat.w[k], __B[k];
          sim.simr.lat.f[k, i, j] += if Fmag != 0.0
                                       (sim.Ar/2 * Fmag * 
                                         (sim.simr.lat.w[k] * 
                                          Fdotcsq/(Fmag*Fmag) - __B[k]));
                                     else
                                       0.0;
                                     end
          @show sim.simr.lat.f[k, i, j]
          sim.simb.lat.f[k, i, j] += if Fmag != 0.0
                                       (sim.Ab/2 * Fmag * 
                                        (sim.simb.lat.w[k] * 
                                         Fdotcsq/(Fmag*Fmag) - __B[k]));
                                     else
                                       0.0;
                                     end
        end
      end
    end
  end
end
