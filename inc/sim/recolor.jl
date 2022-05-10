# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Recoloring step, reduces color diffusion at interface
function recolor!(sim::M2PhaseSim, cbounds::Matrix{Int64}, feq_fa::LBXFunction,
                  feq_fb::LBXFunction)
  nbounds = size(cbounds, 2);
  u0      = Float64[0; 0];

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = sub(cbounds, :, r);
    for j=j_min:j_max, i=i_min:i_max
      # TODO should we have just cached the color gradient somewhere by now?
      F       =   _color_grad(sim, cbounds, i, j);
      R       =   sim.simr.msm.rho[i, j];
      B       =   sim.simb.msm.rho[i, j];
      RB      =   R * B;
      Σ       =   R + B;
      Σ2      =   Σ * Σ;
      if Σ2 > eps() && Σ > eps()
        for k=1:sim.simr.lat.n-1
          ck      =   sub(sim.simr.lat.c, :, k);
          Fmag    =   norm(F, 2);
          if Fmag > eps()
            cosϕ    =   dot(F, ck) / (norm(F, 2) * norm(ck, 2));
            N_k     =   sim.simr.lat.f[k, i, j] + sim.simb.lat.f[k, i, j];
            @show R, B, RB, Σ, Σ2, N_k, cosϕ
            sim.simr.lat.f[k, i, j] = (R / Σ * N_k + sim.β * RB/Σ2 * 
                                       feq_fa(sim.simr.lat, R, u0, k) * cosϕ);
            @show sim.simr.lat.f[k, i, j];
            sim.simb.lat.f[k, i, j] = (B / Σ * N_k - sim.β * RB/Σ2 * 
                                       feq_fb(sim.simb.lat, B, u0, k) * cosϕ);
          end
        end
      end
    end
  end
end

#! Recoloring step, reduces color diffusion at interface
function recolor!(sim::M2PhaseSim, active_cells::Matrix{Bool}, 
                  feq_fa::LBXFunction, feq_fb::LBXFunction)
  u0      = Float64[0; 0];

  for j=1:nj, i=1:ni
    if active_cells[i, j]
      # TODO should we have just cached the color gradient somewhere by now?
      F       =   _color_grad(sim, cbounds, i, j);
      R       =   sim.simr.msm.rho[i, j];
      B       =   sim.simb.msm.rho[i, j];
      RB      =   R * B;
      Σ       =   R + B;
      Σ2      =   Σ * Σ;
      if Σ2 > eps() && Σ > eps()
        for k=1:sim.simr.lat.n
          ck      =   sub(sim.simr.lat.c, :, k);
          cosϕ    =   dot(F, ck);
          N_k     =   sim.simr.lat.f[k, i, j] + sim.simb.lat.f[k, i, j];
          sim.simr.lat.f[k, i, j] = (R / Σ * N_k + sim.β * RB/Σ2 * 
                                     feq_fa(sim.simr.lat, R, u0, k) * cosϕ);
          sim.simb.lat.f[k, i, j] = (B / Σ * N_k - sim.β * RB/Σ2 * 
                                     feq_fb(sim.simb.lat, B, u0, k) * cosϕ);
        end
      end
    end
  end
end
