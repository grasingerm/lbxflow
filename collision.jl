require("lattice.jl");
require("multiscale.jl");

function srt_col_f! (lat::Lattice, msm::MultiscaleMap)
  
  c_s = lat.dx / (sqrt(3) * lat.dt);
  ni, nj = size(lat.f);

  for i=1:ni
    for j=1:nj
      for k=1:9
        f_eq = incomp_f_eq(msm.rho[i,j], lat.w[k], c_s, vec(lat.c[k,:]), 
          vec(msm.u[i,j,:]));
        lat.f[i,j,k] -= 1.0/msm.tau * (lat.f[i,j,k] - f_eq);
      end
    end
  end

end

function incomp_f_eq(rho::Float64, w::Float64, c_s::Float64, 
  c_k::Array{Float64, 1}, u::Array{Float64, 1})
  
  ckdotu = dot(c_k, u);

  return (1.0 + ckdotu/(c_s*c_s) + 0.5*(ckdotu*ckdotu)/(c_s*c_s*c_s*c_s) - 0.5 *
    dot(u, u) / (c_s*c_s));
end
