__collision_root__ = dirname(@__FILE__);
require(abspath(joinpath(__collision_root__, "lattice.jl")));
require(abspath(joinpath(__collision_root__, "multiscale.jl")));

function srt_col_f! (lat::Lattice, msm::MultiscaleMap)

  c_ssq = (lat.dx * lat.dx) / (3 * lat.dt * lat.dt);
  ni, nj = size(lat.f);

  for i=1:ni, j=1:nj, k=1:9
    f_eq = incomp_f_eq(msm.rho[i,j], lat.w[k], c_ssq, vec(lat.c[k,:]),
      vec(msm.u[i,j,:]));
    lat.f[i,j,k] = msm.omega * f_eq + (1.0 - msm.omega) * lat.f[i,j,k];
  end

end

function incomp_f_eq(rho::Float64, w::Float64, c_ssq::Float64,
  c_k::Array{Float64, 1}, u::Array{Float64, 1})

  ckdotu = dot(c_k, u);

  return rho * w * (1.0 + ckdotu/(c_ssq) + 0.5*(ckdotu*ckdotu)/(c_ssq*c_ssq)
    - 0.5 * dot(u, u) / (c_ssq));
end
