const __collision_root__ = dirname(@__FILE__);
require(abspath(joinpath(__collision_root__, "lattice.jl")));
require(abspath(joinpath(__collision_root__, "multiscale.jl")));

#! Single relaxation time collision function for incompressible Newtonian flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
function srt_col_f! (lat::Lattice, msm::MultiscaleMap)

  const c_ssq = @c_ssq(lat.dx, lat.dt);
  const ni, nj = size(lat.f);

  for i=1:ni, j=1:nj, k=1:9
    f_eq = incomp_f_eq(msm.rho[i,j], lat.w[k], c_ssq, vec(lat.c[k,:]),
      vec(msm.u[i,j,:]));
    lat.f[i,j,k] = msm.omega * f_eq + (1.0 - msm.omega) * lat.f[i,j,k];
  end

end

#! Equilibrium frequency distribution for incompressible Newtonian flow
#!
#! \param rho Density at lattice site
#! \param w Weight for lattice direction
#! \param c_ssq Lattice speed of sound squared
#! \param c_k Vector for lattice direction
#! \param u Macroscopic flow at lattice site
#! \return Equilibrium frequency
function incomp_f_eq(rho::Float64, w::Float64, c_ssq::Float64,
  c_k::Array{Float64, 1}, u::Array{Float64, 1})

  const ckdotu = dot(c_k, u);

  return rho * w * (1.0 + ckdotu/(c_ssq) + 0.5*(ckdotu*ckdotu)/(c_ssq*c_ssq)
    - 0.5 * dot(u, u) / (c_ssq));
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
#! \param M Transmation matrix to map f from velocity space to momentum space
#! \param S (Sparse) diagonal relaxation matrix
function mrt_col_f! (lat::Lattice, msm::MultiscaleMap, M::Array{Float64,2},
  S::SparseMatrixCSC)

  c_ssq = @c_ssq(lat.dx, lat.dt);
  ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  for i=1:ni, j=1:nj
    for k=1:9
      f_eq[k] = incomp_f_eq(msm.rho[i,j], lat.w[k], c_ssq, vec(lat.c[k,:]),
        vec(msm.u[i,j,:]));
    end

    f = lat.f[i,j,:];
    m = M * f;
    m_eq = M * f_eq;

    lat.f[i,j,:] = f - inv(M) * S * (m - m_eq); # perform collision
  end
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
#! \param M Transmation matrix to map f from velocity space to momentum space
#! \param S Function for creating a sparse diagonal relaxation matrix
function mrt_col_f! (lat::Lattice, msm::MultiscaleMap, M::Array{Float64,2},
  S::Function)

  c_ssq = @c_ssq(lat.dx, lat.dt);
  ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  for i=1:ni, j=1:nj
    for k=1:9
      f_eq[k] = incomp_f_eq(msm.rho[i,j], lat.w[k], c_ssq, vec(lat.c[k,:]),
        vec(msm.u[i,j,:]));
    end

    f = lat.f[i,j,:];
    m = M * f;
    m_eq = M * f_eq;

    lat.f[i,j,:] = f - inv(M) * S(lat, msm) * (m - m_eq); # perform collision
  end
end

#! Initializes the default multiple relaxation time transformation matrix
macro DEFAULT_MRT_M()
  return :( [1.0    1.0    1.0    1.0    1.0    1.0    1.0    1.0    1.0;
            -4.0   -1.0   -1.0   -1.0   -1.0    2.0    2.0    2.0    2.0;
             4.0   -2.0   -2.0   -2.0   -2.0    1.0    1.0    1.0    1.0;
             0.0    1.0    0.0   -1.0    0.0    1.0   -1.0   -1.0    1.0;
             0.0   -2.0    0.0    2.0    0.0    1.0   -1.0   -1.0    1.0;
             0.0    0.0    1.0    0.0   -1.0    1.0    1.0   -1.0   -1.0;
             0.0    0.0   -2.0    0.0    2.0    1.0    1.0   -1.0   -1.0;
             0.0    1.0   -1.0    1.0   -1.0    0.0    0.0    0.0    0.0;
             0.0    0.0    0.0    0.0    0.0    1.0   -1.0    1.0   -1.0]
          );
end

#! Multiple relaxation time collision function for incompressible flow
#!
#! \param lat Lattice
#! \param msm Multiscale map
#! \param S (Sparse) diagonal relaxation matrix
function mrt_col_f! (lat::Lattice, msm::MultiscaleMap, S::SparseMatrixCSC)
  const M = @DEFAULT_MRT_M();
  const c_ssq = @c_ssq(lat.dx, lat.dt);
  ni, nj = size(lat.f);

  # calc f_eq vector ((f_eq_1, f_eq_2, ..., f_eq_9))
  f_eq = Array(Float64, 9);
  for i=1:ni, j=1:nj
    for k=1:9
      f_eq[k] = incomp_f_eq(msm.rho[i,j], lat.w[k], c_ssq, vec(lat.c[k,:]),
        vec(msm.u[i,j,:]));
    end

    f = lat.f[i,j,:];
    m = M * f;
    m_eq = M * f_eq;

    lat.f[i,j,:] = f - inv(M) * S * (m - m_eq); # perform collision
  end
end
