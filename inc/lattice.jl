abstract Lattice

#! D2Q9 lattice
type LatticeD2Q9 <: Lattice
  dx::FloatingPoint
  dt::FloatingPoint
  cc::FloatingPoint
  cs::FloatingPoint
  cssq::FloatingPoint
  f::Array{Float64, 3}
  c::Matrix{Float64}
  w::Vector{Float64}
  n::Int

  # Default lattice speed vectors and associated weights
  const cdef = [1.0 0.0; 0.0 1.0; -1.0 0.0; 0.0 -1.0; 1.0 1.0; -1.0 1.0;
                -1.0 -1.0; 1.0 -1.0; 0.0 0.0]';
  const wdef = [1.0/9.0; 1.0/9.0; 1.0/9.0; 1.0/9.0; 1.0/36.0; 1.0/36.0;
                1.0/36.0; 1.0/36.0; 4.0/9.0];
  const n = length(wdef);

  # Lattice functions for computing `c`, `cs`, and `cssq`
  cf(dx, dt)    = dx/dt;
  csf(dx, dt)   = cf(dx, dt) / sqrt(3);
  cssqf(dx, dt) = csf(dx, dt) * csf(dx, dt);

  LatticeD2Q9(dx::FloatingPoint, dt::FloatingPoint, ni::Int, nj::Int) =
    new(dx, dt, cf(dx, dt), csf(dx, dt), cssqf(dx, dt),
        zeros(Float64, (n, ni, nj)), cdef, wdef);

  LatticeD2Q9(dx::FloatingPoint, dt::FloatingPoint, f::Array{Float64, 3}) =
    new(dx, dt, cf(dx, dt), csf(dx, dt), cssqf(dx, dt), f, cdef, wdef);

  function LatticeD2Q9(dx::FloatingPoint, dt::FloatingPoint, ni::Int, nj::Int,
                   rho::FloatingPoint)
    f = zeros(Float64, (n, ni, nj));
    for k=1:n
      f[k,:,:] = fill(rho * wdef[k], (ni, nj));
    end
    new(dx, dt, cf(dx, dt), csf(dx, dt), cssqf(dx, dt), f, cdef, wdef);
  end
end

#! D2Q4 lattice
type LatticeD2Q4 <: Lattice
  dx::FloatingPoint
  dt::FloatingPoint
  cc::FloatingPoint
  cs::FloatingPoint
  cssq::FloatingPoint
  f::Array{Float64, 3}
  c::Matrix{Float64}
  w::Vector{Float64} 
  n::Int

  # Default lattice speed vectors and associated weights
  const cdef = [1.0 0.0; -1.0 0.0; 0.0 1.0; 0.0 -1.0]';
  const n = size(cdef, 2);
  const wdef = fill(1.0/4.0, n);

  # Lattice functions for computing `c`, `cs`, and `cssq`
  cf(dx, dt)    = dx/dt;
  csf(dx, dt)   = cf(dx, dt) / sqrt(2);
  cssqf(dx, dt) = csf(dx, dt) * csf(dx, dt);

  LatticeD2Q4(dx::FloatingPoint, dt::FloatingPoint, ni::Int, nj::Int) =
    new(dx, dt, cf(dx, dt), csf(dx, dt), cssqf(dx, dt),
        zeros(Float64, (n, ni, nj)), cdef, wdef);

  LatticeD2Q4(dx::FloatingPoint, dt::FloatingPoint, f::Array{Float64, 3}) =
    new(dx, dt, cf(dx, dt), csf(dx, dt), cssqf(dx, dt), f, cdef, wdef);

  function LatticeD2Q4(dx::FloatingPoint, dt::FloatingPoint, ni::Int, nj::Int,
                   rho::Float64)
    f = zeros(Float64, (n, ni, nj));
    for k=1:n
      f[k,:,:] = fill(rho * wdef[k], (ni, nj));
    end
    new(dx, dt, cf(dx, dt), csf(dx, dt), cssqf(dx, dt), f, cdef, wdef);
  end
end
