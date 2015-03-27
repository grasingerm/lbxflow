type Lattice
  dx::FloatingPoint
  dt::FloatingPoint
  f::Array{Float64, 3}
  c::Array{Float64, 2}
  w::Array{Float64, 1}

  cdef = [0.0 0.0; 1.0 0.0; 0.0 1.0; -1.0 0.0; 0.0 -1.0; 1.0 1.0; -1.0 1.0;
        -1.0 -1.0; 1.0 -1.0];
  wdef = [4.0/9.0; 1.0/9.0; 1.0/9.0; 1.0/9.0; 1.0/9.0; 1.0/36.0; 1.0/36.0;
        1.0/36.0; 1.0/36.0];

  Lattice(dx::FloatingPoint, dt::FloatingPoint, ni::Int, nj::Int) =
    new(dx, dt, zeros(Float64, (ni, nj, 9)), cdef, wdef);

  function Lattice(dx::FloatingPoint, dt::FloatingPoint, ni::Int, nj::Int,
    rho::Float64)

    f = zeros(Float64, (ni, nj, 9));
    for k=1:9
      f[:,:,k] = fill(rho * wdef[k], (ni, nj));
    end

    new(dx, dt, f, cdef, wdef);
  end

end

#! Lattice speed of sound squared
macro c_ssq(dx::Expr, dt::Expr)
  return :(($dx * $dx) / (3 * $dt * $dt));
end
