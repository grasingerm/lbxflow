type Lattice
  dx::FloatingPoint
  dt::FloatingPoint
  f::Array{Float64, 3}

  c = [0.0 0.0; 1.0 0.0; 0.0 1.0; -1.0 0.0; 0.0 -1.0; 1.0 1.0; -1.0 1.0;
        -1.0 -1.0; 1.0 -1.0];
  w = [4.0/9.0; 1.0/9.0; 1.0/9.0; 1.0/9.0; 1.0/9.0; 1.0/36.0; 1.0/36.0; 
        1.0/36.0; 1.0/36.0];

  Lattice(dx::FloatingPoint, dt::FloatingPoint, ni::Int, nj::Int) =
    new(dx, dt, zeros(Float64, (ni, nj, 9)));

  function Lattice(dx::FloatingPoint, dt::FloatingPoint, ni::Int, nj::Int, 
    rho::Float64)

    f = zeros(Float64, (ni, nj, 9));
    for k=1:9
      f[:,:,k] = fill(rho*w[k], (ni, nj));
    end

    new(dx, dt, f);
  end

end

#! Lattice speed of sound squared
macro c_ssq(lat::Lattice)
  return :($lat.dx * $lat.dx) / (3 * $lat.dt * $lat.dt);
end
