const __boundary_root__ = dirname(@__FILE__);
require(abspath(joinpath(__boundary_root__, "lattice.jl")));

#! Bounceback boundary condition for north boundary
function north_bounce_back!(lat::Lattice, i_begin::Int, i_end::Int, j::Int)
  for i=i_begin:i_end
    lat.f[i,j,5] = lat.f[i,j,3];
    lat.f[i,j,8] = lat.f[i,j,6];
    lat.f[i,j,9] = lat.f[i,j,7];
  end
end

#! Bounceback boundary condition for north boundary of domain
function north_bounce_back!(lat::Lattice)
  const ni, nj = size(lat.f);
  north_bounce_back!(lat, 1, ni, nj);
end

#! Bounceback boundary condition for south boundary
function south_bounce_back!(lat::Lattice, i_begin::Int, i_end::Int, j::Int)
  for i=i_begin:i_end
    lat.f[i,j,3] = lat.f[i,j,5];
    lat.f[i,j,6] = lat.f[i,j,8];
    lat.f[i,j,7] = lat.f[i,j,9];
  end
end

#! Bounceback boundary condition for south boundary of domain
function south_bounce_back!(lat::Lattice)
  south_bounce_back!(lat, 1, size(lat.f)[1], 1);
end

#! Bounceback boundary condition for west boundary
function east_bounce_back!(lat::Lattice, i::Int, j_begin::Int, j_end::Int)
  for j=j_begin:j_end
    lat.f[i,j,4] = lat.f[i,j,2];
    lat.f[i,j,8] = lat.f[i,j,6];
    lat.f[i,j,7] = lat.f[i,j,9];
  end
end

#! Bounceback boundary condition for west boundary of domain
function east_bounce_back!(lat::Lattice)
  const ni, nj = size(lat.f);
  west_bounce_back!(lat, ni, 1, nj);
end

#! Bounceback boundary condition for west boundary
function west_bounce_back!(lat::Lattice, i::Int, j_begin::Int, j_end::Int)
  for j=j_begin:j_end
    lat.f[i,j,2] = lat.f[i,j,4];
    lat.f[i,j,6] = lat.f[i,j,8];
    lat.f[i,j,9] = lat.f[i,j,7];
  end
end

#! Bounceback boundary condition for west boundary of domain
function west_bounce_back!(lat::Lattice)
  const nj = size(lat.f)[2];
  west_bounce_back!(lat, 1, 1, nj);
end

#! West inlet boundary condition
function west_inlet!(lat::Lattice, u::FloatingPoint, i::Int, j_begin::Int,
  j_end::Int)

  for j=j_begin:j_end
    rhow = (lat.f[i,j,1] + lat.f[i,j,3] + lat.f[i,j,5] +
            2.0 * (lat.f[i,j,4] + lat.f[i,j,7] + lat.f[i,j,8])) / (1.0 - u);
    lat.f[i,j,2] = lat.f[i,j,4] + 2.0 * rhow * u / 3.0;
    lat.f[i,j,6] = lat.f[i,j,8] + rhow * u / 6.0;
    lat.f[i,j,9] = lat.f[i,j,7] + rhow * u / 6.0;
  end
end

#! West inlet boundary condition
function west_inlet!(lat::Lattice, u::FloatingPoint)
  west_inlet!(lat, u, 1, 1, size(lat.f)[2]);
end

#! East open boundary
function east_open!(lat::Lattice, i::Int, j_begin::Int, j_end::Int)

  for j=j_begin:j_end
    lat.f[i,j,2] = 2.0 * lat.f[i-1,j,2] - lat.f[i-2,j,2];
    lat.f[i,j,6] = 2.0 * lat.f[i-1,j,6] - lat.f[i-2,j,6];
    lat.f[i,j,9] = 2.0 * lat.f[i-1,j,9] - lat.f[i-2,j,9];
  end
end

#! East open boundary
function east_open!(lat::Lattice)
  const ni, nj = size(lat.f);
  east_open!(lat, ni, 1, nj);
end

#! Periodic boundary condition
function periodic!(lat::Lattice, is::Array{Int, 1}, ks::Array{Int, 1},
  j_from::Int, j_to::Int)

  const ni, nj = size(lat.f);

  for i in is, k in ks
    c = lat.c[k,:];

    # where to map frequency distribution to
    i_tok = i + c[1];
    j_tok = j_to + c[2];

    # correct for boundary edges
    if i_tok > ni; i_tok = ni; end;
    if i_tok < 1; i_tok = 1; end;
    if j_tok > nj; j_tok = nj; end;
    if j_tok < 1; j_tok = 1; end;

    # stream
    lat.f[i_tok,j_tok,k] = lat.f[i,j_from,k];
  end
end

#! Periodic boundary condition
function periodic!(lat::Lattice, i_from::Int, i_to::Int, js::Array{Int, 1},
  ks::Array{Int, 1})

  const ni, nj = size(lat.f);

  for j in js, k in ks
    c = lat.c[k,:];

    # where to map frequency distribution to
    i_tok = i_to + c[1];
    j_tok = j + c[2];

    # correct for boundary edges
    if i_tok > ni; i_tok = ni; end;
    if i_tok < 1; i_tok = 1; end;
    if j_tok > nj; j_tok = nj; end;
    if j_tok < 1; j_tok = 1; end;

    # stream
    lat.f[i_tok,j_tok,k] = lat.f[i_from,j,k];
  end
end

#! Periodic east to west
function periodic_east_to_west!(lat::Lattice)
  const ni, nj = size(lat.f);

  for j=1:nj, k in (2, 6, 9)
    lat.f[1,j,k] = lat.f[ni-1,j,k];
  end

  for j=1:nj, k in (4, 7, 8)
    lat.f[ni,j,k] = lat.f[2,j,k];
  end

end

#! Periodic east to west
function periodic_north_to_south!(lat::Lattice)
  const ni, nj = size(lat.f);

  for i=1:ni, k in (3, 6, 7)
    lat.f[i,1,k] = lat.f[i,nj-1,k];
  end

  for i=1:ni, k in (5, 8, 9)
    lat.f[i,nj,k] = lat.f[i,2,k];
  end

end

#! Lid driven flow
function lid_driven!(lat::Lattice, u::FloatingPoint)
  const ni, nj = size(lat.f)

  for i=1:ni
    rho = lat.f[i,nj,1] + lat.f[i,nj,2] + lat.f[i,nj,4] + 2.0 * (lat.f[i,nj,3]
      + lat.f[i,nj,7] + lat.f[i,nj,6]);
    lat.f[i,nj,5] = lat.f[i,nj,3];
    lat.f[i,nj,9] = lat.f[i,nj,7] + rho * u / 6.0;
    lat.f[i,nj,8] = lat.f[i,nj,6] - rho * u / 6.0;
  end
end

#! Pressure west direction
function west_pressure!(lat::Lattice, rho_in::FloatingPoint, i::Int,
  j_begin::Int, j_end::Int)

  for j=j_begin:j_end
    u_x = 1 - (lat.f[i,j,1] + lat.f[i,j,3] + lat.f[i,j,5] +
                2 * (lat.f[i,j,4] + lat.f[i,j,7] + lat.f[i,j,8])) / rho_in;
    lat.f[i,j,2] = lat.f[i,j,4] + 2/3 * rho_in * u_x;

    second_term = 0.5 * (lat.f[i,j,3] - lat.f[i,j,5]);
    third_term = 1/6 * rho_in * u_x;

    lat.f[i,j,6] = lat.f[i,j,8] - second_term + third_term;
    lat.f[i,j,9] = lat.f[i,j,7] + second_term + third_term;
  end
end

#! Pressure west direction
function west_pressure!(lat::Lattice, rho_in::FloatingPoint)
  const ni, nj = size(lat.f);

  # middle
  west_pressure!(lat, rho_in, 1, 2, nj-1);

  # bottom corner
  lat.f[1,1,2] = lat.f[1,1,4];
  lat.f[1,1,3] = lat.f[1,1,5];
  lat.f[1,1,6] = lat.f[1,1,8];
  f68 = 0.5 * (rho_in - (lat.f[1,1,1] + lat.f[1,1,2] + lat.f[1,1,3]
          + lat.f[1,1,4] + lat.f[1,1,5] + lat.f[1,1,6] + lat.f[1,1,8]));
  lat.f[1,1,7] = f68;
  lat.f[1,1,9] = f68;

  # top corner
  lat.f[1,nj,2] = lat.f[1,nj,5];
  lat.f[1,nj,5] = lat.f[1,nj,3];
  lat.f[1,nj,7] = lat.f[1,nj,9];
  f57 = 0.5 * (rho_in - (lat.f[1,nj,1] + lat.f[1,nj,2] + lat.f[1,nj,3]
          + lat.f[1,nj,4] + lat.f[1,nj,5] + lat.f[1,nj,7] + lat.f[1,nj,9]));
  lat.f[1,nj,6] = f57;
  lat.f[1,nj,8] = f57;
end

#! Pressure east direction
function east_pressure!(lat::Lattice, rho_out::FloatingPoint, i::Int,
  j_begin::Int, j_end::Int)

  for j=j_begin:j_end
    u_x = (lat.f[i,j,1] + lat.f[i,j,3] + lat.f[i,j,5] +
                2 * (lat.f[i,j,2] + lat.f[i,j,6] + lat.f[i,j,9])) / rho_out - 1;
    lat.f[i,j,4] = lat.f[i,j,2] - 2/3 * rho_out * u_x;

    second_term = 0.5 * (lat.f[i,j,3] - lat.f[i,j,5]);
    third_term = 1/6 * rho_out * u_x;

    lat.f[i,j,8] = lat.f[i,j,6] + second_term - third_term;
    lat.f[i,j,7] = lat.f[i,j,9] - second_term - third_term;
  end
end

#! Pressure east direction
function east_pressure!(lat::Lattice, rho_out::FloatingPoint)
  const ni, nj = size(lat.f);

  # middle of inlet
  east_pressure!(lat, rho_out, ni, 2, nj-1);

  # bottom corner
  lat.f[ni,1,4] = lat.f[ni,1,2];
  lat.f[ni,1,3] = lat.f[ni,1,5];
  lat.f[ni,1,8] = lat.f[ni,1,6];
  f68 = 0.5 * (rho_out - (lat.f[ni,1,1] + lat.f[ni,1,2] + lat.f[ni,1,3]
          + lat.f[ni,1,4] + lat.f[ni,1,5] + lat.f[ni,1,6] + lat.f[ni,1,8]));
  lat.f[ni,1,7] = f68;
  lat.f[ni,1,9] = f68;

  # top corner
  lat.f[ni,nj,4] = lat.f[ni,nj,2];
  lat.f[ni,nj,5] = lat.f[ni,nj,3];
  lat.f[ni,nj,7] = lat.f[ni,nj,9];
  f57 = 0.5 * (rho_out - (lat.f[ni,nj,1] + lat.f[ni,nj,2] + lat.f[ni,nj,3]
          + lat.f[ni,nj,4] + lat.f[ni,nj,5] + lat.f[ni,nj,7] + lat.f[ni,nj,9]));
  lat.f[ni,nj,6] = f57;
  lat.f[ni,nj,8] = f57;
end

#! Zou and He pressure boundary on north side
function zou_pressure_north!(lat::Lattice, rhoo::FloatingPoint)
  const ni, nj = size(lat.f);

  for i=1:ni
    v = -1. + (lat.f[i,nj,1] + lat.f[i,nj,2] + lat.f[i,nj,4]
        + 2. * (lat.f[i,nj,3] + lat.f[i,nj,6] + lat.f[i,nj,7])) / rhoo;
    ru = rhoo * v;
    lat.f[i,nj,5] = lat.f[i,nj,3] - (2./3.)*ru;
    lat.f[i,nj,8] = lat.f[i,nj,6] - (1./6.)*ru
                      + 0.5 * (lat.f[i,nj,2] - lat.f[i,nj,4]);
    lat.f[i,nj,9] = lat.f[i,nj,7] - (1./6.)*ru
                      + 0.5 * (lat.f[i,nj,4] - lat.f[i,nj,2]);
  end
end

#! Zou and He pressure boundary on south side
function zou_pressure_south!(lat::Lattice, rhoo::FloatingPoint)
  const ni, nj = size(lat.f);

  for i=1:ni
    v = -1. + (lat.f[i,1,1] + lat.f[i,1,2] + lat.f[i,1,4]
        + 2. * (lat.f[i,1,5] + lat.f[i,1,8] + lat.f[i,1,9])) / rhoo;
    ru = rhoo * v;
    lat.f[i,1,3] = lat.f[i,1,5] - (2./3.)*ru;
    lat.f[i,1,6] = lat.f[i,1,8] - (1./6.)*ru
                      + 0.5 * (lat.f[i,1,4] - lat.f[i,1,2]);
    lat.f[i,1,7] = lat.f[i,1,9] - (1./6.)*ru
                      + 0.5 * (lat.f[i,1,2] - lat.f[i,1,4]);
  end
end
