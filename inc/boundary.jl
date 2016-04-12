# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

macro wrap_lat_bc(sim, bc)
  return :(@anon sim -> $bc(sim.lat));
end

# Simulation to lattice wrappers
north_bounce_back!(sim::AbstractSim)       =  north_bounce_back!(sim.lat);
south_bounce_back!(sim::AbstractSim)       =  south_bounce_back!(sim.lat);
east_bounce_back!(sim::AbstractSim)        =  east_bounce_back!(sim.lat);
west_bounce_back!(sim::AbstractSim)        =  west_bounce_back!(sim.lat);

north_reflect!(sim::AbstractSim)           =  north_reflect!(sim.lat);
south_reflect!(sim::AbstractSim)           =  south_reflect!(sim.lat);
east_reflect!(sim::AbstractSim)            =  east_reflect!(sim.lat);
west_reflect!(sim::AbstractSim)            =  west_reflect!(sim.lat);

periodic_east_to_west!(sim::AbstractSim)   =  periodic_east_to_west!(sim.lat); 
periodic_north_to_south!(sim::AbstractSim) =  periodic_north_to_south!(sim.lat); 

north_velocity!(sim::AbstractSim, u::Real) =  north_velocity!(sim.lat, u);
south_velocity!(sim::AbstractSim, u::Real) =  south_velocity!(sim.lat, u);
east_velocity!(sim::AbstractSim, u::Real)  =  east_velocity!(sim.lat, u);
west_velocity!(sim::AbstractSim, u::Real)  =  west_velocity!(sim.lat, u);

north_pressure!(sim::AbstractSim, ρ::Real) =  north_pressure!(sim.lat, ρ);
south_pressure!(sim::AbstractSim, ρ::Real) =  south_pressure!(sim.lat, ρ);
east_pressure!(sim::AbstractSim, ρ::Real)  =  east_pressure!(sim.lat, ρ);
west_pressure!(sim::AbstractSim, ρ::Real)  =  west_pressure!(sim.lat, ρ);

north_open!(sim::AbstractSim)              =  north_open!(sim.lat);
south_open!(sim::AbstractSim)              =  south_open!(sim.lat);
east_open!(sim::AbstractSim)               =  east_open!(sim.lat);
west_open!(sim::AbstractSim)               =  west_open!(sim.lat);

#! Bounceback boundary condition for north boundary
function north_bounce_back!(lat::LatticeD2Q9, i_begin::Int, i_end::Int, j::Int)
  for i=i_begin:i_end
    lat.f[4, i, j] = lat.f[2, i, j];
    lat.f[7, i, j] = lat.f[5, i, j];
    lat.f[8, i, j] = lat.f[6, i, j];
  end
end

#! Bounceback boundary condition for north boundary of domain
function north_bounce_back!(lat::Lattice)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  north_bounce_back!(lat, 1, ni, nj);
end

#! Bounceback boundary condition for south boundary
function south_bounce_back!(lat::LatticeD2Q9, i_begin::Int, i_end::Int, j::Int)
  for i=i_begin:i_end
    lat.f[2, i, j] = lat.f[4, i, j];
    lat.f[5, i, j] = lat.f[7, i, j];
    lat.f[6, i, j] = lat.f[8, i, j];
  end
end

#! Bounceback boundary condition for south boundary of domain
function south_bounce_back!(lat::Lattice)
  south_bounce_back!(lat, 1, size(lat.f, 2), 1);
end

#! Bounceback boundary condition for east boundary
function east_bounce_back!(lat::LatticeD2Q9, i::Int, j_begin::Int, j_end::Int)
  for j=j_begin:j_end
    lat.f[3, i, j] = lat.f[1, i, j];
    lat.f[7, i, j] = lat.f[5, i, j];
    lat.f[6, i, j] = lat.f[8, i, j];
  end
end

#! Bounceback boundary condition for east boundary of domain
function east_bounce_back!(lat::Lattice)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  east_bounce_back!(lat, ni, 1, nj);
end

#! Bounceback boundary condition for west boundary
function west_bounce_back!(lat::LatticeD2Q9, i::Int, j_begin::Int, j_end::Int)
  for j=j_begin:j_end
    lat.f[1, i, j] = lat.f[3, i, j];
    lat.f[5, i, j] = lat.f[7, i, j];
    lat.f[8, i, j] = lat.f[6, i, j];
  end
end

#! Bounceback boundary condition for west boundary of domain
function west_bounce_back!(lat::Lattice)
  const nj = size(lat.f, 3);
  west_bounce_back!(lat, 1, 1, nj);
end

#! Bounceback boundary condition for north boundary
function north_reflect!(lat::LatticeD2Q9, i_begin::Int, i_end::Int, j::Int)
  for i=i_begin:i_end
    lat.f[4, i,   j] = lat.f[2, i, j];
    lat.f[8, i+1, j] = lat.f[5, i, j];
    lat.f[7, i-1, j] = lat.f[6, i, j];
  end
end

#! Bounceback boundary condition for north boundary of domain
function north_reflect!(lat::Lattice)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  north_reflect!(lat, 2, ni-1, nj);
end

#! Bounceback boundary condition for south boundary
function south_reflect!(lat::LatticeD2Q9, i_begin::Int, i_end::Int, j::Int)
  for i=i_begin:i_end
    lat.f[2, i,   j] = lat.f[4, i, j];
    lat.f[5, i+1, j] = lat.f[8, i, j];
    lat.f[6, i-1, j] = lat.f[7, i, j];
  end
end

#! Bounceback boundary condition for south boundary of domain
function south_reflect!(lat::Lattice)
  south_reflect!!(lat, 2, size(lat.f, 2)-1, 1);
end

#! Bounceback boundary condition for east boundary
function east_reflect!(lat::LatticeD2Q9, i::Int, j_begin::Int, j_end::Int)
  for j=j_begin:j_end
    lat.f[3, i, j  ] = lat.f[1, i, j];
    lat.f[7, i, j-1] = lat.f[8, i, j];
    lat.f[6, i, j+1] = lat.f[5, i, j];
  end
end

#! Bounceback boundary condition for east boundary of domain
function east_reflect!(lat::Lattice)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  east_reflect!(lat, ni, 2, nj-1);
end

#! Bounceback boundary condition for west boundary
function west_reflect!(lat::LatticeD2Q9, i::Int, j_begin::Int, j_end::Int)
  for j=j_begin:j_end
    lat.f[1, i, j  ] = lat.f[3, i, j];
    lat.f[5, i, j+1] = lat.f[6, i, j];
    lat.f[8, i, j-1] = lat.f[7, i, j];
  end
end

#! Bounceback boundary condition for west boundary of domain
function west_reflect!(lat::Lattice)
  const nj = size(lat.f, 3);
  west_reflect!(lat, 1, 2, nj-1);
end

# TODO: refactor half bounce back boundary schemes
#! Bounceback boundary condition for north boundary
function north_half_bounce_back!(lat::LatticeD2Q9, i_begin::Int, i_end::Int, j::Int)
  const ks_to   = (4, 7, 8);
  const ks_from = (2, 5, 6);
  const cs      = (0, -1, 1);

  for i=i_begin:i_end, (k_to, k_from, c) in zip(ks_to, ks_from, cs)
    lat.f[k_to, i+c, j-1] = lat.f[k_from, i, j];
  end
end

#! Bounceback boundary condition for north boundary of domain
function north_half_bounce_back!(lat::LatticeD2Q9)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  north_half_bounce_back!(lat, 2, ni-1, nj);

  lat.f[4, 1,    nj-1]     =   lat.f[2, 1,  nj];
  lat.f[8, 2,    nj-1]     =   lat.f[6, 1,  nj];

  lat.f[4, ni,   nj-1]     =   lat.f[2, ni, nj];
  lat.f[7, ni-1, nj-1]     =   lat.f[5, ni, nj];
end

#! Bounceback boundary condition for south boundary
function south_half_bounce_back!(lat::LatticeD2Q9, i_begin::Int, i_end::Int, j::Int)
  const ks_from = (4, 7, 8);
  const ks_to   = (2, 5, 6);
  const cs      = (0, 1, -1);

  for i=i_begin:i_end, (k_to, k_from, c) in zip(ks_to, ks_from, cs)
    lat.f[k_to, i+c, j+1] = lat.f[k_from, i, j];
  end
end

#! Bounceback boundary condition for south boundary of domain
function south_half_bounce_back!(lat::LatticeD2Q9)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  south_half_bounce_back!(lat, 2, ni-1, 1);

  lat.f[2, 1, 2] = lat.f[4, 1, 1];
  lat.f[5, 2, 2] = lat.f[7, 1, 1];

  lat.f[2, ni,   2] = lat.f[4, ni, 1];
  lat.f[6, ni-1, 2] = lat.f[8, ni, 1];
end

#! Bounceback boundary condition for west boundary
function east_half_bounce_back!(lat::LatticeD2Q9, i::Int, j_begin::Int, j_end::Int)
  const ks_to   = (3, 6, 7);
  const ks_from = (1, 8, 5);
  const cs      = (0, 1, -1);

  for j=j_begin:j_end, (k_to, k_from, c) in zip(ks_to, ks_from, cs)
    lat.f[k_to, i-1, j+c] = lat.f[k_from, i, j];
  end
end

#! Bounceback boundary condition for west boundary of domain
function east_half_bounce_back!(lat::LatticeD2Q9)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  east_half_bounce_back!(lat, ni, 2, nj-1);

  lat.f[3, ni-1, 1] = lat.f[1, ni, 1];
  lat.f[6, ni-1, 2] = lat.f[8, ni, 1];

  lat.f[3, ni-1, nj  ] = lat.f[1, ni, nj];
  lat.f[7, ni-1, nj-1] = lat.f[5, ni, nj];
end

#! Bounceback boundary condition for west boundary
function west_half_bounce_back!(lat::LatticeD2Q9, i::Int, j_begin::Int, j_end::Int)
  const ks_to   = (1, 5, 8);
  const ks_from = (2, 7, 6);
  const cs      = (0, 1, -1);

  for j=j_begin:j_end, (k_to, k_from, c) in zip(ks_to, ks_from, cs)
    lat.f[k_to, i+1, j+c] = lat.f[k_from, i, j];
  end
end

#! Bounceback boundary condition for west boundary of domain
function west_half_bounce_back!(lat::LatticeD2Q9)
  const nj = size(lat.f, 2);
  west_half_bounce_back!(lat, 1, 2, nj-1);

  lat.f[1, 2, 1] = lat.f[3, 1, 1];
  lat.f[5, 2, 2] = lat.f[7, 1, 1];

  lat.f[1, 2, nj  ] = lat.f[3, 1, nj];
  lat.f[8, 2, nj-1] = lat.f[6, 1, nj];
end

#! Bounceback boundary condition for north boundary
function north_halfa_bounce_back!(lat::LatticeD2Q9, i_begin::Int, i_end::Int, j::Int)
  const ks_to   = (4, 7, 8);
  const ks_from = (2, 5, 6);

  for i=i_begin:i_end, (k_to, k_from) in zip(ks_to, ks_from)
    lat.f[k_to, i, j] = lat.f[k_from, i, j];
  end

  for i=i_begin:i_end
    (lat.f[1, i, j], lat.f[3, i, j]) = (lat.f[3, i, j], lat.f[1, i, j]);
  end
end

#! Bounceback boundary condition for north boundary of domain
function north_halfa_bounce_back!(lat::Lattice)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  north_halfa_bounce_back!(lat, 1, ni, nj);
end

#! Bounceback boundary condition for south boundary
function south_halfa_bounce_back!(lat::LatticeD2Q9, i_begin::Int, i_end::Int, 
                                  j::Int)
  const ks_from = (4, 7, 8);
  const ks_to   = (2, 5, 6);

  for i=i_begin:i_end, (k_to, k_from) in zip(ks_to, ks_from)
    lat.f[k_to, i, j] = lat.f[k_from, i, j];
  end

  for i=i_begin:i_end
    (lat.f[1, i, j], lat.f[3, i, j]) = (lat.f[3, i, j], lat.f[1, i, j]);
  end
end

#! Bounceback boundary condition for south boundary of domain
function south_halfa_bounce_back!(lat::Lattice)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  south_halfa_bounce_back!(lat, 1, ni, 1);
end

#! Bounceback boundary condition for west boundary
function east_halfa_bounce_back!(lat::LatticeD2Q9, i::Int, j_begin::Int, 
                                 j_end::Int)
  const ks_to   = (3, 6, 7);
  const ks_from = (3, 8, 5);

  for j=j_begin:j_end, (k_to, k_from) in zip(ks_to, ks_from)
    lat.f[k_to, i, j] = lat.f[k_from, i, j];
  end

  for i=i_begin:i_end
    (lat.f[2, i, j], lat.f[4, i, j]) = (lat.f[4, i, j], lat.f[2, i, j]);
  end
end

#! Bounceback boundary condition for west boundary of domain
function east_halfa_bounce_back!(lat::Lattice)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  east_halfa_bounce_back!(lat, ni, 1, nj);
end

#! Bounceback boundary condition for west boundary
function west_halfa_bounce_back!(lat::LatticeD2Q9, i::Int, j_begin::Int, j_end::Int)
  const ks_to   = (1, 5, 8);
  const ks_from = (2, 7, 6);

  for j=j_begin:j_end, (k_to, k_from) in zip(ks_to, ks_from)
    lat.f[k_to, i, j] = lat.f[k_from, i, j];
  end

  for i=i_begin:i_end
    (lat.f[2, i, j], lat.f[4, i, j]) = (lat.f[4, i, j], lat.f[2, i, j]);
  end
end

#! Bounceback boundary condition for west boundary of domain
function west_halfa_bounce_back!(lat::Lattice)
  const nj = size(lat.f, 3);
  west_halfa_bounce_back!(lat, 1, 1, nj);
end

# =========================================================================== #
# ============================ periodic BCs ================================= #
# =========================================================================== #

#! Periodic boundary condition
function periodic!(lat::LatticeD2Q9, is::Vector{Int}, ks::Vector{Int},
                   j_from::Int, j_to::Int)

  const ni, nj = size(lat.f, 2), size(lat.f, 3);

  for i in is, k in ks
    c = lat.c[k, :];

    # where to map frequency distribution to
    i_tok = i + c[1];
    j_tok = j_to + c[2];

    # correct for boundary edges
    if i_tok > ni;  i_tok = ni; end;
    if i_tok < 1;   i_tok = 1;  end;
    if j_tok > nj;  j_tok = nj; end;
    if j_tok < 1;   j_tok = 1;  end;

    # stream
    lat.f[k, i_tok, j_tok] = lat.f[k, i, j_from];
  end
end

#! Periodic boundary condition
function periodic!(lat::LatticeD2Q9, i_from::Int, i_to::Int, js::Vector{Int},
                   ks::Vector{Int})

  const ni, nj = size(lat.f, 2), size(lat.f, 3);

  for j in js, k in ks
    c = lat.c[k, :];

    # where to map frequency distribution to
    i_tok = i_to + c[1];
    j_tok = j + c[2];

    # correct for boundary edges
    if i_tok > ni;  i_tok = ni; end;
    if i_tok < 1;   i_tok = 1;  end;
    if j_tok > nj;  j_tok = nj; end;
    if j_tok < 1;   j_tok = 1;  end;

    # stream
    lat.f[k, i_tok, j_tok] = lat.f[k, i_from, j];
  end
end

#! Periodic east to west
function periodic_east_to_west!(lat::LatticeD2Q9)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);

  for j=1:nj, k in (1, 5, 8)
    lat.f[k, 1, j] = lat.f[k, ni-1, j];
  end

  for j=1:nj, k in (3, 6, 7)
    lat.f[k, ni, j] = lat.f[k, 2, j];
  end
end

#! Periodic north to south 
function periodic_north_to_south!(lat::LatticeD2Q9)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);

  for i=1:ni, k in (2, 5, 6)
    lat.f[k, i, 1] = lat.f[k, i, nj-1];
  end

  for i=1:ni, k in (4, 7, 8)
    lat.f[k, i, nj] = lat.f[k, i, 2];
  end
end

# =========================================================================== #
# ============================ velocity BCs ================================= #
# =========================================================================== #

#! North velocity boundary condition
function north_velocity!(lat::LatticeD2Q9, u::Real, i_begin::Int, i_end::Int,
                         j::Int)
  for i=i_begin:i_end
    const rhon = (lat.f[9, i, j] + lat.f[1, i, j] + lat.f[3, i, j] +
                  2.0 * (lat.f[2, i, j] + lat.f[5, i, j] + lat.f[6, i, j]));
    const jy   = rhon * u;    
    lat.f[4, i, j] = lat.f[2, i, j] - 2.0/3.0 * jy;
    lat.f[7, i, j] = (lat.f[5, i, j] - 1.0/6.0 * jy +
                      0.5 * (lat.f[1, i, j] - lat.f[3, i, j]));
    lat.f[8, i, j] = (lat.f[6, i, j] - 1.0/6.0 * jy +
                      0.5 * (lat.f[3, i, j] - lat.f[1, i, j]));
  end
end

#! North velocity boundary condition
function north_velocity!(lat::Lattice, u::Real)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  north_velocity!(lat, u, 1, ni, nj);
end

#! South velocity boundary condition
function south_velocity!(lat::LatticeD2Q9, u::Real, i_begin::Int,
                         i_end::Int, j::Int)
  for i=i_begin:i_end
    jy = u * ( ( lat.f[9, i, j] + lat.f[1, i, j] + lat.f[3, i, j] +
              2.0 * (lat.f[4, i, j] + lat.f[7, i, j] + lat.f[8, i, j]) ) 
              / (1.0 - u) );
    second_term = 1/6 * jy;
    third_term = 1/2 * (lat.f[1, i, j] - lat.f[3, i, j]);
    lat.f[2, i, j] = lat.f[4, i, j] + 2/3 * jy;
    lat.f[5, i, j] = lat.f[7, i, j] + second_term - third_term;
    lat.f[6, i, j] = lat.f[8, i, j] + second_term + third_term;
  end
end

#! South velocity boundary condition
function south_velocity!(lat::Lattice, u::Real)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  south_velocity!(lat, u, 1, ni, 1);
end

#! East velocity boundary condition
function east_velocity!(lat::LatticeD2Q9, u::Real, i::Int, j_begin::Int,
                        j_end::Int)

  for j=j_begin:j_end
    rhoe = (lat.f[9, i, j] + lat.f[2, i, j] + lat.f[4, i, j] +
            2.0 * (lat.f[1, i, j] + lat.f[5, i, j] + lat.f[8, i, j])) / (1.0 + u);
    lat.f[3, i, j] = lat.f[1, i, j] - 2.0 * rhoe * u / 3.0;
    second_term = rhoe * u / 6.0;
    third_term = 0.5 * (lat.f[2, i, j] - lat.f[4, i, j]);
    lat.f[7, i, j] = lat.f[5, i, j] - second_term + third_term;
    lat.f[6, i, j] = lat.f[6, i, j] - second_term - third_term;
  end
end

#! East velocity boundary condition
function east_velocity!(lat::Lattice, u::Real)
  east_velocity!(lat, u, size(lat.f, 2), 1, size(lat.f, 3));
end

#! West velocity boundary condition
function west_velocity!(lat::LatticeD2Q9, u::Real, i::Int, j_begin::Int,
  j_end::Int)

  for j=j_begin:j_end
    rhow = ((lat.f[9, i, j] + lat.f[2, i, j] + lat.f[4, i, j] +
             2.0 * (lat.f[3, i, j] + lat.f[6, i, j] + lat.f[7, i, j])) 
             / (1.0 - u));
    lat.f[1, i, j] = lat.f[3, i, j] + 2.0 * rhow * u / 3.0;
    second_term = rhow * u / 6.0;
    third_term = 0.5 * (lat.f[2, i, j] - lat.f[4, i, j]);
    lat.f[5, i, j] = lat.f[7, i, j] + second_term - third_term;
    lat.f[8, i, j] = lat.f[6, i, j] + second_term + third_term;
  end
end

#! West velocity boundary condition
function west_velocity!(lat::Lattice, u::Real)
  west_velocity!(lat, u, 1, 1, size(lat.f, 3));
end

#! Lid driven flow
function lid_driven!(lat::LatticeD2Q9, u::Real)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);

  for i=1:ni
    rho = (lat.f[9, i, nj] + lat.f[1, i, nj] + lat.f[3, i, nj] + 
           2.0 * (lat.f[2, i, nj] + lat.f[6, i, nj] + lat.f[5, i, nj]));
    lat.f[4, i, nj] = lat.f[2, i, nj];
    lat.f[8, i, nj] = lat.f[6, i, nj] + rho * u / 6.0;
    lat.f[7, i, nj] = lat.f[5, i, nj] - rho * u / 6.0;
  end
end

# =========================================================================== #
# ============================ pressure BCs ================================= #
# =========================================================================== #

#! Pressure north direction
function north_pressure!(lat::LatticeD2Q9, rhoo::Real)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);

  error("Check implementation details of north pressure BC.");

  for i=1:ni
    v = (-1. + (lat.f[9,i,nj] + lat.f[1,i,nj] + lat.f[3,i,nj]
         + 2. * (lat.f[2,i,nj] + lat.f[5,i,nj] + lat.f[6,i,nj])) / rhoo);
    ru = rhoo * v;
    lat.f[4,i,nj] = lat.f[2,i,nj] - (2./3.)*ru;
    lat.f[7,i,nj] = lat.f[5,i,nj] - (1./6.)*ru
                      + 0.5 * (lat.f[1,i,nj] - lat.f[3,i,nj]);
    lat.f[8,i,nj] = lat.f[6,i,nj] - (1./6.)*ru
                      + 0.5 * (lat.f[3,i,nj] - lat.f[1,i,nj]);
  end
end

#! Pressure north direction
function north_pressure!(lat::Lattice, rho::Real)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  north_pressure!(lat, rho, 1, ni, nj);
end

#! Pressure south direction
function south_pressure!(lat::LatticeD2Q9, rhoo::Real)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);

  error("Check implementation details of south pressure BC.");

  for i=1:ni
    v = (-1. + (lat.f[9,i,1] + lat.f[1,i,1] + lat.f[3,i,1]
         + 2. * (lat.f[4,i,1] + lat.f[7,i,1] + lat.f[8,i,1])) / rhoo);
    ru = rhoo * v;
    lat.f[2,i,1] = lat.f[4,i,1] - (2./3.)*ru;
    lat.f[5,i,1] = (lat.f[7,i,1] - (1./6.)*ru
                    + 0.5 * (lat.f[3,i,1] - lat.f[1,i,1]));
    lat.f[6,i,1] = (lat.f[8,i,1] - (1./6.)*ru
                    + 0.5 * (lat.f[1,i,1] - lat.f[3,i,1]));
  end
end

#! Pressure north direction
function south_pressure!(lat::Lattice, rho::Real)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  south_pressure!(lat, rho, 1, ni, 1);
end

#! Pressure east direction
function east_pressure!(lat::LatticeD2Q9, rho_out::Real, i::Int,
  j_begin::Int, j_end::Int)

  for j=j_begin:j_end
    u_x = (lat.f[9,i,j] + lat.f[2,i,j] + lat.f[4,i,j] +
           2 * (lat.f[1,i,j] + lat.f[5,i,j] + lat.f[8,i,j])) / rho_out - 1;
    lat.f[3,i,j] = lat.f[1,i,j] - 2/3 * rho_out * u_x;

    second_term = 0.5 * (lat.f[2,i,j] - lat.f[5,i,j]);
    third_term = 1/6 * rho_out * u_x;

    lat.f[7,i,j] = lat.f[5,i,j] + second_term - third_term;
    lat.f[6,i,j] = lat.f[8,i,j] - second_term - third_term;
  end
end

#! Pressure east direction
function east_pressure!(lat::LatticeD2Q9, rho_out::Real)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);

  # middle of inlet
  east_pressure!(lat, rho_out, ni, 2, nj-1);

  # bottom corner
  lat.f[3,ni,1] = lat.f[1,ni,1];
  lat.f[2,ni,1] = lat.f[4,ni,1];
  lat.f[7,ni,1] = lat.f[5,ni,1];
  f68 = 0.5 * (rho_out - (lat.f[9,ni,1] + lat.f[1,ni,1] + lat.f[2,ni,1]
          + lat.f[3,ni,1] + lat.f[4,ni,1] + lat.f[5,ni,1] + lat.f[7,ni,1]));
  lat.f[6,ni,1] = f68;
  lat.f[8,ni,1] = f68;

  # top corner
  lat.f[3,ni,nj] = lat.f[1,ni,nj];
  lat.f[4,ni,nj] = lat.f[2,ni,nj];
  lat.f[6,ni,nj] = lat.f[8,ni,nj];
  f57 = 0.5 * (rho_out - (lat.f[9,ni,nj] + lat.f[1,ni,nj] + lat.f[2,ni,nj]
          + lat.f[3,ni,nj] + lat.f[4,ni,nj] + lat.f[6,ni,nj] + lat.f[8,ni,nj]));
  lat.f[5,ni,nj] = f57;
  lat.f[7,ni,nj] = f57;
end

#! Pressure west direction
function west_pressure!(lat::LatticeD2Q9, rho_in::Real, i::Int,
                        j_begin::Int, j_end::Int)

  for j=j_begin:j_end
    u_x = (1 - (lat.f[9,i,j] + lat.f[2,i,j] + lat.f[4,i,j] +
           2 * (lat.f[3,i,j] + lat.f[6,i,j] + lat.f[7,i,j])) / rho_in);
    lat.f[1,i,j] = lat.f[3,i,j] + 2/3 * rho_in * u_x;

    second_term = 0.5 * (lat.f[2,i,j] - lat.f[4,i,j]);
    third_term = 1/6 * rho_in * u_x;

    lat.f[5,i,j] = lat.f[7,i,j] - second_term + third_term;
    lat.f[8,i,j] = lat.f[6,i,j] + second_term + third_term;
  end
end

#! Pressure west direction
function west_pressure!(lat::LatticeD2Q9, rho_in::Real)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);

  # middle
  west_pressure!(lat, rho_in, 1, 2, nj-1);

  # bottom corner
  lat.f[1,1,1] = lat.f[3,1,1];
  lat.f[2,1,1] = lat.f[4,1,1];
  lat.f[5,1,1] = lat.f[7,1,1];
  f68 = 0.5 * (rho_in - (lat.f[9,1,1] + lat.f[1,1,1] + lat.f[2,1,1]
          + lat.f[3,1,1] + lat.f[4,1,1] + lat.f[5,1,1] + lat.f[7,1,1]));
  lat.f[6,1,1] = f68;
  lat.f[8,1,1] = f68;

  # top corner
  lat.f[1,1,nj] = lat.f[4,1,nj];
  lat.f[1,1,nj] = lat.f[2,1,nj];
  lat.f[6,1,nj] = lat.f[8,1,nj];
  f57 = 0.5 * (rho_in - (lat.f[9,1,nj] + lat.f[1,1,nj] + lat.f[2,1,nj]
          + lat.f[3,1,nj] + lat.f[4,1,nj] + lat.f[6,1,nj] + lat.f[8,1,nj]));
  lat.f[5,1,nj] = f57;
  lat.f[7,1,nj] = f57;
end

# =========================================================================== #
# ================================ open BCs ================================= #
# =========================================================================== #

#! North open boundary
function north_open!(lat::LatticeD2Q9, i_begin::Int, i_end::Int, j::Int)
  for i=i_begin:i_end
    lat.f[6, i, j] = 2.0 * lat.f[6, i, j-1] - lat.f[6, i, j-2];
    lat.f[2, i, j] = 2.0 * lat.f[2, i, j-1] - lat.f[2, i, j-2];
    lat.f[5, i, j] = 2.0 * lat.f[5, i, j-1] - lat.f[5, i, j-2];
  end
end

#! North open boundary
function north_open!(lat::Lattice)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  north_open!(lat, 1, ni, nj);
end

#! South open boundary
function south_open!(lat::LatticeD2Q9, i_begin::Int, i_end::Int, j::Int)
  for i=i_begin:i_end
    lat.f[7, i, j] = 2.0 * lat.f[7, i, j+1] - lat.f[7, i, j+2];
    lat.f[4, i, j] = 2.0 * lat.f[4, i, j+1] - lat.f[4, i, j+2];
    lat.f[8, i, j] = 2.0 * lat.f[8, i, j+1] - lat.f[8, i, j+2];
  end
end

#! South open boundary
function south_open!(lat::Lattice)
  const ni = size(lat.f, 2);
  south_open!(lat, 1, ni, 1);
end

#! East open boundary
function east_open!(lat::LatticeD2Q9, i::Int, j_begin::Int, j_end::Int)
  for j=j_begin:j_end
    lat.f[1, i, j] = 2.0 * lat.f[1, i+1, j] - lat.f[1, i+2, j];
    lat.f[5, i, j] = 2.0 * lat.f[5, i+1, j] - lat.f[5, i+2, j];
    lat.f[8, i, j] = 2.0 * lat.f[8, i+1, j] - lat.f[8, i+2, j];
  end
end

#! East open boundary
function east_open!(lat::Lattice)
  const ni, nj = size(lat.f, 2), size(lat.f, 3);
  east_open!(lat, ni, 1, nj);
end

#! West open boundary
function west_open!(lat::LatticeD2Q9, i::Int, j_begin::Int, j_end::Int)
  for j=j_begin:j_end
    lat.f[6, i, j] = 2.0 * lat.f[6, i-1, j] - lat.f[6, i-2, j];
    lat.f[3, i, j] = 2.0 * lat.f[3, i-1, j] - lat.f[3, i-2, j];
    lat.f[7, i, j] = 2.0 * lat.f[7, i-1, j] - lat.f[7, i-2, j];
  end
end

#! West open boundary
function west_open!(lat::Lattice)
  const nj = size(lat.f, 3);
  west_open!(lat, 1, 1, nj);
end

include(joinpath("bcs", "mass.jl"));
