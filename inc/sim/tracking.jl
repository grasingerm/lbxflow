# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

using ArrayViews

#! Check to see if a pair of indices fall inside a bounds
function inbounds(i::Int, j::Int, sbounds::Matrix{Int64})
  const nbounds = size(sbounds, 2);
  for n=1:nbounds
    if i < sbounds[1,n]; return false; end
    if i > sbounds[2,n]; return false; end
    if j < sbounds[3,n]; return false; end
    if j > sbounds[4,n]; return false; end
  end
  return true;
end

#! Simulate mass transfer across interface cells
#!
#! \param sim FreeSurfSim object
#! \param sbounds Bounds where fluid can be streaming
function masstransfer!(sim::FreeSurfSim, sbounds::Matrix{Int64})
  t   = sim.tracker;
  lat = sim.lat;
  msm = sim.msm;

  for node in t.interfacels
    i, j = node.val;
    if !inbounds(i, j, sbounds)
      continue;
    end
    for k=1:lat.n
      i_nbr = i + lat.c[1,k];
      j_nbr = j + lat.c[2,k];
      if !inbounds(i_nbr, j_nbr, sbounds) || t.state[i_nbr,j_nbr] == GAS
        continue;
      elseif t.state[i,j] == FLUID && t.state[i_nbr,j_nbr] == FLUID

        opk = opp_lat_vec(lat, k);
        t.M[i,j] += lat.f[opk,i_nbr,j_nbr] - lat.f[k,i,j] # m in - m out

      elseif ((t.state[i,j] == FLUID ||  && t.state[i_nbr,j_nbr] == INTERFACE) ||
              (t.state[i,j] == INTERFACE && t.state[i_nbr,j_nbr] == FLUID)     ||
              (t.state[i,j] == INTERFACE && t.state[i_nbr,j_nbr] == INTERFACE)) 

        opk = opp_lat_vec(lat, k);
        epsx = t.M[i,j] / msm.rho[i,j];
        epsxe = t.M[i_nbr,j_nbr] / msm.rho[i_nbr,j_nbr];
        @assert(epsx > -0.5 && epsx <= 1.5, "Fluid fraction at ($i,$j) not " *
                "within bounds. ϵ($i,$j) = $epsx");
        @assert(epsxe > -0.5 && epsxe <= 1.5, "Fluid fraction at "
                "($i_nbr,$j_nbr) not within bounds. ϵ($i_nbr,$j_nbr) = $epsxe");
        t.M[i,j] += (epsx + epsxe) / 2 * (lat.f[opk,i_nbr,j_nbr] - lat.f[k,i,j]);

      else
        error("state not understood or invalid");
      end
    end
  end
end

#! Update cell states of tracker
#!
#! \param sim FreeSurfSim object
#! \param sbounds Bounds where fluid can be streaming
function update!(sim::FreeSurfSim, sbounds::Matrix{Int64},
                 threshold::AbstractFloat=1.0e-3)
  t = sim.tracker;
  msm = sim.msm;

  for node in t.interfacels
    i, j = node.val;
    if !inbounds(i, j, sbounds)
      continue;
    end
    
    if t.M[i,j] < -threshold
      println("node $i,$j is becoming a gas cell");
      t.M[i,j] = 0;
      t.state[i,j] = GAS;
      remove!(t.interfacels, node); # no longer an interface cell
      for ci in (-1, 1), cj in (-1, 1)
        i_nbr = i+ci;
        j_nbr = j+cj;
        if !inbounds(i_nbr,j_nbr,sbounds); continue; end
        if t.state[i_nbr,j_nbr] == FLUID
          t.state[i_nbr,j_nbr] = INTERFACE;
          println("node $i_nbr,$j_nbr is becoming an interface cell");
          unshift!(t.interfacels, (i_nbr,j_nbr)); # push into interface list
        end
      end
      for ci in (-1, 1)
        i_nbr = i+ci;
        if !inbounds(i_nbr,j,sbounds); continue; end
        if t.state[i_nbr,j] == FLUID
          t.state[i_nbr,j] = INTERFACE;
          println("node $i_nbr,$j is becoming an interface cell");
          unshift!(t.interfacels, (i_nbr, j)); # push into interface list
        end
      end
      for cj in (-1, 1)
        j_nbr = j+cj;
        if !inbounds(i,j_nbr,sbounds); continue; end
        if t.state[i,j_nbr] == FLUID
          t.state[i,j_nbr] = INTERFACE;
          println("node $i,$j_nbr is becoming an interface cell");
          unshift!(t.interfacels, (i, j_nbr)); # push into interface list
        end
      end
    elseif t.M[i,j] > (msm.rho[i,j] + threshold)
      t.M[i,j] = msm.rho[i,j];
      t.state[i,j] = FLUID;
      remove!(t.interfacels, node); # no longer an interface cell
      println("node $i,$j is becoming a fluid cell");
      for ci in (-1, 1), cj in (-1, 1)
        i_nbr = i+ci;
        j_nbr = j+cj;
        if !inbounds(i_nbr,j_nbr,sbounds); continue; end
        if t.state[i_nbr,j_nbr] == GAS
          t.state[i_nbr,j_nbr] = INTERFACE;
          println("node $i_nbr,$j_nbr is becoming an interface cell");
          unshift!(t.interfacels, (i_nbr,j_nbr)); # push into interface list
        end
      end
      for ci in (-1, 1)
        i_nbr = i+ci;
        if !inbounds(i_nbr,j,sbounds); continue; end
        if t.state[i_nbr,j] == GAS
          t.state[i_nbr,j] = INTERFACE;
          println("node $i_nbr,$j is becoming an interface cell");
          unshift!(t.interfacels, (i_nbr, j)); # push into interface list
        end
      end
      for cj in (-1, 1)
        j_nbr = j+cj;
        if !inbounds(i,j_nbr,sbounds); continue; end
        if t.state[i,j_nbr] == GAS
          t.state[i,j_nbr] = INTERFACE;
          println("node $i,$j_nbr is becoming an interface cell");
          unshift!(t.interfacels, (i, j_nbr)); # push into interface list
        end
      end
    end
  end
end

#! Reconstruct distribution functions at interface
#!
#! \param sim Free surface simulation object
#! \param t Mass tracker
#! \param ij tuple of i and j indices on grid
#! \param sbounds Bounds where fluid can be streaming
#! \param rhog Atmospheric pressure
function f_reconst!(sim::FreeSurfSim, t::Tracker, ij::Tuple{Int64, Int64},
                    sbounds::Matrix{Int64}, rhog::AbstractFloat)
  const n = unit_normal(t, ij);
  const i, j = ij;
  lat = sim.lat;
  msm = sim.msm;

  for k = 1:lat.n
    c = lat.c[:,k]; # TODO: consider using an ArrayView here
    if dot(n, c) >= 0
      opk = opp_lat_vec(lat, k);
      rhoij = msm.rho[i,j];
      uij = msm.u[:,i,j]; # TODO: consider using an ArrayView here
      lat.f[k,i,j] = (feq_incomp(lat, rhog, uij, k) -
                  feq_incomp(lat, rhog, uij, opk) -
                  lat.f[opk,i,j]);
    else
      i_new = i - c[1];
      j_new = j - c[2];
      if !inbounds(i_new, j_new, sbounds); continue; end
      lat.f[k,i,j] = lat.f[k,i_new,j_new];
    end
  end
end

#! Determine unit normal to interface
#!
#! \param t Mass tracker
#! \param ij i and j indices of grid cell
function unit_normal(t::Tracker, ij::Tuple{Int64, Int64})
  const ni, nj = size(t.state);
  const i, j = ij;

  const CORNERS = ( ((-1,0),(-1,-1),(0,-1)),  ((0,-1),(1,-1),(1,0)),
                    ((1,0),(1,1),(0,1)    ),  ((0,1),(-1,1),(-1,0))   );
  const BITS    = (1, 2, 4, 8);
  const NORMALS = (
                    [0; 0],                   [-sqrt(2)/2; -sqrt(2)/2],
                    [sqrt(2)/2; -sqrt(2)/2],  [0;-1],
                    [sqrt(2)/2; sqrt(2)/2],   [0; 0],
                    [1; 0],                   [-sqrt(2)/2; -sqrt(2)/2],
                    [-sqrt(2)/2; sqrt(2)/2],  [1; 0],
                    [0; 0],                   [-sqrt(2)/2; -sqrt(2)/2],
                    [0; 1],                   [-sqrt(2)/2; sqrt(2)/2],
                    [sqrt(2)/2; sqrt(2)/2],   [0; 0]
                  );

  case = 0;
  # are any of the neighbors on a corner a GAS cell??
  for (bit, corner) in zip(BITS, CORNERS)
    is_fluid = true;
    for (ci,cj) in corner
      i_new = i+ci;
      j_new = j+cj;
      if (i_new > ni || j_new > nj || i_new < 1 || j_new < 1
          || t.state[i+ci,j+cj] == GAS)
        is_fluid = false;
        break;
      end
    end

    if is_fluid; case |= bit; end
  end

  case += 1 # correct for indexing starting at 1
  return NORMALS[case]; # TODO: increase accuracy with linear interpolation
end
