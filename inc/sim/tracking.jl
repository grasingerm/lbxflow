# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Simulate mass transfer across interface cells
#!
#! \param sim FreeSurfSim object
#! \param sbounds Bounds where fluid can be streaming
function masstransfer!(sim::FreeSurfSim, sbounds::Matrix{Int64})
  t   = sim.tracker;
  lat = sim.lat;
  msm = sim.msm;
  nbounds = size(sbounds, 2);

  for r = 1:nbounds
    i_min, i_max, j_min, j_max = sbounds[:,r];
    for j=j_min:j_max, i=i_min:i_max
      # TODO inbounds check is kind of redundant
      if t.state[i, j] != GAS && inbounds(i, j, sbounds)
        for k=1:lat.n-1
          i_nbr = i + lat.c[1, k];
          j_nbr = j + lat.c[2, k];
          if !inbounds(i_nbr, j_nbr, sbounds) || t.state[i_nbr, j_nbr] == GAS
            continue;
          elseif  t.state[i, j] == FLUID || t.state[i_nbr, j_nbr] == FLUID

            opk   = opp_lat_vec(lat, k);
            t.M[i, j]  += lat.f[opk, i_nbr, j_nbr] - lat.f[k, i, j] # m in - m out

          elseif  (t.state[i, j] == INTERFACE &&
                   t.state[i_nbr, j_nbr] == INTERFACE)

            opk   = opp_lat_vec(lat, k);
            t.M[i, j]  += ((t.eps[i, j] + t.eps[i_nbr, j_nbr]) / 2 *
                           (lat.f[opk, i_nbr, j_nbr] - lat.f[k, i, j]));

          else
            error("state not understood or invalid");
          end
        end
      end
    end
  end

end

#! Simulate mass transfer across interface cells
#!
#! \param   sim           FreeSurfSim object
#! \param   active_cells  Flags of active cells in the domain
function masstransfer!(sim::FreeSurfSim, active_cells::Matrix{Bool})
  t             =   sim.tracker;
  lat           =   sim.lat;
  msm           =   sim.msm;
  ni, nj  =   size(msm.rho);

  for j=1:nj, i=1:ni
    if t.state[i, j] != GAS && active_cells[i, j]
      for k=1:lat.n-1
        i_nbr = i + lat.c[1, k];
        j_nbr = j + lat.c[2, k];
        if (!inbounds(i_nbr, j_nbr, [1; ni; 1; nj]) ||
            !active_cells[i_nbr, j_nbr] || t.state[i_nbr, j_nbr] == GAS)
          continue;
        elseif  t.state[i, j] == FLUID || t.state[i_nbr, j_nbr] == FLUID

          opk   = opp_lat_vec(lat, k);
          t.M[i, j]  += lat.f[opk, i_nbr, j_nbr] - lat.f[k, i, j];

        elseif  (t.state[i, j] == INTERFACE &&
                 t.state[i_nbr, j_nbr] == INTERFACE)

          opk   = opp_lat_vec(lat, k);
          t.M[i, j]  += ((t.eps[i, j] + t.eps[i_nbr, j_nbr]) / 2 *
                         (lat.f[opk, i_nbr, j_nbr] - lat.f[k, i, j]));

        else
          error("state not understood or invalid");
        end
      end
    end
  end

end

# Update fluid fraction of each interface cell
function _update_fluid_fraction!(sim::FreeSurfSim, kappa::Real)
  t                 = sim.tracker;
  ni, nj      = size(t.state);
  new_empty_cells   = Set{Tuple{Int, Int}}();
  new_fluid_cells   = Set{Tuple{Int, Int}}();

  for (i, j) in t.interfacels
    t.eps[i, j] =   @_safe_calc_eps(t.M[i, j], sim.msm.rho[i, j]);

    if      t.M[i, j] > (1 + kappa) * sim.msm.rho[i, j]
      push!(new_fluid_cells, (i, j));
    elseif  t.M[i, j] < -kappa * sim.msm.rho[i, j]
      push!(new_empty_cells, (i, j));
    end

  end

  return new_empty_cells, new_fluid_cells;
end

#! Kernal function for changing the state of a cell
function _change_state!(t::Tracker, i::Int, j::Int, state::Fluid)
  t.state[i, j]   =   state;
  @assert(delete!(t.interfacels, (i, j)) != nothing,
          "Cell ($i, $j) not found in interface list. Unable to remove.");
end

#! Kernal function for changing the state of a cell
function _change_state!(t::Tracker, i::Int, j::Int, state::Gas)
  t.state[i, j]   =   state;
  @assert(delete!(t.interfacels, (i, j)) != nothing,
          "Cell ($i, $j) not found in interface list. Unable to remove.");
end

#! Kernal function for changing the state of a cell
function _change_state!(t::Tracker, i::Int, j::Int, state::Interface)
  t.state[i, j]   =   state;
  push!(t.interfacels, (i, j));
end

# TODO break this down into more kernal functions...
# Update cell states and redistribute mass to neighborhoods
function _update_cell_states!(sim::FreeSurfSim,
                              feq_f::LBXFunction,
                              new_empty_cells::Set{Tuple{Int, Int}},
                              new_fluid_cells::Set{Tuple{Int, Int}},
                              unorms::Dict{Tuple{Int, Int}, Vector{Float64}})

  lat           = sim.lat;
  msm           = sim.msm;
  t             = sim.tracker;

  ni, nj  = size(t.state);
  nk      = length(lat.w);

  excess_mass_after_redist = 0;

  for (i, j) in new_fluid_cells # First the neighborhood of all filled cells are prepared
    _change_state!(t, i, j, FLUID);

    # Calculate the total density and velocity of the neighborhood
    ρ_sum           =     0.0;
    u_sum           =     Float64[0.0; 0.0];
    counter::UInt   =     0;
    for k=1:nk-1
      i_nbr     =     i + lat.c[1, k];
      j_nbr     =     j + lat.c[2, k];

      if (i_nbr < 1 || i_nbr > ni || j_nbr < 1 || j_nbr > nj ||
          t.state[i_nbr, j_nbr] == GAS)
        continue;
      end

      counter         +=     1;

      ρ_sum           +=     msm.rho[i_nbr, j_nbr];
      u_sum           +=     msm.u[:, i_nbr, j_nbr];
    end

    if counter != 0.0
      ρ_avg            =     ρ_sum / counter;
      u_avg            =     u_sum / counter;
    else
      ρ_avg            =     1.0;
      u_avg            =     zeros(2);
    end

    # Construct interface cells from neighborhood average at equilibrium
    for k=1:nk-1
      i_nbr     =     i + lat.c[1, k];
      j_nbr     =     j + lat.c[2, k];

      if i_nbr < 1 || i_nbr > ni || j_nbr < 1 || j_nbr > nj
        continue;
      end

      # If it is already an interface cell, make sure it is not emptied
      if t.state[i_nbr, j_nbr] == INTERFACE
        delete!(new_empty_cells, (i_nbr, j_nbr));
      elseif t.state[i_nbr, j_nbr] == GAS
        _change_state!(t, i_nbr, j_nbr, INTERFACE);
        for kk=1:nk
          lat.f[kk, i_nbr, j_nbr]   =   feq_f(lat, ρ_avg, u_avg, kk);
        end
      end
    end # end reflag neighbors loop
  end

  for (i, j) in new_empty_cells # convert emptied cells to gas cells
    _change_state!(t, i, j, GAS);

    for k=1:nk-1
      i_nbr     =     i + lat.c[1, k];
      j_nbr     =     j + lat.c[2, k];

      if (i_nbr >= 1 && i_nbr <= ni && j_nbr >= 1 && j_nbr <= nj &&
          t.state[i_nbr, j_nbr] == FLUID)
        _change_state!(t, i_nbr, j_nbr, INTERFACE);
      end
    end
  end

  # Redistribute excess mass from new fluid cells
  for (i, j) in new_fluid_cells
    cells_to_redist_to  =     Set{Tuple{Int, Int, AbstractFloat}}();
    all_int_nbrs        =     Set{Tuple{Int, Int}}();
    #const unit_norm     =     _unit_normal(t, (i, j));
    unit_norm     =     unorms[(i, j)];
    v_sum               =     0.0;

    # Construct interface cells from neighborhood average at equilibrium
    for k=1:nk-1
      i_nbr     =     i + lat.c[1, k];
      j_nbr     =     j + lat.c[2, k];

      if (i_nbr >= 1 && i_nbr <= ni && j_nbr >= 1 && j_nbr <= nj &&
          t.state[i_nbr, j_nbr] == INTERFACE)
        push!(all_int_nbrs, (i_nbr, j_nbr));
        v_i        =     max(dot(unit_norm, lat.c[:, k]), 0.0);
        if v_i != 0.0
          push!(cells_to_redist_to, (i_nbr, j_nbr, v_i));
          v_sum           +=     v_i;
        end
      end
    end # find inteface loop

    # redistribute mass amoung valid neighbors
    mex = t.M[i, j] - msm.rho[i, j];
    if length(cells_to_redist_to) != 0 && v_sum != 0.0

      for (ii, jj, v_i) in cells_to_redist_to
        t.M[ii, jj]      +=   v_i / v_sum * mex;
        t.eps[ii, jj]     =   @_safe_calc_eps(t.M[ii, jj], msm.rho[ii, jj]);
      end

    elseif length(all_int_nbrs) != 0

      for (ii, jj) in all_int_nbrs
        t.M[ii, jj]      +=   mex / length(all_int_nbrs);
        t.eps[ii, jj]     =   @_safe_calc_eps(t.M[ii, jj], msm.rho[ii, jj]);
      end

    else
      excess_mass_after_redist += mex;
    end

    t.M[i, j]   = msm.rho[i, j]; # set mass to local ρ
    t.eps[i, j] = 1.0;

  end

  # TODO consider putting mass distribution in a kernal function
  # Redistribute excess mass from emptied cells
  for (i, j) in new_empty_cells
    cells_to_redist_to  =     Set{Tuple{Int, Int, AbstractFloat}}();
    all_int_nbrs        =     Set{Tuple{Int, Int}}();
    #const unit_norm     =     _unit_normal(t, (i, j));
    unit_norm     =     unorms[(i, j)];
    v_sum               =     0.0;

    # Construct interface cells from neighborhood average at equilibrium
    for k=1:nk-1
      i_nbr     =     i + lat.c[1, k];
      j_nbr     =     j + lat.c[2, k];

      if (i_nbr >= 1 && i_nbr <= ni && j_nbr >= 1 && j_nbr <= nj &&
          t.state[i_nbr, j_nbr] == INTERFACE)
        push!(all_int_nbrs, (i_nbr, j_nbr));
        v_i        =     -min(dot(unit_norm, lat.c[:, k]), 0.0);
        if v_i != 0.0
          push!(cells_to_redist_to, (i_nbr, j_nbr, v_i));
          v_sum           +=     v_i;
        end
      end
    end # find inteface loop
    
    # redistribute mass amoung valid neighbors
    mex = t.M[i, j];
    if length(cells_to_redist_to) != 0 && v_sum != 0.0

      for (ii, jj, v_i) in cells_to_redist_to
        t.M[ii, jj]      +=   v_i / v_sum * mex;
        t.eps[ii, jj]     =   @_safe_calc_eps(t.M[ii, jj], msm.rho[ii, jj]);
      end

    elseif length(all_int_nbrs) != 0

      for (ii, jj) in all_int_nbrs
        t.M[ii, jj]      +=   mex / length(all_int_nbrs);
        t.eps[ii, jj]     =   @_safe_calc_eps(t.M[ii, jj], msm.rho[ii, jj]);
      end

    else
      excess_mass_after_redist += mex;
    end

    t.M[i, j]   = 0.0;
    t.eps[i, j] = 0.0;

  end

  # redistribute remaining mass
  dm  =   excess_mass_after_redist / length(t.interfacels);
  if !isnan(dm)
    for (i, j) in t.interfacels
      t.M[i, j]   +=  dm;
      t.eps[i, j]  =  @_safe_calc_eps(t.M[i, j], msm.rho[i, j]); 
    end
  end
end

#! Update cell states of tracker
#!
#! \param   sim   FreeSurfSim object
#! \param   feq_f Equilibrium distribution function
#! \param   kappa Constant for determining threshold of state change
function update!(sim::FreeSurfSim, feq_f::LBXFunction, 
                 unorms::Dict{Tuple{Int, Int}, Vector{Float64}}, kappa=1.0e-3)
  new_empty_cells, new_fluid_cells  =   _update_fluid_fraction!(sim, kappa);
  _update_cell_states!(sim, feq_f, new_empty_cells, new_fluid_cells, unorms);
end

#! Kernal function for interface f reconstruction
#!
#! \param   lat     Lattice
#! \param   msm     Multiscale map
#! \param   feq_f   Equilibrium distribution function
#! \param   rho_g   Gas pressure (in lattice units)
#! \param   k       Lattice vector direction of neighbor
#! \param   i       Index of node in ith direction
#! \param   j       Index of node in jth direction
function _f_reconst_ij!(lat::Lattice, msm::MultiscaleMap, feq_f::LBXFunction,
                        rho_g::Real, k::Int, i::Int, j::Int)
  opp_k         =   opp_lat_vec(lat, k);
  u_ij          =   msm.u[:, i, j];
  lat.f[opp_k, i, j]  =   (feq_f(lat, rho_g, u_ij, k) +
                           feq_f(lat, rho_g, u_ij, opp_k) -
                           lat.f[k, i, j]);
end

#! Reconstruct distribution functions at interface
#!
#! \param sim       Free surface simulation object
#! \param t         Mass tracker
#! \param ij        Tuple of i and j indices on grid
#! \param feq_f     Equilibrium distribution function
#! \param rho_g     Atmospheric pressure
function f_reconst!(sim::FreeSurfSim, t::Tracker, ij::Tuple{Int64, Int64},
                    feq_f::LBXFunction, rho_g::Real)
  n     = _unit_normal(t, ij);
  i, j  = ij;
  lat         = sim.lat;
  msm         = sim.msm;

  # reconstruct f from gas cells
  for k = 1:lat.n-1
    i_nbr   =   i + lat.c[1,k];
    j_nbr   =   j + lat.c[2,k];

    if (i_nbr < 1 || i_nbr > ni || j_nbr < 1 || j_nbr > nj ||
        t.state[i_nbr, j_nbr] == GAS)
      _f_reconst_ij!(lat, msm, feq_f, rho_g, k, i, j);
    end
  end

  # reconstruct f normal to interface (conservation of linear momentum)
  for k = 1:lat.n-1
    c = lat.c[:,k]; #TODO consider using an ArrayView here
    if dot(n, c) > 0
      # NOTE: reconstructs direction OPPOSITE of k (as desired)
      _f_reconst_ij!(lat, msm, feq_f, rho_g, k, i, j);
    end
  end

  return n;
end

#! Determine unit normal to interface using finite difference
#!
#! \param t Mass tracker
#! \param ij i and j indices of grid cell
function _unit_normal_fd(t::Tracker, ij::Tuple{Int64, Int64})
  ni, nj      =   size(t.state);
  i, j        =   ij;

  # TODO consider using heuristic to tell if we should attach to wall
  eps_il      =   (i <= 1)  ? 0.0 : t.eps[i-1, j];
  eps_ir      =   (i >= ni) ? 0.0 : t.eps[i+1, j];
  eps_jd      =   (j <= 1)  ? 0.0 : t.eps[i, j-1];
  eps_ju      =   (j >= nj) ? 0.0 : t.eps[i, j+1];

  return 1/2 * Float64[eps_il - eps_ir; eps_jd - eps_ju];
end

#! Determine unit normal to interface using a type of marching squares
#!
#! \param t Mass tracker
#! \param ij i and j indices of grid cell
function _unit_normal_ms(t::Tracker, ij::Tuple{Int64, Int64})
  ni, nj = size(t.state);
  i, j = ij;

  CORNERS = ( ((-1,0),(-1,-1),(0,-1)),  ((0,-1),(1,-1),(1,0)),
                  ((1,0),(1,1),(0,1)      ),  ((0,1),(-1,1),(-1,0))   );
  BITS    = (1, 2, 4, 8);
  NORMALS = (
                    Float64[0; 0],                   Float64[-sqrt(2)/2; -sqrt(2)/2],
                    Float64[sqrt(2)/2; -sqrt(2)/2],  Float64[0;-1],
                    Float64[sqrt(2)/2; sqrt(2)/2],   Float64[0; 0],
                    Float64[1; 0],                   Float64[-sqrt(2)/2; -sqrt(2)/2],
                    Float64[-sqrt(2)/2; sqrt(2)/2],  Float64[1; 0],
                    Float64[0; 0],                   Float64[-sqrt(2)/2; -sqrt(2)/2],
                    Float64[0; 1],                   Float64[-sqrt(2)/2; sqrt(2)/2],
                    Float64[sqrt(2)/2; sqrt(2)/2],   Float64[0; 0]
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
  end

  case += 1 # correct for indexing starting at 1
  return NORMALS[case]; # TODO: increase accuracy with linear interpolation
end

_unit_normal = _unit_normal_fd;

function _get_nbrhd(u, i::Int, j::Int)
  ni, nj = size(u);

  return ([ (i-1 > 0 && j+1 <= nj) ? u[i-1,j+1] : "N/A"
            (j+1 <= nj) ? u[i,j+1] : "N/A"
            (i+1 <= ni && j+1 <= nj) ? u[i+1,j+1] : "N/A";
            (i-1 > 0) ? u[i-1, j] : "N/A"
            u[i, j]
            (i+1 <= ni) ? u[i+1, j] : "N/A";
            (i-1 > 0 && j-1 > 0) ? u[i-1, j-1] : "N/A"
            (j-1 > 0) ? u[i, j-1] : "N/A"
            (i+1 <= ni && j-1 > 0) ? u[i+1, j-1] : "N/A"
          ]);
end
