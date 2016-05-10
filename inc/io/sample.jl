# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

function sample_u_uniform_sqr(i_begin::Int, j_begin::Int, nside::Int,
                              nsidepts::Int)
  return (sim::AbstractSim) -> begin
    dn    =   convert(Int, floor(nside / nsidepts));
    vals  =   Vector{Float64}(nsidepts * nsidepts);

    i       =   i_begin;
    j       =   j_begin;
    ins_at  =   1;

    for n = 1:nsidepts
      for n = 1:nsidepts
        vals[ins_at]  =   sim.msm.u[1, i, j];
        i       +=  dn;
        ins_at  +=  1;
      end
      i    =   i_begin;
      j   +=   dn;
    end

    return vals;
  end;
end

function sample_v_uniform_sqr(i_begin::Int, j_begin::Int, nside::Int,
                              nsidepts::Int)
  return (sim::AbstractSim) -> begin
    dn    =   convert(Int, floor(nside / nsidepts));
    vals  =   Vector{Float64}(nsidepts * nsidepts);

    i       =   i_begin;
    j       =   j_begin;
    ins_at  =   1;

    for n = 1:nsidepts
      for n = 1:nsidepts
        vals[ins_at]  =   sim.msm.u[2, i, j];
        i       +=  dn;
        ins_at  +=  1;
      end
      i    =   i_begin;
      j   +=   dn;
    end

    return vals;
  end;
end

function sample_ρ_uniform_sqr(i_begin::Int, j_begin::Int, nside::Int,
                              nsidepts::Int)
  return (sim::AbstractSim) -> begin
    dn    =   convert(Int, floor(nside / nsidepts));
    vals  =   Vector{Float64}(nsidepts * nsidepts);

    i       =   i_begin;
    j       =   j_begin;
    ins_at  =   1;

    for n = 1:nsidepts
      for n = 1:nsidepts
        vals[ins_at]  =   sim.msm.rho[i, j];
        i       +=  dn;
        ins_at  +=  1;
      end
      i    =   i_begin;
      j   +=   dn;
    end

    return vals;
  end;
end
