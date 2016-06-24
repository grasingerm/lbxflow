# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Test if flow is steady state
#!
#! \param msm Current multiscale map
#! \param prev_msm Previous multiscale map
#! \param tol Threshold for determining steady state
#! \return Flag for sim termination
function is_steadystate(msm::MultiscaleMap, prev_msm::MultiscaleMap,
                        tol::AbstractFloat = 5.0e-7)

  sum_diff = 0.0;
  sum_u = 0.0;

  for (u, u_prev) in zip(msm.u, prev_msm.u)
    sum_diff += abs(u - u_prev);
    sum_u += abs(u);
  end

  if sum_diff / sum_u <= tol
    return true;
  end

  return false;

end

#! Initialize is steadystate function with a specified tolerance
function init_is_steadystate(tol::AbstractFloat)
  return (current, prev) -> is_steadystate(current, prev, tol);
end

#! Test if flow is steady state
#!
#! \param msm Current multiscale map
#! \param prev_msms Previous multiscale map
#! \param tol Threshold for determining steady state
#! \return Flag for sim termination
function is_steadystate(msm::MultiscaleMap, prev_msms::Vector{MultiscaleMap},
                        tol::AbstractFloat = 5.0e-7)

  sum_diff = 0.0;
  sum_u = 0.0;

  for prev_msm in prev_msms
    for (u, u_prev) in zip(msm.u, prev_msm.u)
      sum_diff += abs(u - u_prev);
      sum_u += abs(u);
    end
  end

  if sum_diff / sum_u <= tol
    return true;
  end

  return false;

end

#! Test if flow is steady state in the x-direction
#!
#! \param msm Current multiscale map
#! \param prev_msm Previous multiscale map
#! \param tol Threshold for determining steady state
#! \return Flag for sim termination
function is_steadystate_x(msm::MultiscaleMap, prev_msm::MultiscaleMap,
                          tol::AbstractFloat = 5.0e-7)

  sum_diff = 0.0;
  sum_u = 0.0;

  for (u, u_prev) in zip(msm.u[1,:,:], prev_msm.u[1,:,:])
    sum_diff += abs(u - u_prev);
    sum_u += abs(u);
  end

  if sum_diff / sum_u <= tol
    return true;
  end

  return false;

end

#! Test if flow is steady state
#!
#! \param msm Current multiscale map
#! \param prev_msms Previous multiscale map
#! \param tol Threshold for determining steady state
#! \return Flag for sim termination
function is_steadystate_x(msm::MultiscaleMap, prev_msms::Vector{MultiscaleMap},
                          tol::AbstractFloat = 5.0e-7)

  sum_diff = 0.0;
  sum_u = 0.0;

  for prev_msm in prev_msms
    for (u, u_prev) in zip(msm.u[1,:,:], prev_msm.u[1,:,:])
      sum_diff += abs(u - u_prev);
      sum_u += abs(u);
    end
  end

  if sum_diff / sum_u <= tol
    return true;
  end

  return false;

end

#! Test if flow is steady state in the y-direction
#!
#! \param msm Current multiscale map
#! \param prev_msm Previous multiscale map
#! \param tol Threshold for determining steady state
#! \return Flag for sim termination
function is_steadystate_y(msm::MultiscaleMap, prev_msm::MultiscaleMap,
                          tol::AbstractFloat = 5.0e-7)

  sum_diff = 0.0;
  sum_u = 0.0;

  for (u, u_prev) in zip(msm.u[2,:,:], prev_msm.u[2,:,:])
    sum_diff += abs(u - u_prev);
    sum_u += abs(u);
  end

  if sum_diff / sum_u <= tol
    return true;
  end

  return false;

end

#! Test if flow is steady state
#!
#! \param msm Current multiscale map
#! \param prev_msms Previous multiscale map
#! \param tol Threshold for determining steady state
#! \return Flag for sim termination
function is_steadystate_y(msm::MultiscaleMap, prev_msms::Vector{MultiscaleMap},
                          tol::AbstractFloat = 5.0e-7)

  sum_diff = 0.0;
  sum_u = 0.0;

  for prev_msm in prev_msms
    for (u, u_prev) in zip(msm.u[2,:,:], prev_msm.u[2,:,:])
      sum_diff += abs(u - u_prev);
      sum_u += abs(u);
    end
  end

  if sum_diff / sum_u <= tol
    return true;
  end

  return false;

end
