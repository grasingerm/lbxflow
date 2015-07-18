# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Compute the double dot product of two matrices
function ddot(A::Array{Float64,2}, B::Array{Float64,2})
  sum = 0;
  for (a,b) in zip(A,B)
    sum += a*b;
  end
  return sum;
end
