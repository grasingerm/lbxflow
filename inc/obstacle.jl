# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

abstract Obstacle;

type BounceBack <: Obstacle
  i_min::Int;
  i_max::Int;
  j_min::Int;
  j_max::Int;

  BounceBack(i_min::Int, i_max::Int, j_min::Int, j_max::Int) = (
                                            new(i_min, i_max, j_min, j_max));
end

# Boundary condition for bounce back obstacle
function call(bb::BounceBack, lat::Lattice)
  north_bounce_back!(lat, bb.i_min, bb.i_max, bb.j_min);
  south_bounce_back!(lat, bb.i_min, bb.i_max, bb.j_max);
  east_bounce_back!(lat, bb.i_min, bb.j_min, bb.j_max);
  west_bounce_back!(lat, bb.i_max, bb.j_min, bb.j_max);
end

type Reflector <: Obstacle
  i_min::Int;
  i_max::Int;
  j_min::Int;
  j_max::Int;

  Reflector(i_min::Int, i_max::Int, j_min::Int, j_max::Int) = (
                                            new(i_min, i_max, j_min, j_max));
end

# Boundary condition for reflecting obstacle
function call(r::Reflector, lat::Lattice)
  north_reflect!(lat, r.i_min, r.i_max, r.j_min);
  south_reflect!(lat, r.i_min, r.i_max, r.j_max);
  east_reflect!(lat, r.i_min, r.j_min, r.j_max);
  west_reflect!(lat, r.i_max, r.j_min, r.j_max);
end
