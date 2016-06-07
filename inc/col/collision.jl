# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

macro _scale_omega(omega, st)
  return :(1.0 / ($st * (1/$omega - 0.5) + 0.5));
end

macro _scale_f(f, st)
  return :($st*$st * $f);
end

macro _scale_mu(mu, st)
  return :($st * $mu);
end

include("forcing.jl");
include("equilibrium.jl");
include("constitutive.jl");
include("mrt_matrices.jl");
include("freecol.jl"); # free surface collisions
include("modcol.jl");  # "modular" collisions
include("pmodcol.jl"); # parallel mod collisions
include("bgk.jl");
include("mrt.jl");
include("entropic.jl");
include("filtering.jl");
