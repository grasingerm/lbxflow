# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

_NDEBUG       =   false;
_NDEBUG_MASS  =   true;

function turn_off_debugging(); global _NDEBUG; _NDEBUG = true; end
function turn_off_mass_cons_debugging()
  global _NDEBUG_MASS;
  _NDEBUG_MASS = true;
end

macro mdebug(message)
  if _NDEBUG
    return;
  else
    return :(@warn($message));
  end
end

macro checkdebug(condition, message)
  if _NDEBUG
    return;
  else
    return quote
      if !$condition
        @warn($message);
      end
    end;
  end
end

# Debugging macros
macro _mdebug_mass_cons(opname, M, block)
  if _NDEBUG_MASS
    return block;
  else
    return quote
      _init_mass = sum($M);
      $block;
      @mdebug($opname * ": ΔM = $(sum($M) - _init_mass)");
    end
  end
end

# Debugging macros
macro _checkdebug_mass_cons(opname, M, block, eps)
  if _NDEBUG_MASS
    return block;
  else
    return quote
      _init_mass = sum($M);
      $block;
      @checkdebug(abs(sum($M) - _init_mass) < $eps,
                  $opname * ": ΔM = $(sum($M) - _init_mass)");
    end
  end
end
