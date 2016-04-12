# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

macro NDEBUG()
  return false;
end

macro NDEBUG_MASS()
  return false;
end

macro mdebug(message)
  if @NDEBUG()
    return;
  else
    return :(warn($message));
  end
end

macro checkdebug(condition, message)
  if @NDEBUG()
    return;
  else
    return quote
      if !$condition
        warn($message);
      end
    end;
  end
end

# Debugging macros
macro _mdebug_mass_cons(opname, M, block)
  if @NDEBUG_MASS()
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
  if @NDEBUG_MASS()
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
