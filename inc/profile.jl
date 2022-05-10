# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Profile;

macro profif(condition, statement)
  quote
    if $condition
      @profile $statement
    else
      $statement
    end
  end
end
