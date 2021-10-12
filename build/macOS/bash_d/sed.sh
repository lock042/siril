# SPDX-License-Identifier: GPL-2.0-or-later
# https://github.com/dehesselle/bash_d

include_guard

### includes ###################################################################

# Nothing here.

### variables ##################################################################

# Nothing here.

### functions ##################################################################

# Escape slashes, backslashes and ampersands in strings to be used s as
# replacement strings when using 'sed'. Newlines are not taken into
# consideartion here.
# https://stackoverflow.com/a/2705678
function sed_escape_str
{
  local string="$*"

  echo "$string" | sed -e 's/[\/&]/\\&/g'
}

### aliases ####################################################################

# Nothing here.

### main #######################################################################

# Nothing here.