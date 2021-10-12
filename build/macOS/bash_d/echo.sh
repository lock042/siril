# SPDX-License-Identifier: GPL-2.0-or-later
# https://github.com/dehesselle/bash_d

include_guard

### includes ###################################################################

include_file ansi.sh

### variables ##################################################################

# Nothing here.

### functions ##################################################################

function _echo_message
{
  local name=$1   # can be empty
  local type=$2
  local color=$3
  local args=${@:4}

  if [ ! -z $name ]; then
    name=" $name"
  fi

  if ansi_is_usable; then
    echo -e "$color[$type]$name$ANSI_FG_RESET $args"
  else
    echo "[$type]$name $args"
  fi
}

### aliases ####################################################################

alias echo_d='>&2 _echo_message "$FUNCNAME" " debug " "$ANSI_FG_BLACK_BOLD"'
alias echo_e='>&2 _echo_message "$FUNCNAME" " error " "$ANSI_FG_RED_BOLD"'
alias echo_i='>&2 _echo_message "$FUNCNAME" " info  " "$ANSI_FG_BLUE_BOLD"'
alias echo_o='>&2 _echo_message "$FUNCNAME" "  ok   " "$ANSI_FG_GREEN_BOLD"'
alias echo_w='>&2 _echo_message "$FUNCNAME" "warning" "$ANSI_FG_YELLOW_BOLD"'

### main #######################################################################

# Nothing here.
