# SPDX-License-Identifier: GPL-2.0-or-later
# https://github.com/dehesselle/bash_d

include_guard

### description ################################################################

# This file contains ANSI control codes to control terminal output.

### includes ###################################################################

# Nothing here.

### variables ##################################################################

ANSI_ENABLED=true       # allow ANSI escape sequences for e.g. colors...
ANSI_TERM_ONLY=true     # ...but only when running in a terminal

ANSI_CURSOR_LEFT="\033[1D"

#            0: normal
#            1: bold
#            2: faint
#            |
#       \033[1;32m
#               |
#               color

ANSI_FG_BLACK="\033[0;30m"
ANSI_FG_BLACK_BOLD="\033[1;30m"
ANSI_FG_BLACK_BRIGHT="\033[0;90m"
ANSI_FG_BLUE="\033[0;34m"
ANSI_FG_BLUE_BOLD="\033[1;34m"
ANSI_FG_BLUE_BRIGHT="\033[0;94m"
ANSI_FG_CYAN="\033[0;36m"
ANSI_FG_CYAN_BOLD="\033[1;36m"
ANSI_FG_CYAN_BRIGHT="\033[0;96m"
ANSI_FG_GREEN="\033[0;32m"
ANSI_FG_GREEN_BOLD="\033[1;32m"
ANSI_FG_GREEN_BRIGHT="\033[0;92m"
ANSI_FG_MAGENTA="\033[0;35m"
ANSI_FG_MAGENTA_BOLD="\033[1;35m"
ANSI_FG_MAGENTA_BRIGHT="\033[0;95m"
ANSI_FG_RED="\033[0;31m"
ANSI_FG_RED_BOLD="\033[1;31m"
ANSI_FG_RED_BRIGHT="\033[0;91m"
ANSI_FG_RESET="\033[0;0m"
ANSI_FG_WHITE="\033[0;37m"
ANSI_FG_WHITE_BOLD="\033[1;37m"
ANSI_FG_WHITE_BRIGHT="\033[0;97m"
ANSI_FG_YELLOW="\033[0;33m"
ANSI_FG_YELLOW_BOLD="\033[1;33m"
ANSI_FG_YELLOW_BRIGHT="\033[0;93m"

### functions ##################################################################

# Nothing here.

### aliases ####################################################################

# This performs the following check:
#   - usage of ANSI is generally enabled AND
#   - we are either running in a terminal OR we don't care at all
alias ansi_is_usable='$ANSI_ENABLED && ([ -t 1 ] || ! $ANSI_TERM_ONLY)'

### main #######################################################################

# Nothing here.