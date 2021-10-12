# SPDX-License-Identifier: GPL-2.0-or-later
# https://github.com/dehesselle/bash_d

include_guard

### description ################################################################

# This file provides a mechanism to restrict sourcing a script to specific
# platforms.

### includes ###################################################################

include_file echo.sh

### variables ##################################################################

[ "$(uname)" = "Darwin" ] && PLATFORM_DARWIN=true || PLATFORM_DARWIN=false
[ "$(uname)" = "Linux"  ] && PLATFORM_LINUX=true  || PLATFORM_LINUX=false

### functions ##################################################################

# Nothing here.

### aliases ####################################################################

alias platform_darwin_only='$PLATFORM_DARWIN && true || return'
alias platform_linux_only='$PLATFORM_LINUX && true || return'

alias platform_bash4_and_above_only='if [ ${BASH_VERSINFO[0]} -ge 4 ]; then : ; else echo_e "bash 4.x or above required"; return 1; fi'

### main #######################################################################

# Nothing here.
