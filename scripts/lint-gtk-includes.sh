#!/bin/sh
# lint-gtk-includes.sh — enforce the GTK-free boundary in non-GUI source files.
#
# Fails if any file outside the GUI compilation units includes <gtk/gtk.h> or
# <gdk/gdk.h> directly.  These headers pull in the entire GTK/GDK toolkit and
# must stay confined to:
#   src/gui/          — all GUI widget code
#   src/main.c        — GTK application entry point
#   src/livestacking/gui.c — livestacking GUI (in src_files_gui)
#   src/core/siril_cmd_help.c — keyboard shortcuts dialog (in src_files_gui)
#
# Files in the allowed list are GUI compilation units despite not living under
# src/gui/.  The list must be updated whenever a new GUI-only .c file is placed
# outside src/gui/.
#
# Usage: scripts/lint-gtk-includes.sh [src-root]
#   src-root defaults to "src" relative to the script's parent directory.

set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
SRC_ROOT="${1:-$REPO_ROOT/src}"

# Files outside src/gui/ that are intentionally in src_files_gui (GUI-only).
ALLOWED="
$SRC_ROOT/main.c
$SRC_ROOT/livestacking/gui.c
$SRC_ROOT/core/siril_cmd_help.c
"

bad=0

for pattern in '#include <gtk/gtk.h>' '#include <gdk/gdk.h>'; do
    # Find all matches, then filter out the allowed list and src/gui/
    while IFS=: read -r file line content; do
        # Skip files under src/gui/
        case "$file" in
            "$SRC_ROOT/gui/"*) continue ;;
        esac

        # Skip explicitly allowed files
        allowed=0
        for a in $ALLOWED; do
            if [ "$file" = "$a" ]; then
                allowed=1
                break
            fi
        done
        [ "$allowed" -eq 1 ] && continue

        echo "GTK-include violation: $file:$line: $content"
        bad=1
    done <<EOF
$(grep -rn --include='*.c' --include='*.h' -F "$pattern" "$SRC_ROOT" 2>/dev/null || true)
EOF
done

if [ "$bad" -eq 1 ]; then
    echo ""
    echo "FAIL: GTK/GDK headers found outside GUI compilation units."
    echo "Move the code to src/gui/ or route calls through gui_iface."
    exit 1
fi

echo "OK: no GTK/GDK header violations found."
exit 0
