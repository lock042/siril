#!/bin/sh
# lint-gtk-includes.sh — enforce the GUI boundary in non-GUI source files.
#
# Siril has two parallel GUI trees:
#   src/gui/       — the GTK3 build
#   src/gui-gtk4/  — the GTK4 build (one of these is compiled per build,
#                    selected via the `gtk_version` meson option)
# Both trees are GUI compilation units and may include GTK/GDK headers.
#
# CHECK 1: Direct GTK/GDK header includes
#   Fails if any file outside the GUI compilation units includes <gtk/gtk.h>
#   or <gdk/gdk.h> directly.
#
# CHECK 2: gui*/ header includes from outside the GUI trees
#   Fails if any non-GUI source file includes a "gui/..." or "gui-gtk4/..."
#   header, unless the (file, header) pair is in the whitelist below.
#
#   Rationale: headers under src/gui*/ are part of a GUI compilation unit
#   and may pull in GTK types even when they are themselves GTK-free.
#   Keeping them out of non-GUI code prevents future GTK leakage and
#   documents which cross-boundary dependencies still exist.
#
#   Whitelist entries must include a justification comment explaining why
#   the include is acceptable and what would be needed to remove it.
#
# Usage: scripts/lint-gtk-includes.sh [src-root]
#   src-root defaults to "src" relative to the script's parent directory.

set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
SRC_ROOT="${1:-$REPO_ROOT/src}"

# ── Check 1: direct GTK/GDK headers ────────────────────────────────────────

# Files outside src/gui/ that are intentionally GUI-only compilation units.
GTK_ALLOWED="
$SRC_ROOT/main.c
"

bad=0

for pattern in '#include <gtk/gtk.h>' '#include <gdk/gdk.h>'; do
    while IFS=: read -r file line content; do
        case "$file" in
            "$SRC_ROOT/gui/"*|"$SRC_ROOT/gui-gtk4/"*) continue ;;
        esac

        allowed=0
        for a in $GTK_ALLOWED; do
            [ "$file" = "$a" ] && allowed=1 && break
        done
        [ "$allowed" -eq 1 ] && continue

        echo "GTK-include violation: $file:$line: $content"
        bad=1
    done <<EOF
$(grep -rn --include='*.c' --include='*.h' --include='*.cpp' --include='*.hpp' -F "$pattern" "$SRC_ROOT" 2>/dev/null || true)
EOF
done

# ── Check 2: gui/ header includes from outside src/gui/ ────────────────────

# Whitelist of (file_suffix:gui_header) pairs that are intentionally allowed.
# Format: "FILE_SUFFIX:GUI_HEADER"
#   FILE_SUFFIX — path suffix after $SRC_ROOT/ (used for substring match)
#   GUI_HEADER  — the included header path after "gui/"
#
# To add an entry: document WHY it is acceptable and what clean-up would
# look like, then add "file_suffix:header" to the list.
GUI_INCLUDE_WHITELIST="
io/siril_plot.h:plot.h
io/siril_plot.c:plot.h
io/siril_pythoncommands.c:user_polygons.h
io/siril_pythonmodule.c:user_polygons.h
"
# Justifications:
#   io/siril_plot.{h,c} : plot.h — siril_plot_data embeds plot_draw_data_t by
#     value; moving that type out of gui/plot.h requires splitting the kplot
#     rendering API.  gui/plot.h is Cairo-only (no GTK).
#   io/siril_python*.c  : user_polygons.h — the Python IPC bridge manages
#     polygon overlay data (pure data structures, no GTK).  The header and
#     its implementation should eventually move out of gui/ to algos/ or io/.

while IFS=: read -r file line content; do
    # Skip files under src/gui/ and src/gui-gtk4/ (both are GUI trees)
    case "$file" in
        "$SRC_ROOT/gui/"*|"$SRC_ROOT/gui-gtk4/"*) continue ;;
    esac
    # Skip main.c and main-cli.c (entry points)
    case "$file" in
        "$SRC_ROOT/main.c"|"$SRC_ROOT/main-cli.c") continue ;;
    esac
    # Skip test files
    case "$file" in
        "$SRC_ROOT/tests/"*) continue ;;
    esac

    # Extract the included gui*/ header name (everything after "gui<…>/").
    gui_header=$(printf '%s' "$content" | sed -E 's|.*"gui(-gtk4)?/([^"]*)".*|\2|')

    # Build the file suffix (path relative to SRC_ROOT)
    file_suffix="${file#$SRC_ROOT/}"

    # Check whitelist
    whitelisted=0
    for entry in $GUI_INCLUDE_WHITELIST; do
        wl_file="${entry%%:*}"
        wl_hdr="${entry##*:}"
        if [ "$file_suffix" = "$wl_file" ] && [ "$gui_header" = "$wl_hdr" ]; then
            whitelisted=1
            break
        fi
    done
    [ "$whitelisted" -eq 1 ] && continue

    echo "gui-include violation: $file:$line: $content"
    bad=1
done <<EOF
$(grep -rn --include='*.c' --include='*.h' --include='*.cpp' --include='*.hpp' \
    -E '#include "gui(-gtk4)?/[^"]+\.h"' "$SRC_ROOT" 2>/dev/null || true)
EOF

if [ "$bad" -eq 1 ]; then
    echo ""
    echo "FAIL: GUI boundary violations found."
    echo "Move the code to src/gui/ (GTK3) or src/gui-gtk4/ (GTK4),"
    echo "route calls through gui_iface, or add an entry to the whitelist"
    echo "in scripts/lint-gtk-includes.sh with a justification comment."
    exit 1
fi

echo "OK: no GUI boundary violations found."
exit 0
