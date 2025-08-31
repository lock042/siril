#ifndef GUI_PYTHON
#define GUI_PYTHON
#include <gtksourceview/gtksource.h>

typedef struct {
	GtkTextTag *search_tag;
	GtkSourceView *source_view;
	GtkTextBuffer *buffer;
	GtkEntry *search_entry;
	GtkTextIter current_match_start;
	GtkTextIter current_match_end;
	GtkLabel *info_label;
	gint current_match;
	gint total_matches;
	GSList *match_positions;
} SearchData;

typedef struct {
    GtkWidget *tab_widget;
    GtkLabel *tab_label;
    GtkSourceView *source_view;
    GtkSourceBuffer *source_buffer;  // Changed from GtkTextBuffer
    GtkSourceMap *minimap;
    GFile *file;
    gboolean modified;
    gint language;
    SearchData *search_data;
    gchar *title;
    GtkSourceSpaceDrawer *space_drawer;
} TabInfo;

typedef struct {
	GtkTextMark *start_mark;
	GtkTextMark *end_mark;
} MatchPosition;

// Add new function declarations for tab management
void setup_search_for_tab(TabInfo *tab);
void toggle_find_overlay(gboolean show);
void hide_find_box();

// Existing function declarations
int on_open_pythonpad(GtkMenuItem *menuitem, gpointer user_data);
void on_pythondebug_toggled(GtkCheckMenuItem *item, gpointer user_data);
void set_code_view_theme();
gboolean code_view_exists();
void new_script(const gchar *content, gint length, const char *ext);
gboolean script_editor_has_unsaved_changes();
gboolean get_python_debug_mode();

#endif
