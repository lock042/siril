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
	GtkTextMark *start_mark;
	GtkTextMark *end_mark;
} MatchPosition;


int on_open_pythonpad(GtkMenuItem *menuitem, gpointer user_data);
void on_pythondebug_toggled(GtkCheckMenuItem *item, gpointer user_data);
void set_code_view_theme();
gboolean code_view_exists();
#endif
