#ifndef SRC_FILTERS_REMIXER_H_
#define SRC_FILTERS_REMIXER_H_

#define CALL_FROM_STARNET 1
#define CALL_FROM_MENU 2

gboolean on_remix_close_clicked(GtkButton *button, gpointer user_data); // callback needed
int toggle_remixer_window_visibility(int _invocation, fits* _fit_left, fits* _fit_right);


#endif /* SRC_FILTERS_REMIXER_H_ */
