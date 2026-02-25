/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#include <locale.h>

#include "core/siril.h"
#include "core/proto.h"
#include "gui/message_dialog.h"
#include "core/siril_log.h"
#include "gui/utils.h"


#define GET_DOCUMENTATION_URL "https://siril.readthedocs.io"

void siril_get_documentation(const gchar *page_path) {
	gboolean ret;
	const char *locale;
	const char *supported_languages[] = { "de", "fr", "it", "ru", NULL };  // en is NULL: default language
	gchar *lang = NULL;
	int i = 0;

	if (!com.pref.lang || !g_strcmp0(com.pref.lang, "")) {
		locale = setlocale(LC_MESSAGES, NULL);
	} else {
		locale = com.pref.lang;
	}

	if (locale) {
		while (supported_languages[i]) {
			if (!strncmp(locale, supported_languages[i], 2)) {
				lang = g_strndup(locale, 2);
				break;
			}
			i++;
		}
	}
	if (!lang) {
		lang = g_strdup_printf("en"); // Last gasp fallback in case there is an error with the locale
	}
	const gchar *version = NULL;
#ifdef SIRIL_UNSTABLE
	version = "latest";
#else
	version = "stable";
#endif
	gchar *url = g_build_path (G_DIR_SEPARATOR_S, GET_DOCUMENTATION_URL, lang, version, page_path, NULL);
	control_window_switch_to_tab(OUTPUT_LOGS);
	siril_log_message(_("Siril documentation URL: %s\n"), url);

#if GTK_CHECK_VERSION(3, 22, 0)
	GtkWidget* win = lookup_widget("control_window");
	ret = gtk_show_uri_on_window(GTK_WINDOW(GTK_APPLICATION_WINDOW(win)), url,
			gtk_get_current_event_time(), NULL);
#else
	ret = gtk_show_uri(gdk_screen_get_default(), url,
			gtk_get_current_event_time(), NULL);
#endif
	if (!ret) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Could not show link"),
				_("Please go to <a href=\""GET_DOCUMENTATION_URL"\">"GET_DOCUMENTATION_URL"</a> "
								"by copying the link."));
	}
	g_free(url);
	g_free(lang);
}
