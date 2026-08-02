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

/*
 * Language lists that we want to generate only once at program startup:
 * @l10n_lang_list: all available localizations self-localized;
 * @all_lang_list: all known languages, in the user-selected language.
 */

#include <locale.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "gui/utils.h"

#ifdef OS_OSX
#include <CoreFoundation/CoreFoundation.h>
#endif

#include "siril_language.h"

static GHashTable *l10n_lang_list = NULL;
static GHashTable *full_lang_list = NULL;

parsed_code locale_str[] = {
	{"ar_DZ", "العربية"},
	{"ca", "Català"},
	{"da", "Dansk"},
	{"de", "Deutsch"},
	{"el", "Ελληνικά"},
	{"en", "English"},
	{"es_ES", "Espanol"},
	{"fr", "Français"},
	{"it_IT", "Italiano"},
	{"ja_JP", "日本語"},
	{"ko_KR", "한국어"},
	{"nl_BE", "Nederlands"},
	{"pl_PL", "Polish"},
	{"pt_BR", "Português brasileiro"},
	{"pt_PT", "Português"},
	{"ru", "русский"},
	{"tl_PH", "Tagalog"},
	{"uk", "Українська"},
	{"zh_CN", "汉语"},
	{"zh_TW", "正體中文"},
	{NULL, NULL}
};

static GHashTable *parse_locale_codes(GHashTable *table) {
	GHashTableIter lang_iter;
	gpointer key;

	GHashTable *string_lang = g_hash_table_new_full(g_str_hash, g_str_equal,
			(GDestroyNotify) g_free, (GDestroyNotify) g_free);

	g_hash_table_iter_init(&lang_iter, table);
	while (g_hash_table_iter_next(&lang_iter, &key, NULL)) {
		gchar *code = (gchar*) key;
		gchar *str_name = NULL;
		int i = 0;

		while (locale_str[i].locale) {
			if (!g_strcmp0(code, locale_str[i].locale)) {
				str_name = g_strdup(locale_str[i].language_name);
				break;
			}
			i++;
		}
		g_hash_table_insert(string_lang,
				g_strdup_printf("%s [%s]", str_name ? str_name : "???", code),
				NULL);
		g_free(str_name);
	}
	return string_lang;
}

/* extract locale from a string following this pattern:
 * xxxxxxxxx [locale]
 */
static gchar *extract_locale_from_string(const gchar *str) {
	gchar *locale = g_strstr_len(str, -1, "[");
	return g_strndup(locale + 1, strlen(locale) - 2);
}

void siril_language_parser_init() {
	GDir *locales_dir;

	if (l10n_lang_list != NULL) {
		g_warning("siril_language_parser_init() must be run only once.");
		return;
	}

	l10n_lang_list = g_hash_table_new_full(g_str_hash, g_str_equal,
			(GDestroyNotify) g_free, (GDestroyNotify) g_free);

	/* By default Siril is written in english */
	g_hash_table_insert(l10n_lang_list, g_strdup("en"), NULL);
	/* Check all locales we have translations for. */
	locales_dir = g_dir_open(siril_get_locale_dir(), 0, NULL);
	if (locales_dir) {
		const gchar *locale;

		while ((locale = g_dir_read_name(locales_dir)) != NULL) {
			gchar *filename = g_build_filename(siril_get_locale_dir(), locale,
					"LC_MESSAGES",
					PACKAGE ".mo",
					NULL);
			if (g_file_test(filename, G_FILE_TEST_EXISTS)) {
				/* Save the full language code. */
				g_hash_table_insert(l10n_lang_list, g_strdup(locale), NULL);
			}

			g_free(filename);
		}

		g_dir_close(locales_dir);
	}
	full_lang_list = parse_locale_codes(l10n_lang_list);
}

static gint locale_compare(gconstpointer *a, gconstpointer *b) {
	const gchar *s1 = (const gchar *) a;
	const gchar *s2 = (const gchar *) b;

	gchar *k1 = extract_locale_from_string(s1);
	gchar *k2 = extract_locale_from_string(s2);

	gint result = g_strcmp0(k1, k2);
	g_free(k1);
	g_free(k2);

	return result;
}

void siril_language_fill_combo(const gchar *language) {
	GtkComboBoxText *lang_combo = GTK_COMBO_BOX_TEXT(lookup_widget("combo_language"));
	GList *list = g_hash_table_get_keys(full_lang_list);
	gboolean lang_changed = FALSE;
	int i = 1;

	gtk_combo_box_text_remove_all(lang_combo);
	gtk_combo_box_text_append(lang_combo, 0, _("System Language"));

	list = g_list_sort(list, (GCompareFunc) locale_compare);

	for (GList *l = list; l; l = l->next) {
		gtk_combo_box_text_append_text(lang_combo, l->data);
		gchar *locale = extract_locale_from_string(l->data);
		if (!g_strcmp0(language, locale)) {
			gtk_combo_box_set_active(GTK_COMBO_BOX(lang_combo), i);
			lang_changed = TRUE;
		}
		g_free(locale);
		i++;
	}
	if (!lang_changed) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(lang_combo), 0);
	}

	g_list_free(list);
}

#ifdef OS_OSX
/* TRUE if we ship a message catalog for this exact locale code. */
static gboolean have_catalog(const gchar *code) {
	gchar *filename = g_build_filename(siril_get_locale_dir(), code,
			"LC_MESSAGES", PACKAGE ".mo", NULL);
	gboolean exists = g_file_test(filename, G_FILE_TEST_EXISTS);
	g_free(filename);
	return exists;
}

/* Find the catalog that serves @code, which is a locale code as the platform
 * spells it ("nl_NL", "de_CH", "zh_Hant_TW"). An exact match wins; failing
 * that any catalog for the same language will do, because gettext strips the
 * territory but never substitutes a different one: a nl_NL locale would find
 * GTK's nl catalog but not our nl_BE one, leaving the interface half
 * translated. Returns the catalog name, or NULL if we have none.
 */
static gchar *find_catalog_for(const gchar *code) {
	if (have_catalog(code))
		return g_strdup(code);

	gchar **parts = g_strsplit(code, "_", -1);
	const gchar *lang = parts[0];
	const gchar *region = NULL;	// last subtag, if any: "zh_Hant_TW" -> "TW"
	for (guint i = 1; parts[i]; i++)
		region = parts[i];

	gchar *found = NULL;
	if (have_catalog(lang)) {
		found = g_strdup(lang);		// "de" serves "de_CH"
	} else {
		GSList *candidates = NULL;
		GDir *locales_dir = g_dir_open(siril_get_locale_dir(), 0, NULL);
		if (locales_dir) {
			const gchar *name;
			gsize len = strlen(lang);

			while ((name = g_dir_read_name(locales_dir)) != NULL) {
				if (strncmp(name, lang, len) || (name[len] && name[len] != '_'))
					continue;
				if (have_catalog(name))	// sorted, g_dir_read_name order is not
					candidates = g_slist_insert_sorted(candidates,
							g_strdup(name), (GCompareFunc) g_strcmp0);
			}
			g_dir_close(locales_dir);
		}
		/* prefer the catalog for the region that was actually asked for, so
		 * that zh_Hant_TW picks zh_TW rather than zh_CN */
		for (GSList *l = candidates; l && region; l = l->next) {
			const gchar *sep = strchr((const gchar *) l->data, '_');
			if (sep && !g_ascii_strcasecmp(sep + 1, region)) {
				found = g_strdup((const gchar *) l->data);
				break;
			}
		}
		if (!found && candidates)
			found = g_strdup((const gchar *) candidates->data);
		g_slist_free_full(candidates, g_free);
	}
	g_strfreev(parts);
	return found;
}

/* The user's preferred interface languages, most wanted first, as BCP 47 tags
 * rewritten into the POSIX form gettext uses ("nl-NL" -> "nl_NL"). This is the
 * raw preference list, not the one macOS filters through the bundle's
 * CFBundleLocalizations, so it is unaffected by which localizations the .app
 * happens to advertise. Free with g_strfreev(). */
static gchar **macos_preferred_languages(void) {
	CFPropertyListRef value = CFPreferencesCopyAppValue(CFSTR("AppleLanguages"),
			kCFPreferencesCurrentApplication);
	if (!value)
		return NULL;
	if (CFGetTypeID(value) != CFArrayGetTypeID()) {
		CFRelease(value);
		return NULL;
	}

	CFArrayRef languages = (CFArrayRef) value;
	CFIndex count = CFArrayGetCount(languages);
	GPtrArray *codes = g_ptr_array_new();

	for (CFIndex i = 0; i < count; i++) {
		CFTypeRef tag = CFArrayGetValueAtIndex(languages, i);
		char buffer[64];

		if (!tag || CFGetTypeID(tag) != CFStringGetTypeID())
			continue;
		if (!CFStringGetCString((CFStringRef) tag, buffer, sizeof(buffer),
				kCFStringEncodingASCII))
			continue;
		for (char *p = buffer; *p; p++)
			if (*p == '-')
				*p = '_';
		g_ptr_array_add(codes, g_strdup(buffer));
	}
	CFRelease(value);

	g_ptr_array_add(codes, NULL);
	return (gchar **) g_ptr_array_free(codes, FALSE);
}

/* Pick a language for the "System Language" setting.
 *
 * A bundled application launched from the Finder gets no LANG or LC_* in its
 * environment at all, so gettext falls back to CoreFoundation and uses
 * AppleLocale, which is the first preferred language glued to the region,
 * regardless of what the application can actually offer. English is not a
 * catalog, it is only what a failed lookup falls back to, so it can never win
 * that comparison: an English-first user whose second language is Dutch gets
 * our strings in English (we have no nl_NL catalog) and GTK's in Dutch (it has
 * an nl one). Do the negotiation ourselves instead, walking the preferences in
 * order and stopping at the first language we can serve, counting English as
 * always available since it is the language the sources are written in.
 *
 * Returns the language to use, or NULL to leave the environment alone.
 */
static gchar *negotiate_system_language(void) {
	static const gchar * const locale_vars[] = { "LANGUAGE", "LC_ALL", "LC_MESSAGES", "LANG" };

	for (gsize i = 0; i < G_N_ELEMENTS(locale_vars); i++) {
		const gchar *value = g_getenv(locale_vars[i]);
		if (value && value[0] != '\0')
			return NULL;	// started from a shell that asked for something
	}

	gchar **preferred = macos_preferred_languages();
	if (!preferred)
		return NULL;

	gchar *language = NULL;
	for (gsize i = 0; preferred[i] && !language; i++) {
		if (!g_ascii_strncasecmp(preferred[i], "en", 2)
				&& (preferred[i][2] == '\0' || preferred[i][2] == '_'))
			language = g_strdup("en");
		else
			language = find_catalog_for(preferred[i]);
	}
	g_strfreev(preferred);

	/* nothing we can serve: English is every application's last resort here,
	 * and picking it keeps GTK's half of the interface in step with ours */
	return language ? language : g_strdup("en");
}
#endif

void language_init(const gchar *language) {
	gchar *negotiated = NULL;

	if ((!language) || (language[0] == '\0')) {
#ifdef OS_OSX
		negotiated = negotiate_system_language();
		if (!negotiated)
			return;
		language = negotiated;
#else
		return;
#endif
	}

	g_mutex_lock(&com.env_mutex);
	/* This is default language */
	if (!g_ascii_strcasecmp(language, "en")) {
		if (!g_setenv("LANGUAGE", "C", TRUE))
			siril_debug_print("Error setting LANGUAGE to C\n");
	} else {
		if (!g_setenv("LANGUAGE", language, TRUE))
			siril_debug_print("Error setting LANGUAGE\n");
	}
	g_mutex_unlock(&com.env_mutex);
	setlocale(LC_ALL, "");
	setlocale(LC_NUMERIC, "C");
	if (negotiated) {
		/* worth a line in the log: this is the case users report against */
		siril_log_message(_("Interface language taken from the system preferences: %s\n"),
				negotiated);
		g_free(negotiated);
	}
}

gchar *get_interface_language() {
	GtkComboBoxText *lang_combo = GTK_COMBO_BOX_TEXT(lookup_widget("combo_language"));

	if (gtk_combo_box_get_active(GTK_COMBO_BOX(lang_combo)) == 0) {
		return g_strdup("");
	} else {
		gchar *str = gtk_combo_box_text_get_active_text(lang_combo);
		return extract_locale_from_string(str);
	}
}
