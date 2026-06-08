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

#include <glib.h>
#include <gtk/gtk.h>
#include "algos/astrometry_solver.h"
#include "algos/siril_wcs.h"
#include "registration/registration.h"
#include "io/annotation_catalogues.h"
#include "algos/search_objects.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/callbacks.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/image_display.h"
#include "gui-gtk4/image_interactions.h"
#include "gui-gtk4/PSF_list.h"
#include "gui-gtk4/photometric_cc.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/siril_catalogues.h"
#include "io/local_catalogues.h"

/* forward-declared here after being removed from algos/astrometry_solver.h (plan 1.4) */
void on_GtkButton_IPS_metadata_clicked(GtkButton *button, gpointer user_data);

enum {
	COLUMN_RESOLVER,// string
	COLUMN_NAME,	// string
	N_COLUMNS
};

/* ---- Per-row GObject for the IPS resolver list ------------------- */
/* Forward-declared statics so the factory callbacks below can reference
 * the live model and selection. */
static GListStore        *list_IPS;
static GtkSingleSelection *selection_model;
static GtkColumnView     *columnview_IPS;
static GtkScrolledWindow *scrolled_IPS;
static void on_IPS_selection_changed(GObject *sm, GParamSpec *pspec, gpointer user_data);

#define SIRIL_TYPE_IPS_ROW (siril_ips_row_get_type())
G_DECLARE_FINAL_TYPE(SirilIPSRow, siril_ips_row, SIRIL, IPS_ROW, GObject)
struct _SirilIPSRow {
	GObject parent_instance;
	gchar *resolver;
	gchar *name;
};
G_DEFINE_TYPE(SirilIPSRow, siril_ips_row, G_TYPE_OBJECT)
static void siril_ips_row_finalize(GObject *obj) {
	SirilIPSRow *self = SIRIL_IPS_ROW(obj);
	g_clear_pointer(&self->resolver, g_free);
	g_clear_pointer(&self->name,     g_free);
	G_OBJECT_CLASS(siril_ips_row_parent_class)->finalize(obj);
}
static void siril_ips_row_class_init(SirilIPSRowClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_ips_row_finalize;
}
static void siril_ips_row_init(SirilIPSRow *self) { (void)self; }
static SirilIPSRow *siril_ips_row_new(const gchar *resolver, const gchar *name) {
	SirilIPSRow *self = g_object_new(SIRIL_TYPE_IPS_ROW, NULL);
	self->resolver = g_strdup(resolver);
	self->name     = g_strdup(name);
	return self;
}

static void ips_setup_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkWidget *lbl = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(lbl), 0.0);
	gtk_list_item_set_child(li, lbl);
}
static void ips_bind_resolver_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilIPSRow *row = SIRIL_IPS_ROW(gtk_list_item_get_item(li));
	gtk_label_set_text(lbl, row->resolver);
}
static void ips_bind_name_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	(void)f; (void)u;
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilIPSRow *row = SIRIL_IPS_ROW(gtk_list_item_get_item(li));
	gtk_label_set_text(lbl, row->name);
}

/* Build the IPS column view if not already built; idempotent. */
static void ensure_ips_columnview(void) {
	if (columnview_IPS) return;
	if (!scrolled_IPS) return;
	if (!list_IPS)
		list_IPS = g_list_store_new(SIRIL_TYPE_IPS_ROW);
	selection_model = gtk_single_selection_new(G_LIST_MODEL(g_object_ref(list_IPS)));
	gtk_single_selection_set_can_unselect(selection_model, TRUE);
	gtk_single_selection_set_autoselect(selection_model, FALSE);
	g_signal_connect(selection_model, "notify::selected", G_CALLBACK(on_IPS_selection_changed), NULL);

	columnview_IPS = GTK_COLUMN_VIEW(gtk_column_view_new(GTK_SELECTION_MODEL(g_object_ref(selection_model))));
	gtk_widget_add_css_class(GTK_WIDGET(columnview_IPS), "siril-dense-rows");

	GtkSignalListItemFactory *fr = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(fr, "setup", G_CALLBACK(ips_setup_cb),         NULL);
	g_signal_connect(fr, "bind",  G_CALLBACK(ips_bind_resolver_cb), NULL);
	GtkColumnViewColumn *cr = gtk_column_view_column_new(N_("Resolver"), GTK_LIST_ITEM_FACTORY(fr));
	gtk_column_view_append_column(columnview_IPS, cr);
	g_object_unref(cr);

	GtkSignalListItemFactory *fn = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(fn, "setup", G_CALLBACK(ips_setup_cb),     NULL);
	g_signal_connect(fn, "bind",  G_CALLBACK(ips_bind_name_cb), NULL);
	GtkColumnViewColumn *cn = gtk_column_view_column_new(N_("Name"), GTK_LIST_ITEM_FACTORY(fn));
	gtk_column_view_column_set_expand(cn, TRUE);
	gtk_column_view_append_column(columnview_IPS, cn);
	g_object_unref(cn);

	gtk_scrolled_window_set_child(scrolled_IPS, GTK_WIDGET(columnview_IPS));
}

// caching all UI elements
static GtkBox *IPSbox_seq_info = NULL;
static GtkCheckButton *flipbutton = NULL, *automagbutton = NULL, *DEC_S = NULL,
	*manualbutton = NULL, *downsamplebutton = NULL, *autocropbutton = NULL, *autocatbutton = NULL,
	*nonearbutton = NULL, *blindposbutton = NULL, *blindresbutton = NULL,
	*seqsolvebutton = NULL, *seqnocache = NULL, *seqskipsolved = NULL,
	*sequseheadercoords = NULL, *sequseheaderpixel = NULL, *sequseheaderfocal = NULL,
	*sequseforreg = NULL, *masterbutton = NULL;
static GtkButton *distomaster_save_button = NULL;
static GtkSpinButton *magspin = NULL, *RA_h = NULL, *RA_m = NULL, *DEC_d = NULL, *DEC_m = NULL, *radiusspin = NULL;
static GtkDropDown *catalogbox = NULL, *orderbox = NULL, *solverbox = NULL, *serverbox = NULL;
static GtkLabel *cataloglabel = NULL, *radiuslabel = NULL;
static GtkEntry *focalentry = NULL, *pixelentry = NULL, *resolutionentry = NULL,
	*RA_s = NULL, *DEC_s = NULL, *distomaster_entry = NULL;
/* GTK4: GtkSearchEntry no longer derives from GtkEntry — both implement
 * GtkEditable.  Hold this one as GtkEditable so the editable API works. */
static GtkEditable *searchentry = NULL;
/* Phase 11: GTK4 model + view, replacing GtkListStore + GtkTreeView.
 * Definitions for the forward-declared statics above. */
static GtkExpander *cataloguesexp = NULL, *stardetectionexp = NULL, *sequenceexp = NULL;
static GtkWindow *astrometry_dialog = NULL;
static gboolean have_local_cat = FALSE, radius_set = FALSE, order_set = FALSE, have_asnet = FALSE,
				has_coords = FALSE, has_pixel = FALSE, has_focal = FALSE; // those bools tell if the metadata was present in the header of gfit
static gboolean use_local_catalogue();
void on_comboastro_catalog_changed(GObject *obj, GParamSpec *pspec, gpointer user_data);
void on_comboastro_solver_changed(GObject *obj, GParamSpec *pspec, gpointer user_data);
void on_comboastro_order_changed(GObject *obj, GParamSpec *pspec, gpointer user_data);
void on_GtkCheckButton_solveseq_toggled(GtkCheckButton *button, gpointer user);
static int get_order();
extern struct sky_object platedObject[RESOLVER_NUMBER];

static void unselect_all_items();
static void on_IPS_selection_changed(GObject *sm, GParamSpec *pspec, gpointer user_data);

void reset_astrometry_checks() {
	have_local_cat = local_catalogues_available();
	have_asnet = asnet_is_available();
	radius_set = FALSE;
	order_set= FALSE;
	if (!have_asnet) {
		gtk_drop_down_set_selected(solverbox, 0);
		gtk_widget_set_sensitive(GTK_WIDGET(solverbox), FALSE);
	}
}
// this one is called by init_astrometry()
// which is called first when loading the app
static void load_all_ips_statics() {
	if (!flipbutton) {
		// toggles
		flipbutton = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkButton_IPS_flip"));
		automagbutton = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "GtkCheckButton_Mag_Limit"));
		DEC_S = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "GtkCheckButtonIPS_S"));
		manualbutton = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkButton_IPS_manual"));
		downsamplebutton = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "downsample_ips_button"));
		autocropbutton = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "autocrop_ips_button"));
		autocatbutton = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "GtkCheckButton_OnlineCat"));
		nonearbutton = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_IPS_nonear"));
		blindposbutton = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_IPS_ignorepos"));
		blindresbutton = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_IPS_ignoreres"));
		seqsolvebutton = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_IPS_seqsolve"));
		seqnocache = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_IPS_nocache"));
		seqskipsolved = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_IPS_skipsolved"));
		sequseheadercoords = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_IPS_sequseheadercoords"));
		sequseheaderpixel = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_IPS_sequseheaderpixel"));
		sequseheaderfocal = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_IPS_sequseheaderfocal"));
		sequseforreg = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutton_IPS_useforreg"));
		masterbutton = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "master_ips_button"));
		distomaster_save_button = GTK_BUTTON(gtk_builder_get_object(gui.builder, "distomaster_save_button"));
		// combos
		catalogbox = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "ComboBoxIPSCatalog"));
		orderbox = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "ComboBoxIPS_order"));
		solverbox = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "ComboBoxIPS_solver"));
		serverbox = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combo_server_ips"));
		// labels
		cataloglabel = GTK_LABEL(gtk_builder_get_object(gui.builder, "astrometry_catalog_label"));
		radiuslabel = GTK_LABEL(gtk_builder_get_object(gui.builder, "label_IPS_radius"));
		// spins
		magspin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "GtkSpinIPS_Mag_Limit"));
		radiusspin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "GtkSpinIPS_radius"));
		RA_h = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "GtkSpinIPS_RA_h"));
		RA_m = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "GtkSpinIPS_RA_m"));
		DEC_d = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "GtkSpinIPS_Dec_deg"));
		DEC_m = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "GtkSpinIPS_Dec_m"));
		// entries
		focalentry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "GtkEntry_IPS_focal"));
		pixelentry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "GtkEntry_IPS_pixels"));
		resolutionentry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "GtkEntry_IPS_resolution"));
		RA_s = GTK_ENTRY(gtk_builder_get_object(gui.builder, "GtkEntryIPS_RA_s"));
		DEC_s = GTK_ENTRY(gtk_builder_get_object(gui.builder, "GtkEntryIPS_Dec_s"));
		searchentry = GTK_EDITABLE(gtk_builder_get_object(gui.builder, "GtkSearchIPS"));
		distomaster_entry = GTK_ENTRY(gtk_builder_get_object(gui.builder, "distomaster_entry"));
		// scrolled window holding the column view (Phase 11)
		scrolled_IPS = GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "GtkScrollIPS"));
		ensure_ips_columnview();
		// expanders
		cataloguesexp = GTK_EXPANDER(gtk_builder_get_object(gui.builder, "labelIPSCatalogParameters"));
		stardetectionexp = GTK_EXPANDER(gtk_builder_get_object(gui.builder, "Frame_IPS_star_detection"));
		sequenceexp = GTK_EXPANDER(gtk_builder_get_object(gui.builder, "Frame_IPS_sequence"));
		// box
		IPSbox_seq_info = GTK_BOX(gtk_builder_get_object(gui.builder, "IPSbox_seq_info"));
		// window
		astrometry_dialog = GTK_WINDOW(gtk_builder_get_object(gui.builder, "astrometry_dialog"));
	}
}

void initialize_ips_dialog() {
	if (!order_set) {
		gtk_drop_down_set_selected(orderbox, com.pref.astrometry.sip_correction_order - 1);
		order_set = TRUE;
	}
	if (!radius_set) {
		gtk_spin_button_set_value(radiusspin, com.pref.astrometry.radius_degrees);
		radius_set = TRUE;
	}
	// loads image metadata
	on_GtkButton_IPS_metadata_clicked(NULL, NULL);	// fill it automatically
	// sequence related controls
	gboolean isseq = sequence_is_loaded() && com.seq.current != RESULT_IMAGE;
	gboolean is_bayer = !isseq && gfit->keywords.bayer_pattern[0] != '\0';
	gtk_widget_set_visible(GTK_WIDGET(flipbutton), !isseq && !is_bayer);
	gtk_expander_set_expanded(sequenceexp, isseq);
	gtk_widget_set_visible(GTK_WIDGET(sequenceexp), isseq);
	gtk_widget_set_visible(GTK_WIDGET(stardetectionexp), !isseq);
	gtk_widget_set_visible(GTK_WIDGET(seqsolvebutton), isseq);
	on_GtkCheckButton_solveseq_toggled(NULL, NULL);
	siril_toggle_set_active(GTK_WIDGET(seqsolvebutton), FALSE);
	if (isseq) {
		siril_toggle_set_active(GTK_WIDGET(sequseheadercoords), has_coords);
		siril_toggle_set_active(GTK_WIDGET(sequseheaderpixel), has_pixel);
		siril_toggle_set_active(GTK_WIDGET(sequseheaderfocal), has_focal);
		gtk_widget_set_sensitive(GTK_WIDGET(sequseheadercoords), has_coords);
		gtk_widget_set_sensitive(GTK_WIDGET(sequseheaderpixel), has_pixel);
		gtk_widget_set_sensitive(GTK_WIDGET(sequseheaderfocal), has_focal);
	}
	siril_toggle_set_active(GTK_WIDGET(nonearbutton), isseq);
	on_comboastro_order_changed(NULL, NULL, NULL);
	// solver-related controls
	on_comboastro_catalog_changed(NULL, NULL, NULL);
}

static void add_style_class(GtkWidget *widget, const char *class) {
	gtk_widget_add_css_class(widget, class);
}

static void remove_style_class(GtkWidget *widget, const char *class) {
	gtk_widget_remove_css_class(widget, class);
}

static void change_coords_colors_to_unset() {
	add_style_class(GTK_WIDGET(RA_h), "header_unset");
	add_style_class(GTK_WIDGET(RA_m), "header_unset");
	add_style_class(GTK_WIDGET(RA_s), "header_unset");
	add_style_class(GTK_WIDGET(DEC_d), "header_unset");
	add_style_class(GTK_WIDGET(DEC_m), "header_unset");
	add_style_class(GTK_WIDGET(DEC_s), "header_unset");
}

static void change_coords_colors_to_set() {
	remove_style_class(GTK_WIDGET(RA_h), "header_unset");
	remove_style_class(GTK_WIDGET(RA_m), "header_unset");
	remove_style_class(GTK_WIDGET(RA_s), "header_unset");
	remove_style_class(GTK_WIDGET(DEC_d), "header_unset");
	remove_style_class(GTK_WIDGET(DEC_m), "header_unset");
	remove_style_class(GTK_WIDGET(DEC_s), "header_unset");
}

static void change_entry_colors_to_unset() {
	change_coords_colors_to_unset();
	add_style_class(GTK_WIDGET(focalentry), "header_unset");
	add_style_class(GTK_WIDGET(pixelentry), "header_unset");
}

static void change_entry_colors_to_set() {
	change_coords_colors_to_set();
	remove_style_class(GTK_WIDGET(focalentry), "header_unset");
	remove_style_class(GTK_WIDGET(pixelentry), "header_unset");
}

static gboolean use_local_catalogue() {
	int cat = gtk_drop_down_get_selected(catalogbox);
	gboolean autocat = siril_toggle_get_active(GTK_WIDGET(autocatbutton));
	return have_local_cat && (autocat || (cat != CAT_PPMXL && cat != CAT_APASS));
}

static void get_mag_settings_from_GUI(limit_mag_mode *mag_mode, double *magnitude_arg) {
	gboolean autob = siril_toggle_get_active(GTK_WIDGET(automagbutton));
	if (autob)
		*mag_mode = LIMIT_MAG_AUTO;
	else {
		*magnitude_arg = gtk_spin_button_get_value(magspin);
		*mag_mode = LIMIT_MAG_ABSOLUTE;
	}
}

gboolean has_any_keywords() {
	return (gfit->keywords.focal_length > 0.0 ||
			gfit->keywords.pixel_size_x > 0.f ||
			gfit->keywords.pixel_size_y > 0.f ||
			(has_wcs(gfit) && gfit->keywords.wcslib->crval[0] != 0.0 && gfit->keywords.wcslib->crval[1] != 0.0) ||
			(gfit->keywords.wcsdata.objctra[0] != '\0' && gfit->keywords.wcsdata.objctdec[0] != '\0') ||
			(gfit->keywords.wcsdata.ra > DEFAULT_DOUBLE_VALUE && gfit->keywords.wcsdata.dec > DEFAULT_DOUBLE_VALUE));
}

/* effective focal length in mm */
static double get_focal() {
	const gchar *value = gtk_editable_get_text(GTK_EDITABLE(focalentry));
	return g_ascii_strtod(value, NULL);	// 0 is parse error
}

/* get pixel in µm */
static double get_pixel() {
	const gchar *value = gtk_editable_get_text(GTK_EDITABLE(pixelentry));
	return g_ascii_strtod(value, NULL);	// 0 is parse error
}

static int get_order() {
	return gtk_drop_down_get_selected(orderbox) + 1;
}

static int get_server_from_combobox() {
	return gtk_drop_down_get_selected(serverbox);
}

static SirilWorldCS *get_center_of_catalog() {
	/* get alpha center */
	gdouble hour = gtk_spin_button_get_value_as_int(RA_h);
	gdouble min = gtk_spin_button_get_value_as_int(RA_m);
	gdouble sec = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(RA_s)), NULL);
	/* get Dec center */
	gdouble deg = gtk_spin_button_get_value_as_int(DEC_d);
	gdouble m = gtk_spin_button_get_value_as_int(DEC_m);
	gdouble s = g_ascii_strtod(gtk_editable_get_text(GTK_EDITABLE(DEC_s)), NULL);
	if (siril_toggle_get_active(GTK_WIDGET(DEC_S))) {
		deg = -deg;
	}
	return siril_world_cs_new_from_ra_dec(hour, min, sec, deg, m, s);;
}

static gboolean is_detection_manual() {
	return siril_toggle_get_active(GTK_WIDGET(manualbutton));
}

static gboolean flip_image_after_ps() {
	return gtk_widget_get_visible(GTK_WIDGET(flipbutton)) && siril_toggle_get_active(GTK_WIDGET(flipbutton));
}

static gboolean is_downsample_activated() {
	return siril_toggle_get_active(GTK_WIDGET(downsamplebutton));
}

static gboolean is_autocrop_activated() {
	return siril_toggle_get_active(GTK_WIDGET(autocropbutton));
}

static gboolean is_save_disto_activated() {
	return gtk_widget_get_sensitive(GTK_WIDGET(masterbutton)) && siril_toggle_get_active(GTK_WIDGET(masterbutton)) && strlen(gtk_editable_get_text(GTK_EDITABLE(distomaster_entry))) > 0;
}

static void update_pixel_size() {
	double pixel = gfit->keywords.pixel_size_x > gfit->keywords.pixel_size_y ? gfit->keywords.pixel_size_x : gfit->keywords.pixel_size_y;
	if (com.pref.binning_update && gfit->keywords.binning_x > 1) {
		pixel *= gfit->keywords.binning_x;
	}

	if (pixel > 0.0) {
		gchar *cpixels = g_strdup_printf("%.2lf", pixel);
		gtk_editable_set_text(GTK_EDITABLE(pixelentry), cpixels);
		g_free(cpixels);
		has_pixel = gfit->pixelkey;
	} else
		has_pixel = FALSE;
	if (!has_pixel)
		add_style_class(GTK_WIDGET(pixelentry), "header_unset");
	else
		remove_style_class(GTK_WIDGET(pixelentry), "header_unset");
}

static void update_focal() {
	double focal = gfit->keywords.focal_length;

	if (focal > 0.0) {
		gchar *cfocal = g_strdup_printf("%.1lf", focal);
		gtk_editable_set_text(GTK_EDITABLE(focalentry), cfocal);
		g_free(cfocal);
		has_focal = gfit->focalkey;
	} else
		has_focal = FALSE;
	if (!has_focal)
		add_style_class(GTK_WIDGET(focalentry), "header_unset");
	else
		remove_style_class(GTK_WIDGET(focalentry), "header_unset");
}

static void update_resolution_field() {
	double res = get_resolution(get_focal(), get_pixel());
	gchar *cres;

	cres = g_strdup_printf("%1.3lf", res);
	gtk_editable_set_text(GTK_EDITABLE(resolutionentry), cres);
	g_free(cres);
}

static void update_coordinates(SirilWorldCS *world_cs) {
	gchar *RA_sec, *Dec_sec;
	gint ra_h, ra_m;
	gint dec_deg, dec_m;
	gdouble ra_s, dec_s;
	if (world_cs) {
		siril_world_cs_get_ra_hour_min_sec(world_cs, &ra_h, &ra_m, &ra_s);
		siril_world_cs_get_dec_deg_min_sec(world_cs, &dec_deg, &dec_m, &dec_s);
	} else {
		ra_h = 0; ra_m = 0; ra_s = 0.0;
		dec_deg = 0; dec_m = 0; dec_s = 0.0;
	}

	RA_sec = g_strdup_printf("%6.4lf", ra_s);
	Dec_sec = g_strdup_printf("%6.4lf", dec_s);

	siril_toggle_set_active(GTK_WIDGET(DEC_S), dec_deg < 0);

	gtk_spin_button_set_value(RA_h, ra_h);
	gtk_spin_button_set_value(RA_m, ra_m);
	gtk_editable_set_text(GTK_EDITABLE(RA_s), RA_sec);

	gtk_spin_button_set_value(DEC_d, abs(dec_deg));
	gtk_spin_button_set_value(DEC_m, dec_m);
	gtk_editable_set_text(GTK_EDITABLE(DEC_s), Dec_sec);

	g_free(RA_sec);
	g_free(Dec_sec);
	if (world_cs)
		change_coords_colors_to_set();
	else
		change_coords_colors_to_unset();
}

void update_coords() {
	SirilWorldCS *world_cs = get_eqs_from_header(gfit);
	update_coordinates(world_cs);
	if (world_cs) {
		unselect_all_items();
		siril_world_cs_unref(world_cs);
		has_coords = TRUE;
	} else
		has_coords = FALSE;
}

static void update_image_parameters_GUI() {
	/* update all fields. Resolution is updated as well thanks to the Entry
	 * and combo changed signal */
	update_focal();
	update_pixel_size();
	update_coords();
}

gboolean end_process_catsearch(gpointer p) {
	sky_object_query_args *args = (sky_object_query_args*)p;
	if (!com.script && !args->retval) {
		GtkToggleButton *button = GTK_TOGGLE_BUTTON(gtk_builder_get_object(gui.builder, "annotate_button"));
		// purge_user_catalogue(CAT_AN_USER_TEMP);
		refresh_annotation_visibility();
		if (!siril_toggle_get_active(GTK_WIDGET(button))) {
			siril_toggle_set_active(GTK_WIDGET(button), TRUE);
		} else {
			refresh_found_objects();
			redraw(REDRAW_OVERLAY);
		}
	}
	set_cursor_waiting(FALSE);
	return(FALSE);
}

gboolean end_plate_solver(gpointer p) {
	struct astrometry_data *args = (struct astrometry_data *) p;
	stop_processing_thread();

	set_cursor_waiting(FALSE);

	if (args->ret > 0) {
		char *title = siril_log_error(_("Plate Solving failed. "
					"The image could not be aligned with the reference stars.\n"));
		gchar *msg = NULL;
		if (args->ret == SOLVE_NO_MATCH)
			msg = g_strdup(_("This is usually because the initial parameters (pixel size, focal length, initial coordinates) "
					"are too far from the real metadata of the image.\n\n"
					"You could also try to look into another catalogue, or try to click on the \"Downsampling\" button, especially for image done with Drizzle.\n\n"
					"Finally, keep in mind that plate solving algorithm should only be applied on linear image."));
		else
			msg = platesolve_msg(args);
		siril_message_dialog(GTK_MESSAGE_ERROR, title, msg);
		g_free(msg);
	} else {
		if (args->ret) {
			gchar *msg = platesolve_msg(args);
			siril_log_warning(msg);
			g_free(msg);
		}
		/* update UI */
		update_image_parameters_GUI();
		set_GUI_CAMERA();
		update_coords();
		delete_selected_area();
		refresh_annotations(FALSE);
		close_astrometry_dialog();
		/* ****** */
		if (args->image_flipped) {
			redraw(REDRAW_ALL);
		} else {
			redraw(REDRAW_OVERLAY);
		}
	}
	if (args->image_flipped)
		clear_stars_list(TRUE);
	gui_function(update_MenuItem, NULL);
	free_astrometry_data(args);
	return FALSE;
}

static void start_image_plate_solve() {
	struct astrometry_data *args = calloc(1, sizeof(struct astrometry_data));
	set_cursor_waiting(TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);
	if (!fill_plate_solver_structure_from_GUI(args)) {
		if (!args->for_sequence) {
			if (!start_in_new_thread(plate_solver, args)) {
				free_astrometry_data(args);
			}
		} else {
			start_sequence_astrometry(&com.seq, args);
		}
	} else {
		free_astrometry_data(args);
		set_cursor_waiting(FALSE);
	}
}

static void clear_all_objects() {
	if (list_IPS) g_list_store_remove_all(list_IPS);
}

static void append_ips_row(const gchar *resolver, const gchar *name) {
	SirilIPSRow *row = siril_ips_row_new(resolver, name);
	g_list_store_append(list_IPS, row);
	g_object_unref(row);
}

static void add_object_to_list() {
	gboolean anyobject = FALSE;

	clear_all_objects();

	if (platedObject[RESOLVER_NED].name) {
		append_ips_row("NED",    platedObject[RESOLVER_NED].name);
		anyobject = TRUE;
	}
	if (platedObject[RESOLVER_SIMBAD].name) {
		append_ips_row("Simbad", platedObject[RESOLVER_SIMBAD].name);
		anyobject = TRUE;
	}
	if (platedObject[RESOLVER_VIZIER].name) {
		append_ips_row("VizieR", platedObject[RESOLVER_VIZIER].name);
		anyobject = TRUE;
	}
	if (platedObject[RESOLVER_LOCAL].name) {
		append_ips_row("Local",  platedObject[RESOLVER_LOCAL].name);
		anyobject = TRUE;
	}
	if (!anyobject) {
		append_ips_row("N/A", "No object found");
	}
}

static void unselect_all_items() {
	if (selection_model)
		gtk_single_selection_set_selected(selection_model, GTK_INVALID_LIST_POSITION);
}

static void add_object_in_tree_view(const gchar *object) {
	if (!object || object[0] == '\0')
		return;

	set_cursor_waiting(TRUE);
	gboolean found = FALSE;

	cat_item *local_obj = search_in_annotations_by_name(object, NULL);
	if (local_obj) {	// always search for local first
		add_plated_from_annotations(local_obj);
		g_signal_handlers_block_by_func(selection_model, on_IPS_selection_changed, NULL);
		add_object_to_list();
		g_signal_handlers_unblock_by_func(selection_model, on_IPS_selection_changed, NULL);
		found = TRUE;
	} else {
		sky_object_query_args *args= init_sky_object_query();
		args->name = g_strdup(object);
		args->server = get_server_from_combobox();
		gchar *result = search_in_online_catalogs(args);
		if (result) {
			free_Platedobject();
			struct sky_object obj;
			parse_resolver_buffer(result, &obj);
			if (!has_nonzero_coords()) {
				free(result);
				set_cursor_waiting(FALSE);
				// the list is empty, it will just write "No object found" as the first entry
				g_signal_handlers_block_by_func(selection_model, on_IPS_selection_changed, NULL);
				add_object_to_list();
				g_signal_handlers_unblock_by_func(selection_model, on_IPS_selection_changed, NULL);
				siril_log_error(_("No object found\n"));
				return;
			}
			g_signal_handlers_block_by_func(selection_model, on_IPS_selection_changed, NULL);
			add_object_to_list();
			g_signal_handlers_unblock_by_func(selection_model, on_IPS_selection_changed, NULL);
			free(result);
			found = TRUE;
		}
		free_sky_object_query(args);
	}

	if (found) {
		/* select first object found in the list */
		if (selection_model && g_list_model_get_n_items(G_LIST_MODEL(list_IPS)) > 0)
			gtk_single_selection_set_selected(selection_model, 0);
		siril_toggle_set_active(GTK_WIDGET(sequseheadercoords), FALSE);
	}
	if (local_obj) {
		siril_catalog_free_item(local_obj);
		free(local_obj);
	}
	else control_window_switch_to_tab(OUTPUT_LOGS);
	set_cursor_waiting(FALSE);
}

/*****
 * CALLBACKS FUNCTIONS
 */

void on_GtkEntry_IPS_focal_changed(GtkEditable *editable, gpointer user_data) {
	update_resolution_field();
	has_focal = FALSE;
	remove_style_class(GTK_WIDGET(focalentry), "header_unset");
	siril_toggle_set_active(GTK_WIDGET(sequseheaderfocal), FALSE);
}

void on_GtkEntry_IPS_pixels_changed(GtkEditable *editable, gpointer user_data) {
	update_resolution_field();
	has_pixel = FALSE;
	remove_style_class(GTK_WIDGET(pixelentry), "header_unset");
	siril_toggle_set_active(GTK_WIDGET(sequseheaderpixel), FALSE);
}

void on_GtkEntry_IPS_insert_text(GtkEntry *entry, const gchar *text, gint length,
		gint *position, gpointer data) {
	GtkEditable *editable = GTK_EDITABLE(entry);
	int i, count = 0;

	gchar *result = g_strndup(text, length);

	for (i = 0; i < length; i++) {
		if (!g_ascii_isdigit(text[i]) && text[i] != '.')
			continue;
		result[count++] = text[i];
	}

	if (count > 0) {
		g_signal_handlers_block_by_func(G_OBJECT (editable),
				G_CALLBACK (on_GtkEntry_IPS_insert_text), data);
		gtk_editable_insert_text(editable, result, count, position);
		g_signal_handlers_unblock_by_func(G_OBJECT (editable),
				G_CALLBACK (on_GtkEntry_IPS_insert_text), data);
	}
	g_signal_stop_emission_by_name(G_OBJECT(editable), "insert_text");
	remove_style_class(GTK_WIDGET(entry), "header_unset");
	g_free(result);
}

void on_buttonIPS_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("astrometry_dialog");
}

static void on_IPS_selection_changed(GObject *sm, GParamSpec *pspec, gpointer user_data) {
	(void)pspec; (void)user_data;
	GtkSingleSelection *sel = GTK_SINGLE_SELECTION(sm);
	guint pos = gtk_single_selection_get_selected(sel);
	if (pos == GTK_INVALID_LIST_POSITION)
		return;
	SirilIPSRow *row = SIRIL_IPS_ROW(g_list_model_get_item(G_LIST_MODEL(list_IPS), pos));
	if (!row) return;
	{
		const gchar *res = row->resolver;
		int selected_item;

		if (!g_strcmp0(res, "NED")) {
			selected_item = 0;
		} else if (!g_strcmp0(res, "Simbad")) {
			selected_item = 1;
		} else if (!g_strcmp0(res, "VizieR")) {
			selected_item = 2;
		} else if (!g_strcmp0(res, "Local")) {
			selected_item = 3;
		} else {
			selected_item = -1;
		}

		if (selected_item >= 0) {
			if (platedObject[selected_item].world_cs) {
				update_coordinates(platedObject[selected_item].world_cs);
				has_coords = FALSE; // the user changed the header coords, we assume it's because they were wrong
			} else {
				char *msg = siril_log_message(_("There are no available coordinates with this name, try with another name\n"));
				siril_message_dialog(GTK_MESSAGE_WARNING, _("No coordinates"), msg);
			}
		}
	}
	g_object_unref(row);
}

void on_GtkButton_IPS_metadata_clicked(GtkButton *button, gpointer user_data) {
	change_entry_colors_to_set();
	if (!has_any_keywords()) {
		char *msg = siril_log_message(_("There are no keywords stored in the FITS header.\n"));
		siril_message_dialog(GTK_MESSAGE_WARNING, _("No metadata"), msg);
		change_entry_colors_to_unset();
	} else {
		update_image_parameters_GUI();
	}
	siril_log_debug("metadata loaded\n");
}

void on_GtkButtonIPS_clicked(GtkButton *button, gpointer user_data) {
	add_object_in_tree_view(gtk_editable_get_text(searchentry));
}

void on_buttonIPS_ok_clicked(GtkButton *button, gpointer user_data) {
	start_image_plate_solve();
}

void on_GtkSearchIPS_activate(GtkEntry *entry, gpointer user_data) {
	add_object_in_tree_view(gtk_editable_get_text(GTK_EDITABLE(entry)));
}

void on_GtkCheckButton_Mag_Limit_toggled(GtkCheckButton *button, gpointer user) {
	gtk_widget_set_sensitive(GTK_WIDGET(magspin), !siril_toggle_get_active(GTK_WIDGET(button)));
}

void on_GtkCheckButton_OnlineCat_toggled(GtkCheckButton *button, gpointer user) {
	gtk_drop_down_set_selected(catalogbox, local_gaia_available() ? 2 : 1);
	gtk_widget_set_sensitive(GTK_WIDGET(catalogbox), !siril_toggle_get_active(GTK_WIDGET(button)));
	on_comboastro_catalog_changed(NULL, NULL, NULL);
}

void on_GtkCheckButton_solveseq_toggled(GtkCheckButton *button, gpointer user) {
	gboolean solveseq = gtk_widget_get_visible(GTK_WIDGET(seqsolvebutton)) && siril_toggle_get_active(GTK_WIDGET(seqsolvebutton));
	gboolean shownocache = FALSE;
	if (!gtk_drop_down_get_selected(solverbox)) { // SOLVER_SIRIL
		gboolean uselocal = use_local_catalogue();
		shownocache = (!uselocal) && (has_coords || has_pixel || has_focal);
	}
	gtk_widget_set_visible(GTK_WIDGET(IPSbox_seq_info), solveseq);
	gtk_widget_set_visible(GTK_WIDGET(seqnocache), solveseq && shownocache);
	gtk_widget_set_visible(GTK_WIDGET(sequseforreg), solveseq && com.seq.type == SEQ_REGULAR);
}

void on_GtkCheckButton_nonear_toggled(GtkCheckButton *button, gpointer user) {
	gtk_widget_set_sensitive(GTK_WIDGET(radiusspin), !siril_toggle_get_active(GTK_WIDGET(nonearbutton)));
}

void on_spinbuttoncoords_changed(GtkSpinButton* button, gpointer user) {
	has_coords = FALSE;
	remove_style_class(GTK_WIDGET(button), "header_unset");
}

void on_entrycoords_changed(GtkEditable *entry, gpointer user) {
	has_coords = FALSE;
	remove_style_class(GTK_WIDGET(entry), "header_unset");
}

void on_togglecoords_changed(GtkCheckButton *button, gpointer user) {
	has_coords = FALSE;
}

void on_GtkCheckButton_blindpos_toggled(GtkCheckButton *button, gpointer user) {
	gtk_widget_set_sensitive(GTK_WIDGET(radiusspin), !siril_toggle_get_active(GTK_WIDGET(blindposbutton)));
}

static void set_filter(SirilFileChooser *fc) {
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, _("wcs files (*.wcs)"));
	gtk_file_filter_add_pattern(f, "*.wcs");
	siril_fc_add_filter(fc, f, TRUE);
	g_object_unref(f);
}

void on_distomaster_save_button_clicked(GtkButton *button, gpointer user_data) {
	gchar *filename = NULL;

	if (sequence_is_loaded()) {
		filename = g_strdup_printf("%s.wcs", com.seq.seqname);
	} else {
		gchar *basename = g_path_get_basename(com.uniq->filename);
		char *root = remove_ext_from_filename(basename);
		filename = g_strdup_printf("%s.wcs", root);
		g_free(basename);
		free(root);
	}

	SirilFileChooser *fc = siril_fc_save(astrometry_dialog, GTK_FILE_CHOOSER_ACTION_SAVE);
	siril_fc_set_current_folder_path(fc, com.wd);
	siril_fc_set_select_multiple(fc, FALSE);
	siril_fc_set_current_name(fc, filename);
	set_filter(fc);

	if (siril_fc_run(fc) == GTK_RESPONSE_ACCEPT) {
		gchar *file = siril_fc_get_filename(fc);
		if (file) {
			gtk_editable_set_text(GTK_EDITABLE(distomaster_entry), file);
			gtk_editable_set_position(GTK_EDITABLE(distomaster_entry), -1);
			g_free(file);
			siril_toggle_set_active(GTK_WIDGET(masterbutton), TRUE);
		}
	}
	siril_fc_destroy(fc);
	g_free(filename);
}

void open_astrometry_dialog() {
	if (single_image_is_loaded() || sequence_is_loaded()) {
		initialize_ips_dialog();
		siril_open_dialog("astrometry_dialog");
	}
}

void close_astrometry_dialog() {
	siril_close_dialog("astrometry_dialog");
}

int fill_plate_solver_structure_from_GUI(struct astrometry_data *args) {
	args->solver = gtk_drop_down_get_selected(solverbox);
	gboolean is_siril = !args->solver;
	args->for_sequence = siril_toggle_get_active(GTK_WIDGET(seqsolvebutton)) && sequence_is_loaded();
	if (!args->for_sequence) {
		args->fit = gfit;
		args->manual = is_detection_manual();
		args->verbose = TRUE;
		args->flip_image = flip_image_after_ps() && (!sequence_is_loaded() || com.seq.current == RESULT_IMAGE);
		args->numthreads = com.max_thread;
	} else {
		args->force = !siril_toggle_get_active(GTK_WIDGET(seqskipsolved));
		args->update_reg = siril_toggle_get_active(GTK_WIDGET(sequseforreg)) && gtk_widget_get_visible(GTK_WIDGET(sequseforreg)); // not visible for FITSEQ and SER
		args->sfargs = calloc(1, sizeof(struct starfinder_data));
		args->sfargs->im.from_seq = &com.seq;
		args->sfargs->layer = (gfit->naxes[2] == 1) ? RLAYER : GLAYER;
		args->sfargs->keep_stars = TRUE;
		args->sfargs->save_to_file = com.selection.w == 0 || com.selection.h == 0; // TODO make this a pref
		args->sfargs->max_stars_fitted = BRIGHTEST_STARS;
	}
	args->downsample = is_downsample_activated();
	args->trans_order = get_order();
	args->pixel_size = get_pixel();
	args->focal_length = get_focal();
	if (is_save_disto_activated())
		args->distofilename = g_strdup(gtk_editable_get_text(GTK_EDITABLE(distomaster_entry)));
	SirilWorldCS *catalog_center = get_center_of_catalog();
	gboolean no_coords = siril_world_cs_get_alpha(catalog_center) == 0.0 &&
			siril_world_cs_get_delta(catalog_center) == 0.0;
	if (is_siril) { // SOLVER_SIRIL
		if (no_coords) {
			siril_message_dialog(GTK_MESSAGE_WARNING, _("No coordinates"),
					_("Please enter object coordinates."));
			siril_world_cs_unref(catalog_center);
			return 1;
		}
		int cat = gtk_drop_down_get_selected(catalogbox);
		gboolean autocat = siril_toggle_get_active(GTK_WIDGET(autocatbutton));
		gboolean uselocal = use_local_catalogue();
		args->autocrop = is_autocrop_activated();
		get_mag_settings_from_GUI(&args->mag_mode, &args->magnitude_arg);
		if (uselocal && !siril_toggle_get_active(GTK_WIDGET(nonearbutton))) {
			args->searchradius = gtk_spin_button_get_value(radiusspin);
		}
		args->cat_center = catalog_center;
		siril_cat_index cat_index =  autocat ? CAT_AUTO : cat;
		args->ref_stars = siril_catalog_new(cat_index);
		args->ref_stars->center_ra = siril_world_cs_get_alpha(catalog_center);
		args->ref_stars->center_dec = siril_world_cs_get_delta(catalog_center);
		if (com.selection.w != 0 && com.selection.h != 0 && (!args->for_sequence || !com.seq.is_variable)) {
			// we can't use selection for variable size sequences
			args->solvearea = com.selection;
			args->autocrop = FALSE; // we force the selection instead of autocrop
		}
		if (args->for_sequence) {
			// we solve each image individually if:
			// - we use local catalogues
			// - or user has selected nocache and there's at least one of the 3 metadata present in gfit header
			args->nocache = uselocal || (siril_toggle_get_active(GTK_WIDGET(seqnocache)) && (has_coords || has_focal || has_pixel));
			args->forced_metadata[FORCED_CENTER] = !siril_toggle_get_active(GTK_WIDGET(sequseheadercoords));
			args->forced_metadata[FORCED_PIXEL] = !siril_toggle_get_active(GTK_WIDGET(sequseheaderpixel));
			args->forced_metadata[FORCED_FOCAL] = !siril_toggle_get_active(GTK_WIDGET(sequseheaderfocal));
		}
	} else { //SOLVER_LOCALASNET
		args->asnet_checked = TRUE;
		args->autocrop = FALSE;
		if (com.selection.w != 0 && com.selection.h != 0)
			siril_log_message(_("Selection is not used with the astrometry.net solver\n"));
		args->asnet_blind_pos = siril_toggle_get_active(GTK_WIDGET(blindposbutton)) || no_coords;
		args->asnet_blind_res = siril_toggle_get_active(GTK_WIDGET(blindresbutton));
		if (!args->asnet_blind_pos) {
			args->searchradius = gtk_spin_button_get_value(radiusspin);
			args->cat_center = catalog_center;
		} else {
			args->cat_center = NULL;
			siril_world_cs_unref(catalog_center);
		}
		if (single_image_is_loaded() && com.uniq && com.uniq->filename) {
			args->filename = g_strdup(com.uniq->filename);
		} else if (sequence_is_loaded()) {
			args->filename = g_strdup_printf("%s%.5d", com.seq.seqname, com.seq.current + 1);
		} else {
			args->filename = g_strdup_printf("image");
		}
		if (args->for_sequence) {
			// we solve each image individually with asnet
			args->nocache = TRUE;
			if (!args->asnet_blind_pos)
				args->forced_metadata[FORCED_CENTER] = !has_coords;
			if (!args->asnet_blind_res) {
				args->forced_metadata[FORCED_PIXEL] = !has_pixel;
				args->forced_metadata[FORCED_FOCAL] = !has_focal;
			}
		}
	}
	if (!args->for_sequence) // this will be done in the prepare_hook
		process_plate_solver_input(args);
	return 0;
}

gboolean confirm_delete_wcs_keywords(fits *fit) {
	gboolean erase = TRUE;

	if (has_wcs(fit)) {
		erase = siril_confirm_dialog(_("Astrometric solution detected"),
				_("The astrometric solution contained in "
				"the image will be erased by the geometric transformation and no undo "
				"will be possible."), _("Process"));
	}
	return erase;
}

void init_astrometry() {
	load_all_ips_statics();
	reset_astrometry_checks();
	// Prefer Gaia to NOMAD and local to remote
	gtk_drop_down_set_selected(catalogbox, local_gaia_available() ? 2 : have_local_cat ? 1 : 2);
}

void on_comboastro_catalog_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	(void)obj;
	(void)pspec;
	if (!use_local_catalogue())
		gtk_label_set_text(cataloglabel, _("(online catalogue)"));
	else
		gtk_label_set_text(cataloglabel, _("(local catalogue)"));
	on_comboastro_solver_changed(NULL, NULL, NULL);
}

void on_comboastro_order_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	(void)obj;
	(void)pspec;
	gboolean enable = get_order() > 1;
	gtk_widget_set_sensitive(GTK_WIDGET(masterbutton), enable);
	gtk_widget_set_sensitive(GTK_WIDGET(distomaster_entry), enable);
	gtk_widget_set_sensitive(GTK_WIDGET(distomaster_save_button), enable);
	if (com.pref.prepro.use_disto_lib && com.pref.prepro.disto_lib && com.pref.prepro.disto_lib[0] != '\0') {
		gtk_editable_set_text(GTK_EDITABLE(distomaster_entry), com.pref.prepro.disto_lib);
	}
}

void on_comboastro_solver_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	(void)obj;
	(void)pspec;
	gboolean is_siril = !gtk_drop_down_get_selected(solverbox);
	gtk_widget_set_visible(GTK_WIDGET(autocropbutton), is_siril);
	gtk_widget_set_visible(GTK_WIDGET(blindposbutton), !is_siril);
	gtk_widget_set_visible(GTK_WIDGET(blindresbutton), !is_siril);
	gtk_widget_set_visible(GTK_WIDGET(cataloguesexp), is_siril);
	if (is_siril) {
		gboolean uselocal = use_local_catalogue();
		gtk_widget_set_visible(GTK_WIDGET(radiuslabel), uselocal);
		gtk_widget_set_visible(GTK_WIDGET(radiusspin), uselocal);
		gtk_widget_set_visible(GTK_WIDGET(nonearbutton), uselocal);
		on_GtkCheckButton_nonear_toggled(NULL, NULL);
	} else {
		gtk_widget_set_visible(GTK_WIDGET(nonearbutton), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(radiuslabel), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(radiusspin), TRUE);
		on_GtkCheckButton_blindpos_toggled(NULL, NULL);
	}
	on_GtkCheckButton_solveseq_toggled(NULL, NULL);
}

gboolean end_platesolve_sequence(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *) p;
	if (args->has_output && args->load_new_sequence &&
			args->new_seq_prefix && !args->retval) {
		gchar *basename = g_path_get_basename(args->seq->seqname);
		gchar *seqname = g_strdup_printf("%s%s.seq", args->new_seq_prefix, basename);
		check_seq();
		update_sequences_list(seqname);
		g_free(seqname);
		g_free(basename);
	}
	if (!check_seq_is_comseq(args->seq))
		free_sequence(args->seq, TRUE);
	else if (!args->retval && !args->has_output) {
		gchar *seqname = NULL;
		if (g_str_has_suffix(com.seq.seqname, ".seq"))
			seqname = g_strdup(com.seq.seqname);
		else
			seqname = g_strdup_printf("%s.seq", com.seq.seqname);
		set_seq(seqname);
		g_free(seqname);
	}
	free(p);	
	return end_generic(NULL);
}

gboolean astrometry_hide_on_delete(GtkWidget *widget) {
	mark_imgproc_dialog_closed();
	gtk_widget_set_visible(widget, FALSE);
	return TRUE;
}
