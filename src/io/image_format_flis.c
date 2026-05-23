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
 * FLIS (FITS Layered Image Specification) I/O
 *
 * Loads and saves FLIS files to/from com.uniq.  The FLIS format extends
 * standard FITS with a binary metadata table (FLIS_META) and per-layer
 * image HDUs, allowing multi-layer astrophotography editing state to be
 * round-tripped through a single FITS-compatible file.
 *
 * Functions from image_format_fits.c are reused wherever possible to
 * avoid duplication:
 *   read_fits_with_convert()   — layer pixel data reading
 *   save_opened_fits()         — layer pixel data writing
 *   write_icc_profile_to_fptr()— ICC profile HDU writing
 *   read_icc_profile_from_fptr()— ICC profile HDU reading
 *   siril_fits_compress()      — tile compression setup
 *   siril_fits_create_diskfile()— file creation
 *   manage_bitpix()            — BZERO/BSCALE handling
 *   report_fits_error()        — CFITSIO error reporting
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_log.h"

/* siril_debug_print no longer exists in mainline master; the FLIS
 * code from the flis branch relied on it for diagnostic output.
 * Provide a local no-op shim so the call sites remain syntactically
 * intact and can be promoted to a real logger later. */
#define siril_debug_print(...) ((void)0)

#include "core/icc_profile.h"
#include "core/masks.h"
#include "core/undo.h"
#include "algos/statistics.h"
#include "core/gui_iface.h"
#include "image_format_fits.h"
#include "image_format_flis.h"
#include "flis_compose.h"

/* -----------------------------------------------------------------------
 * FLIS_META binary table column definitions.
 *
 * The METADATA column uses the 1PA variable-length ASCII array format as
 * described in the FLIS spec.  Each row stores a semicolon-delimited
 * key=value string of arbitrary length in the FITS heap.  On read the
 * actual byte count is queried with fits_read_descript and a correctly
 * sized buffer is g_malloc'd per row; callers must free each
 * row->metadata before freeing the rows array.
 * ----------------------------------------------------------------------- */
#define FLIS_META_NCOLS      14

static const char *FLIS_COL_NAMES[FLIS_META_NCOLS] = {
    "ITEM_ID", "ITEM_TYPE", "HDU_INDEX", "PARENT_ID", "LAYER_ORDER",
    "LAYER_NAME", "COLOR_MDL", "BLEND_MODE", "OPACITY", "VISIBLE",
    "POSITION_X", "POSITION_Y", "METADATA", "GROUP_ID"
};
static const char *FLIS_COL_FMTS[FLIS_META_NCOLS] = {
    "1J", "8A", "1J", "1J", "1J",
    "32A", "4A", "16A", "1E", "1L",
    "1J", "1J", "1PA", "1J"
};
static const char *FLIS_COL_UNITS[FLIS_META_NCOLS] = {
    "", "", "", "", "",
    "", "", "", "", "",
    "pix", "pix", "", ""
};

/* Column index constants (1-based for CFITSIO) */
#define COL_ITEM_ID      1
#define COL_ITEM_TYPE    2
#define COL_HDU_INDEX    3
#define COL_PARENT_ID    4
#define COL_LAYER_ORDER  5
#define COL_LAYER_NAME   6
#define COL_COLOR_MDL    7
#define COL_BLEND_MODE   8
#define COL_OPACITY      9
#define COL_VISIBLE     10
#define COL_POSITION_X  11
#define COL_POSITION_Y  12
#define COL_METADATA    13
#define COL_GROUP_ID    14

/* ITEM_TYPE string constants */
#define FLIS_TYPE_LAYER  "LAYER"
#define FLIS_TYPE_MASK   "MASK"
#define FLIS_TYPE_LMASK  "LMASK"
#define FLIS_TYPE_GROUP  "GROUP"

/* -----------------------------------------------------------------------
 * Scratch structure for one FLIS_META table row during load.
 * ----------------------------------------------------------------------- */
typedef struct {
    gint    item_id;
    char    item_type[9];    /* 8A + NUL */
    gint    hdu_index;       /* 1-based CFITSIO HDU number */
    gint    parent_id;
    gint    layer_order;
    char    layer_name[33];  /* 32A + NUL */
    char    color_mdl[5];    /* 4A + NUL  */
    char    blend_mode[17];  /* 16A + NUL */
    gfloat  opacity;
    gboolean visible;
    gint    position_x;
    gint    position_y;
    char   *metadata;        /* 1PA variable-length string, g_malloc'd */
    gint    group_id;
} flis_meta_row_t;

/* ===================================================================== */
/* Blend mode string conversion                                          */
/* ===================================================================== */

static const struct { flis_blend_mode_t mode; const char *name; }
    BLEND_TABLE[] = {
        { FLIS_BLEND_NORMAL,      "NORMAL"      },
        { FLIS_BLEND_MULTIPLY,    "MULTIPLY"    },
        { FLIS_BLEND_SCREEN,      "SCREEN"      },
        { FLIS_BLEND_OVERLAY,     "OVERLAY"     },
        { FLIS_BLEND_SOFT_LIGHT,  "SOFT_LIGHT"  },
        { FLIS_BLEND_HARD_LIGHT,  "HARD_LIGHT"  },
        { FLIS_BLEND_COLOR_DODGE, "COLOR_DODGE" },
        { FLIS_BLEND_COLOR_BURN,  "COLOR_BURN"  },
        { FLIS_BLEND_DARKEN,      "DARKEN"      },
        { FLIS_BLEND_LIGHTEN,     "LIGHTEN"     },
        { FLIS_BLEND_DIFFERENCE,  "DIFFERENCE"  },
        { FLIS_BLEND_EXCLUSION,   "EXCLUSION"   },
        { FLIS_BLEND_HUE,         "HUE"         },
        { FLIS_BLEND_SATURATION,  "SATURATION"  },
        { FLIS_BLEND_COLOR,       "COLOR"       },
        { FLIS_BLEND_LUMINOSITY,  "LUMINOSITY"  },
        { FLIS_BLEND_CHROMA,        "CHROMA"    },
        { FLIS_BLEND_PASS_THROUGH,  "PASSTHRU"  },
    };

static const char *blend_mode_to_str(flis_blend_mode_t mode) {
    for (int i = 0; i < G_N_ELEMENTS(BLEND_TABLE); i++) {
        if (BLEND_TABLE[i].mode == mode)
            return BLEND_TABLE[i].name;
    }
    return "NORMAL";
}

static flis_blend_mode_t str_to_blend_mode(const char *s) {
    for (int i = 0; i < G_N_ELEMENTS(BLEND_TABLE); i++) {
        if (!g_ascii_strcasecmp(s, BLEND_TABLE[i].name))
            return BLEND_TABLE[i].mode;
    }
    siril_log_message(_("FLIS: unknown blend mode '%s', defaulting to NORMAL\n"), s);
    return FLIS_BLEND_NORMAL;
}

/* ===================================================================== */
/* Metadata string encoding / decoding                                   */
/* ===================================================================== */

/*
 * Build the semicolon-delimited METADATA string for one layer.
 * Format:  KEY=VALUE;KEY=VALUE;...
 * Returns a g_malloc'd string that the caller must g_free().
 */
static gchar *build_metadata_string(const flis_layer_t *layer) {
    GString *s = g_string_new(NULL);

    if (layer->created)
        g_string_append_printf(s, "CREATED=%s;", layer->created);
    if (layer->modified)
        g_string_append_printf(s, "MODIFIED=%s;", layer->modified);
    if (layer->locked)
        g_string_append(s, "LOCKED=T;");
    if (layer->has_tint)
        g_string_append_printf(s, "LAYER_COLOR=%.6f,%.6f,%.6f;",
                               layer->layer_tint.r,
                               layer->layer_tint.g,
                               layer->layer_tint.b);

    /* Strip trailing semicolon */
    if (s->len > 0 && s->str[s->len - 1] == ';')
        g_string_truncate(s, s->len - 1);

    return g_string_free(s, FALSE);
}

/*
 * Parse a METADATA string into the relevant fields of a flis_layer_t.
 * Unknown keys are silently ignored (per spec Section 8).
 */
static void parse_metadata_string(const char *meta, flis_layer_t *layer) {
    if (!meta || meta[0] == '\0')
        return;

    gchar **pairs = g_strsplit(meta, ";", -1);
    for (int i = 0; pairs[i]; i++) {
        gchar **kv = g_strsplit(pairs[i], "=", 2);
        if (!kv[0] || !kv[1]) { g_strfreev(kv); continue; }

        const char *k = g_strstrip(kv[0]);
        const char *v = g_strstrip(kv[1]);

        if (!g_ascii_strcasecmp(k, "CREATED")) {
            g_free(layer->created);
            layer->created = g_strdup(v);
        } else if (!g_ascii_strcasecmp(k, "MODIFIED")) {
            g_free(layer->modified);
            layer->modified = g_strdup(v);
        } else if (!g_ascii_strcasecmp(k, "LOCKED")) {
            layer->locked = (v[0] == 'T' || v[0] == 't');
        } else if (!g_ascii_strcasecmp(k, "LAYER_COLOR")) {
            double r, g, b;
            if (sscanf(v, "%lf,%lf,%lf", &r, &g, &b) == 3) {
                layer->has_tint = TRUE;
                layer->layer_tint.r = r;
                layer->layer_tint.g = g;
                layer->layer_tint.b = b;
            }
        }
        /* MASKGEN and other keys are preserved on save but not acted upon here */
        g_strfreev(kv);
    }
    g_strfreev(pairs);
}

/* ===================================================================== */
/* Thumbnail generation                                                  */
/* ===================================================================== */

/*
 * Generate an 8-bit RGB or greyscale thumbnail from the base layer.
 *
 * A proper composite thumbnail (all visible layers, layer masks applied)
 * requires a compositor that does not yet exist.  For now we use the base
 * layer as the thumbnail source.  The thumbnail is written bottom-to-top
 * (FITS convention) so that FITS viewers display it correctly.
 *
 * @fit:       source layer (base layer, lowest layer_order)
 * @out_w:     populated with thumbnail width
 * @out_h:     populated with thumbnail height
 * Returns:    heap-allocated uint8_t array (caller must free()), or NULL.
 */
static uint8_t *generate_thumbnail(const fits *fit,
                                   long *out_w, long *out_h) {
    if (!fit) return NULL;

    /* Compute scale factor to fit within FLIS_THUMBNAIL_MAX x FLIS_THUMBNAIL_MAX */
    double scale = 1.0;
    if (fit->rx > FLIS_THUMBNAIL_MAX || fit->ry > FLIS_THUMBNAIL_MAX) {
        double sx = (double)FLIS_THUMBNAIL_MAX / fit->rx;
        double sy = (double)FLIS_THUMBNAIL_MAX / fit->ry;
        scale = MIN(sx, sy);
    }

    long tw = (long)round(fit->rx * scale);
    long th = (long)round(fit->ry * scale);
    if (tw < 1) tw = 1;
    if (th < 1) th = 1;

    long nchans = fit->naxes[2] > 0 ? fit->naxes[2] : 1;
    /* FLIS thumbnail is always RGB (3 channels) in FITS convention:
     * planes stored as separate planes in a 3D array. */
    long out_chans = (nchans >= 3) ? 3 : 1;
    size_t npix_thumb = (size_t)(tw * th);

    uint8_t *thumb = malloc(npix_thumb * out_chans);
    if (!thumb) { PRINT_ALLOC_ERR; return NULL; }

    /* Compute min/max for auto-stretch — sample every pixel */
    double fmin = DBL_MAX, fmax = -DBL_MAX;
    size_t npix_src = (size_t)(fit->rx * fit->ry);

    if (fit->type == DATA_FLOAT && fit->fdata) {
        for (size_t i = 0; i < npix_src; i++) {
            if (fit->fdata[i] < fmin) fmin = fit->fdata[i];
            if (fit->fdata[i] > fmax) fmax = fit->fdata[i];
        }
    } else if (fit->type == DATA_USHORT && fit->data) {
        for (size_t i = 0; i < npix_src; i++) {
            if (fit->data[i] < fmin) fmin = fit->data[i];
            if (fit->data[i] > fmax) fmax = fit->data[i];
        }
    } else {
        free(thumb);
        return NULL;
    }
    if (fmax <= fmin) fmax = fmin + 1.0;
    double range = fmax - fmin;

    /* Nearest-neighbour downscale into output buffer */
    for (long c = 0; c < out_chans; c++) {
        /* Choose source channel: if mono source and RGB output, use ch 0 */
        long src_chan = (nchans >= 3) ? c : 0;
        size_t chan_off_src = (size_t)(src_chan * npix_src);
        size_t chan_off_dst = (size_t)(c * npix_thumb);

        for (long y = 0; y < th; y++) {
            long sy = (long)((y / scale) + 0.5);
            if (sy >= fit->ry) sy = fit->ry - 1;
            for (long x = 0; x < tw; x++) {
                long sx = (long)((x / scale) + 0.5);
                if (sx >= fit->rx) sx = fit->rx - 1;

                double v;
                size_t src_idx = chan_off_src + (size_t)(sy * fit->rx + sx);
                if (fit->type == DATA_FLOAT)
                    v = fit->fdata[src_idx];
                else
                    v = (double)fit->data[src_idx];

                /* Stretch to [0,255] */
                double stretched = (v - fmin) / range * 255.0;
                if (stretched < 0.0)   stretched = 0.0;
                if (stretched > 255.0) stretched = 255.0;
                thumb[chan_off_dst + (size_t)(y * tw + x)] = (uint8_t)stretched;
            }
        }
    }

    *out_w = tw;
    *out_h = th;
    return thumb;
}

/* ===================================================================== */
/* HDU writing helpers                                                   */
/* ===================================================================== */

/*
 * Write HDU 0: the 8-bit composite thumbnail.
 * The fptr must be positioned at HDU 1 (the primary HDU) on entry;
 * it remains there on successful return.
 */
static int write_thumbnail_hdu(fitsfile *fptr, const fits *base_fit,
                               long canvas_w, long canvas_h,
                               gboolean icc_present) {
    int status = 0;
    long tw = 0, th = 0;

    uint8_t *thumb = generate_thumbnail(base_fit, &tw, &th);
    if (!thumb) {
        siril_log_warning(_("FLIS: thumbnail generation failed, writing empty HDU\n"));
        /* Fall through to write a degenerate 1x1 placeholder */
        tw = 1; th = 1;
        thumb = calloc(1, 1);
        if (!thumb) return 1;
    }

    long nchans = (base_fit && base_fit->naxes[2] >= 3) ? 3 : 1;
    int naxis = (nchans == 3) ? 3 : 2;
    long naxes[3] = { tw, th, nchans };

    /* Create the primary HDU.  fits_resize_img() would be wrong here because
     * the file was just created by fits_create_diskfile() and has no primary
     * HDU yet — calling resize on an empty file produces the "Extension
     * doesn't start with SIMPLE or XTENSION" error.  fits_create_img() is
     * the correct call for writing the very first HDU into a new file. */
    if (fits_create_img(fptr, BYTE_IMG, naxis, naxes, &status)) {
        report_fits_error(status);
        free(thumb);
        return 1;
    }

    /* Write FLIS identification keywords */
    int ival;
    fits_write_key(fptr, TLOGICAL, "FLIS",    &(int){1},       "FITS Layered Image Specification", &status);
    fits_write_key(fptr, TSTRING,  "FLISVER", FLIS_VERSION_STRING, "FLIS version",                 &status);
    fits_write_key(fptr, TLOGICAL, "FLISICC", &(int){icc_present}, "ICC profile present",          &status);
    fits_write_key(fptr, TLOGICAL, "FLISTHMB",&(int){1},       "Primary HDU contains thumbnail",   &status);
    ival = (int)canvas_w;
    fits_write_key(fptr, TINT,     "FLISSIZX",&ival,           "Original canvas width in pixels",  &status);
    ival = (int)canvas_h;
    fits_write_key(fptr, TINT,     "FLISSIZY",&ival,           "Original canvas height in pixels", &status);

    /* Capability flags */
    fits_write_key(fptr, TSTRING,  "FLISIMPL","Siril",         "FLIS implementation",              &status);
    fits_write_key(fptr, TLOGICAL, "FLISCORE",&(int){1},       "Core FLIS features supported",     &status);
    ival = FLIS_BLND_ALL;
    fits_write_key(fptr, TINT,     "FLISBLND",&ival,           "Supported blend modes bitmask",    &status);
    fits_write_key(fptr, TLOGICAL, "FLISLMSK",&(int){1},       "Layer mask (LMASK) supported",     &status);
    fits_write_key(fptr, TLOGICAL, "FLISEFF", &(int){0},       "Effects metadata supported",       &status);
    fits_write_key(fptr, TLOGICAL, "FLISICC", &(int){1},       "ICC profile supported",            &status);

    if (status) {
        report_fits_error(status);
        free(thumb);
        return 1;
    }

    long fpixel = 1;
    size_t npix = (size_t)(tw * th * nchans);
    if (fits_write_img(fptr, TBYTE, fpixel, npix, thumb, &status)) {
        report_fits_error(status);
        free(thumb);
        return 1;
    }
    free(thumb);
    return 0;
}

/*
 * Write the FLIS_META binary table HDU.
 *
 * @fptr:         open FITS file (positioned at last HDU on entry)
 * @layers:       GSList of flis_layer_t* sorted by layer_order
 * @first_data_hdu: 1-based HDU number of the first data (layer/mask) HDU
 *
 * The FLIS_META HDU is appended to the file.  HDU_INDEX values in the
 * table are pre-computed from first_data_hdu and the layer/mask structure.
 */
static int write_flis_meta_hdu(fitsfile *fptr, GSList *layers,
                               int first_data_hdu) {
    int status = 0;
    int nrows = 0;

    /* Count total rows: each layer + its optional lmask + its optional pmask */
    for (GSList *l = layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        nrows++;
        if (lay->lmask) nrows++;
        if (lay->fit && lay->fit->mask) nrows++;
    }
    /* Also count group rows */
    if (com.uniq && com.uniq->groups)
        nrows += (int)g_slist_length(com.uniq->groups);

    /* Create the binary table — cast away const for cfitsio API */
    if (fits_create_tbl(fptr, BINARY_TBL, nrows, FLIS_META_NCOLS,
                        (char **)FLIS_COL_NAMES, (char **)FLIS_COL_FMTS,
                        (char **)FLIS_COL_UNITS, "FLIS_META", &status)) {
        report_fits_error(status);
        return 1;
    }

    /* Walk layers in order, assigning HDU indices as we go */
    int row = 1;
    int next_hdu = first_data_hdu;

    for (GSList *l = layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        const char *color_mdl = (lay->fit && lay->fit->naxes[2] >= 3) ? "RGB" : "MONO";
        int layer_hdu = next_hdu++;
        int lmask_hdu = lay->lmask ? next_hdu++ : 0;
        int pmask_hdu = (lay->fit && lay->fit->mask) ? next_hdu++ : 0;

        /* --- Layer row --- */
        {
            gint id  = lay->item_id;
            gint idx = layer_hdu;
            gint pid = 0;
            gint ord = lay->layer_order;
            gint px  = lay->position_x;
            gint py  = lay->position_y;
            gfloat op = lay->opacity;
            int vis  = lay->visible ? 1 : 0;
            const char *blend = blend_mode_to_str(lay->blend_mode);
            gchar *meta = build_metadata_string(lay);
            const char *name = lay->layer_name ? lay->layer_name : "Layer";

            fits_write_col(fptr, TINT,     COL_ITEM_ID,     row, 1, 1, &id,              &status);
            fits_write_col(fptr, TSTRING,  COL_ITEM_TYPE,   row, 1, 1, &(char*){FLIS_TYPE_LAYER}, &status);
            fits_write_col(fptr, TINT,     COL_HDU_INDEX,   row, 1, 1, &idx,             &status);
            fits_write_col(fptr, TINT,     COL_PARENT_ID,   row, 1, 1, &pid,             &status);
            fits_write_col(fptr, TINT,     COL_LAYER_ORDER, row, 1, 1, &ord,             &status);
            fits_write_col(fptr, TSTRING,  COL_LAYER_NAME,  row, 1, 1, &name,            &status);
            fits_write_col(fptr, TSTRING,  COL_COLOR_MDL,   row, 1, 1, &color_mdl,       &status);
            fits_write_col(fptr, TSTRING,  COL_BLEND_MODE,  row, 1, 1, &blend,           &status);
            fits_write_col(fptr, TFLOAT,   COL_OPACITY,     row, 1, 1, &op,              &status);
            fits_write_col(fptr, TLOGICAL, COL_VISIBLE,     row, 1, 1, &vis,             &status);
            fits_write_col(fptr, TINT,     COL_POSITION_X,  row, 1, 1, &px,              &status);
            fits_write_col(fptr, TINT,     COL_POSITION_Y,  row, 1, 1, &py,              &status);
            if (meta && meta[0]) {
                fits_write_col(fptr, TSTRING, COL_METADATA, row, 1, 1, &meta,            &status);
            } else {
                const char *empty = "";
                fits_write_col(fptr, TSTRING, COL_METADATA, row, 1, 1, &empty,           &status);
            }
            g_free(meta);
            {
                gint gid = lay->group_id;
                fits_write_col(fptr, TINT, COL_GROUP_ID, row, 1, 1, &gid, &status);
            }
            row++;
        }

        /* --- Layer mask row (LMASK) --- */
        if (lay->lmask) {
            gint id  = lay->item_id + 10000; /* convention: lmask ID = layer ID + 10000 */
            gint idx = lmask_hdu;
            gint pid = lay->item_id;
            gint ord = 0;
            gint px  = 0, py = 0;
            gfloat op = 1.0f;
            int vis = 1;
            gchar *lname = g_strdup_printf("%s Layer Mask",
                           lay->layer_name ? lay->layer_name : "Layer");
            const char *empty = "";

            fits_write_col(fptr, TINT,     COL_ITEM_ID,     row, 1, 1, &id,              &status);
            fits_write_col(fptr, TSTRING,  COL_ITEM_TYPE,   row, 1, 1, &(char*){FLIS_TYPE_LMASK}, &status);
            fits_write_col(fptr, TINT,     COL_HDU_INDEX,   row, 1, 1, &idx,             &status);
            fits_write_col(fptr, TINT,     COL_PARENT_ID,   row, 1, 1, &pid,             &status);
            fits_write_col(fptr, TINT,     COL_LAYER_ORDER, row, 1, 1, &ord,             &status);
            fits_write_col(fptr, TSTRING,  COL_LAYER_NAME,  row, 1, 1, &lname,           &status);
            fits_write_col(fptr, TSTRING,  COL_COLOR_MDL,   row, 1, 1, &(char*){"MONO"}, &status);
            fits_write_col(fptr, TSTRING,  COL_BLEND_MODE,  row, 1, 1, &(char*){"NORMAL"}, &status);
            fits_write_col(fptr, TFLOAT,   COL_OPACITY,     row, 1, 1, &op,              &status);
            fits_write_col(fptr, TLOGICAL, COL_VISIBLE,     row, 1, 1, &vis,             &status);
            fits_write_col(fptr, TINT,     COL_POSITION_X,  row, 1, 1, &px,              &status);
            fits_write_col(fptr, TINT,     COL_POSITION_Y,  row, 1, 1, &py,              &status);
            fits_write_col(fptr, TSTRING,  COL_METADATA,    row, 1, 1, &empty,           &status);
            {
                gint gid = 0;
                fits_write_col(fptr, TINT, COL_GROUP_ID, row, 1, 1, &gid, &status);
            }
            g_free(lname);
            row++;
        }

        /* --- Processing mask row (MASK) --- */
        if (lay->fit && lay->fit->mask) {
            gint id  = lay->item_id + 20000; /* convention: pmask ID = layer ID + 20000 */
            gint idx = pmask_hdu;
            gint pid = lay->item_id;
            gint ord = 0;
            gint px  = 0, py = 0;
            gfloat op = 1.0f;
            int vis = 1;
            gchar *mname = g_strdup_printf("%s Processing Mask",
                           lay->layer_name ? lay->layer_name : "Layer");
            const char *empty = "";

            fits_write_col(fptr, TINT,     COL_ITEM_ID,     row, 1, 1, &id,              &status);
            fits_write_col(fptr, TSTRING,  COL_ITEM_TYPE,   row, 1, 1, &(char*){FLIS_TYPE_MASK}, &status);
            fits_write_col(fptr, TINT,     COL_HDU_INDEX,   row, 1, 1, &idx,             &status);
            fits_write_col(fptr, TINT,     COL_PARENT_ID,   row, 1, 1, &pid,             &status);
            fits_write_col(fptr, TINT,     COL_LAYER_ORDER, row, 1, 1, &ord,             &status);
            fits_write_col(fptr, TSTRING,  COL_LAYER_NAME,  row, 1, 1, &mname,           &status);
            fits_write_col(fptr, TSTRING,  COL_COLOR_MDL,   row, 1, 1, &(char*){"MONO"}, &status);
            fits_write_col(fptr, TSTRING,  COL_BLEND_MODE,  row, 1, 1, &(char*){"NORMAL"}, &status);
            fits_write_col(fptr, TFLOAT,   COL_OPACITY,     row, 1, 1, &op,              &status);
            fits_write_col(fptr, TLOGICAL, COL_VISIBLE,     row, 1, 1, &vis,             &status);
            fits_write_col(fptr, TINT,     COL_POSITION_X,  row, 1, 1, &px,              &status);
            fits_write_col(fptr, TINT,     COL_POSITION_Y,  row, 1, 1, &py,              &status);
            fits_write_col(fptr, TSTRING,  COL_METADATA,    row, 1, 1, &empty,           &status);
            {
                gint gid = 0;
                fits_write_col(fptr, TINT, COL_GROUP_ID, row, 1, 1, &gid, &status);
            }
            g_free(mname);
            row++;
        }

        if (status) {
            report_fits_error(status);
            return 1;
        }
    }

    /* Write GROUP rows */
    if (com.uniq && com.uniq->groups) {
        for (GSList *g = com.uniq->groups; g; g = g->next) {
            flis_group_t *grp = (flis_group_t *)g->data;
            if (!grp) continue;
            gint id   = grp->item_id;
            gint zero = 0;
            gfloat op = grp->opacity;
            int vis   = grp->visible ? 1 : 0;
            const char *name = grp->name ? grp->name : "Group";
            /* Build metadata string for group: COLLAPSED and timestamps */
            gchar *meta = g_strdup_printf("COLLAPSED=%d;CREATED=%s;MODIFIED=%s",
                grp->collapsed ? 1 : 0,
                grp->created  ? grp->created  : "",
                grp->modified ? grp->modified : "");
            fits_write_col(fptr, TINT,     COL_ITEM_ID,     row, 1, 1, &id,                        &status);
            fits_write_col(fptr, TSTRING,  COL_ITEM_TYPE,   row, 1, 1, &(char*){FLIS_TYPE_GROUP},  &status);
            fits_write_col(fptr, TINT,     COL_HDU_INDEX,   row, 1, 1, &zero,                      &status);
            fits_write_col(fptr, TINT,     COL_PARENT_ID,   row, 1, 1, &zero,                      &status);
            fits_write_col(fptr, TINT,     COL_LAYER_ORDER, row, 1, 1, &zero,                      &status);
            fits_write_col(fptr, TSTRING,  COL_LAYER_NAME,  row, 1, 1, &name,                      &status);
            fits_write_col(fptr, TSTRING,  COL_COLOR_MDL,   row, 1, 1, &(char*){"N/A"},            &status);
            const char *bm_str = blend_mode_to_str(grp->blend_mode);
            fits_write_col(fptr, TSTRING,  COL_BLEND_MODE,  row, 1, 1, &bm_str,                    &status);
            fits_write_col(fptr, TFLOAT,   COL_OPACITY,     row, 1, 1, &op,                        &status);
            fits_write_col(fptr, TLOGICAL, COL_VISIBLE,     row, 1, 1, &vis,                       &status);
            fits_write_col(fptr, TINT,     COL_POSITION_X,  row, 1, 1, &zero,                      &status);
            fits_write_col(fptr, TINT,     COL_POSITION_Y,  row, 1, 1, &zero,                      &status);
            if (meta && meta[0]) {
                fits_write_col(fptr, TSTRING, COL_METADATA, row, 1, 1, &meta,                      &status);
            } else {
                const char *empty = "";
                fits_write_col(fptr, TSTRING, COL_METADATA, row, 1, 1, &empty,                     &status);
            }
            fits_write_col(fptr, TINT,     COL_GROUP_ID,    row, 1, 1, &zero,                      &status);
            g_free(meta);

            if (status) { report_fits_error(status); return 1; }
            row++;
        }
    }

    return 0;
}

/*
 * Write one layer's pixel data as a new image HDU.
 * Temporarily borrows fptr into layer->fit to reuse save_opened_fits().
 * The flis_id and flis_type keywords are added after the pixel write.
 */
static int write_layer_hdu(fitsfile *fptr, flis_layer_t *layer) {
    int status = 0;
    fits *f = layer->fit;

    f->naxes[0] = f->rx;
    f->naxes[1] = f->ry;

    if (fits_create_img(fptr, f->bitpix, f->naxis, f->naxes, &status)) {
        report_fits_error(status);
        return 1;
    }

    /* BZERO/BSCALE for unsigned 16-bit: CFITSIO sets these automatically
     * when USHORT_IMG is used, but be explicit for FLIS conformance. */
    if (f->bitpix == USHORT_IMG) {
        double bzero = 32768.0, bscale = 1.0;
        fits_write_key(fptr, TDOUBLE, "BZERO",  &bzero,  "Unsigned 16-bit offset", &status);
        fits_write_key(fptr, TDOUBLE, "BSCALE", &bscale, "Scale factor",           &status);
        status = 0;
    }

    /* Temporarily borrow fptr so save_opened_fits() writes to this HDU */
    fitsfile *saved_fptr = f->fptr;
    f->fptr = fptr;
    int ret = save_opened_fits(f);
    f->fptr = saved_fptr;

    if (ret) return 1;

    /* Add FLIS-specific keywords */
    gint flis_id = layer->item_id;
    const char *name = layer->layer_name ? layer->layer_name : "Layer";
    fits_write_key(fptr, TINT,    "FLIS_ID",   &flis_id,        "FLIS item ID",          &status);
    fits_write_key(fptr, TSTRING, "FLIS_TYPE", FLIS_TYPE_LAYER, "FLIS HDU type",         &status);
    fits_write_key(fptr, TSTRING, "EXTNAME",   (void *)name,    "Extension name",        &status);

    if (status) { report_fits_error(status); return 1; }
    return 0;
}

/*
 * Write a layermask_t as a new greyscale image HDU.
 * @bitpix_fits: BYTE_IMG (8) for LMASK, FLOAT_IMG (-32) for MASK.
 * @cfitsio_type: TBYTE or TFLOAT, matching bitpix_fits.
 */
static int write_mask_hdu(fitsfile *fptr, const layermask_t *mask,
                          gint flis_id, const char *flis_type,
                          const char *extname,
                          int bitpix_fits, int cfitsio_type,
                          gboolean mask_active) {
    int status = 0;
    long naxes[2] = { (long)mask->w, (long)mask->h };

    if (fits_create_img(fptr, bitpix_fits, 2, naxes, &status)) {
        report_fits_error(status);
        return 1;
    }

    gint id = flis_id;
    fits_write_key(fptr, TINT,    "FLIS_ID",   &id,             "FLIS item ID",  &status);
    fits_write_key(fptr, TSTRING, "FLIS_TYPE", (void *)flis_type,"FLIS HDU type",&status);
    fits_write_key(fptr, TSTRING, "EXTNAME",   (void *)extname,  "Extension name",&status);

    /* MASK_ACT: mask active flag — present in both LMASK and MASK HDUs */
    {
        int act = mask_active ? 1 : 0;
        fits_write_key(fptr, TLOGICAL, "MASK_ACT", &act,
                       "Mask active in compositing", &status);
    }

    if (status) { report_fits_error(status); return 1; }

    long npix = (long)(mask->w * mask->h);
    if (fits_write_img(fptr, cfitsio_type, 1, npix, mask->data, &status)) {
        report_fits_error(status);
        return 1;
    }
    return 0;
}

/* ===================================================================== */
/* HDU reading helpers                                                   */
/* ===================================================================== */

/*
 * Read all rows from the FLIS_META binary table into a heap-allocated
 * array.  The file's active HDU is expected to be the FLIS_META table
 * on entry; it is left unchanged on return.
 *
 * @nrows_out: populated with the number of rows read.
 * Returns: heap-allocated flis_meta_row_t array, or NULL on error.
 *          Caller must free() the array.
 */
static flis_meta_row_t *read_flis_meta(fitsfile *fptr, long *nrows_out) {
    int status = 0;
    long nrows = 0;
    fits_get_num_rows(fptr, &nrows, &status);
    if (status || nrows <= 0) {
        report_fits_error(status);
        return NULL;
    }

    flis_meta_row_t *rows = calloc(nrows, sizeof(flis_meta_row_t));
    if (!rows) { PRINT_ALLOC_ERR; return NULL; }

    for (long r = 1; r <= nrows; r++) {
        flis_meta_row_t *row = &rows[r - 1];
        char *strptr;
        char anull = '\0';

        fits_read_col(fptr, TINT,    COL_ITEM_ID,     r, 1, 1, &(int){0},  &row->item_id,     NULL, &status);
        strptr = row->item_type;
        fits_read_col(fptr, TSTRING, COL_ITEM_TYPE,   r, 1, 1, &anull,     &strptr,            NULL, &status);
        fits_read_col(fptr, TINT,    COL_HDU_INDEX,   r, 1, 1, &(int){0},  &row->hdu_index,    NULL, &status);
        fits_read_col(fptr, TINT,    COL_PARENT_ID,   r, 1, 1, &(int){0},  &row->parent_id,    NULL, &status);
        fits_read_col(fptr, TINT,    COL_LAYER_ORDER, r, 1, 1, &(int){0},  &row->layer_order,  NULL, &status);
        strptr = row->layer_name;
        fits_read_col(fptr, TSTRING, COL_LAYER_NAME,  r, 1, 1, &anull,     &strptr,            NULL, &status);
        strptr = row->color_mdl;
        fits_read_col(fptr, TSTRING, COL_COLOR_MDL,   r, 1, 1, &anull,     &strptr,            NULL, &status);
        strptr = row->blend_mode;
        fits_read_col(fptr, TSTRING, COL_BLEND_MODE,  r, 1, 1, &anull,     &strptr,            NULL, &status);
        fits_read_col(fptr, TFLOAT,  COL_OPACITY,     r, 1, 1, &(float){1.0f}, &row->opacity,  NULL, &status);

        char vis_char = 'T';
        fits_read_col(fptr, TLOGICAL,COL_VISIBLE,     r, 1, 1, &anull,     &vis_char,          NULL, &status);
        row->visible = (vis_char != 'F' && vis_char != 0);

        fits_read_col(fptr, TINT,    COL_POSITION_X,  r, 1, 1, &(int){0},  &row->position_x,   NULL, &status);
        fits_read_col(fptr, TINT,    COL_POSITION_Y,  r, 1, 1, &(int){0},  &row->position_y,   NULL, &status);
        {
            /* Variable-length 1PA column: query actual byte count first,
             * then allocate an exact-fit buffer and read into it. */
            long vl_repeat = 0, vl_offset = 0;
            int vl_status = 0;
            fits_read_descript(fptr, COL_METADATA, r, &vl_repeat, &vl_offset, &vl_status);
            if (vl_status || vl_repeat < 0) vl_repeat = 0;
            row->metadata = g_malloc(vl_repeat + 1);
            row->metadata[0] = '\0';
            if (vl_repeat > 0) {
                strptr = row->metadata;
                fits_read_col(fptr, TSTRING, COL_METADATA, r, 1, 1, &anull, &strptr, NULL, &status);
                row->metadata[vl_repeat] = '\0'; /* guarantee NUL termination */
            }
        }
        fits_read_col(fptr, TINT, COL_GROUP_ID, r, 1, 1, &(int){0}, &row->group_id, NULL, &status);
        if (status == COL_NOT_FOUND || status == 219) {
            row->group_id = 0;
            status = 0;
        }

        if (status) {
            siril_log_warning(_("FLIS: error reading metadata row %ld: %d\n"), r, status);
            status = 0; /* continue with next row */
        }
    }

    *nrows_out = nrows;
    return rows;
}

/*
 * Load pixel data from the HDU at hdu_index (1-based) into a newly
 * allocated fits*.  The file's active HDU is restored on return.
 */
static fits *load_layer_from_hdu(fitsfile *fptr, int hdu_index) {
    int status = 0;
    int orig_hdu;
    fits_get_hdu_num(fptr, &orig_hdu);

    if (fits_movabs_hdu(fptr, hdu_index, NULL, &status)) {
        report_fits_error(status);
        return NULL;
    }

    fits *f = calloc(1, sizeof(fits));
    if (!f) { PRINT_ALLOC_ERR; goto restore; }

    /* Read geometry and BZERO/BSCALE */
    f->naxes[2] = 1;
    fits_get_img_param(fptr, 3, &f->bitpix, &f->naxis, f->naxes, &status);
    if (status) { report_fits_error(status); free(f); f = NULL; goto restore; }

    manage_bitpix(fptr, &f->bitpix, &f->orig_bitpix);

    if (f->naxis == 2 && f->naxes[2] == 0)
        f->naxes[2] = 1;

    f->rx = (unsigned int)f->naxes[0];
    f->ry = (unsigned int)f->naxes[1];
    f->fptr = fptr;  /* temporary borrow for read_fits_with_convert */

    if (read_fits_with_convert(f, "FLIS layer", FALSE)) {
        clearfits(f);
        free(f);
        f = NULL;
        goto restore;
    }
    f->fptr = NULL;
    f->top_down = FALSE;

    /* Read FITS keywords (exposure, WCS, etc.) into the layer's fits struct.
     * This preserves per-layer science metadata for round-tripping. */
    f->fptr = fptr;
    read_fits_header(f);
    f->fptr = NULL;

restore:
    status = 0;
    fits_movabs_hdu(fptr, orig_hdu, NULL, &status);
    return f;
}

/*
 * Load a mask from the HDU at hdu_index into a newly allocated layermask_t.
 * The file's active HDU is restored on return.
 * @expect_float: TRUE for processing masks (FLOAT_IMG), FALSE for layer masks (BYTE_IMG).
 */
/* Read the MASK_ACT logical keyword from the given HDU index.
 * Returns TRUE (active) if the keyword is absent or unreadable. */
static gboolean read_mask_act_from_hdu(fitsfile *fptr, int hdu_index) {
    int hdu_status = 0;
    int orig_hdu;
    fits_get_hdu_num(fptr, &orig_hdu);
    int act_val = 1; /* default: active */
    if (!fits_movabs_hdu(fptr, hdu_index, NULL, &hdu_status)) {
        fits_read_key(fptr, TLOGICAL, "MASK_ACT", &act_val, NULL, &hdu_status);
        hdu_status = 0;
        fits_movabs_hdu(fptr, orig_hdu, NULL, &hdu_status);
    }
    return (act_val != 0);
}

static layermask_t *load_mask_from_hdu(fitsfile *fptr, int hdu_index,
                                       gboolean expect_float) {
    int status = 0;
    int orig_hdu;
    fits_get_hdu_num(fptr, &orig_hdu);

    if (fits_movabs_hdu(fptr, hdu_index, NULL, &status)) {
        report_fits_error(status);
        return NULL;
    }

    int bitpix, naxis;
    long naxes[2] = { 0, 0 };
    fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status);
    if (status || naxis != 2) {
        report_fits_error(status);
        fits_movabs_hdu(fptr, orig_hdu, NULL, &status);
        return NULL;
    }

    layermask_t *mask = calloc(1, sizeof(layermask_t));
    if (!mask) { PRINT_ALLOC_ERR; goto restore_mask; }

    mask->w = (size_t)naxes[0];
    mask->h = (size_t)naxes[1];
    mask->bitpix = expect_float ? 32 : 8;
    size_t nbytes = mask->w * mask->h * (mask->bitpix / 8);

    mask->data = malloc(nbytes);
    if (!mask->data) {
        PRINT_ALLOC_ERR;
        free(mask);
        mask = NULL;
        goto restore_mask;
    }

    int zero = 0;
    int cftype = expect_float ? TFLOAT : TBYTE;
    long npix = (long)(mask->w * mask->h);
    fits_read_img(fptr, cftype, 1, npix, &zero, mask->data, &zero, &status);
    if (status) {
        report_fits_error(status);
        layermask_free(mask);
        mask = NULL;
    }

restore_mask:
    status = 0;
    fits_movabs_hdu(fptr, orig_hdu, NULL, &status);
    return mask;
}

/* ===================================================================== */
/* Public API — lifecycle helpers                                        */
/* ===================================================================== */

void layermask_free(layermask_t *mask) {
    if (!mask) return;
    free(mask->data);
    free(mask);
}

flis_layer_t *flis_layer_new(fits *fit, const gchar *name) {
    flis_layer_t *layer = g_new0(flis_layer_t, 1);
    if (!layer) return NULL;

    layer->fit          = fit;
    layer->layer_name   = g_strdup(name ? name : "Layer");
    layer->blend_mode   = FLIS_BLEND_NORMAL;
    layer->opacity      = 1.0f;
    layer->visible      = TRUE;
    layer->has_tint     = FALSE;
    layer->layer_tint   = (flis_tint_t){ 1.0, 1.0, 1.0 };
    layer->lmask_active = TRUE;
    /* item_id and layer_order must be set by the caller */
    return layer;
}

void flis_layer_free(flis_layer_t *layer) {
    if (!layer) return;
    if (layer->fit) {
        clearfits(layer->fit);
        free(layer->fit);
    }
    layermask_free(layer->lmask);
    g_free(layer->layer_name);
    g_free(layer->created);
    g_free(layer->modified);
    g_free(layer);
}

void flis_free_layers(single *uniq) {
    if (!uniq) return;
    g_slist_free_full(uniq->layers, (GDestroyNotify)flis_layer_free);
    uniq->layers = NULL;
    uniq->fit    = NULL;
    uniq->chans  = 0;
    uniq->active_layer = 0;
    flis_free_groups(uniq);
}

void uniq_set_active_layer(single *uniq, gint index) {
    flis_layer_t *layer = (flis_layer_t *)g_slist_nth_data(uniq->layers, index);
    g_return_if_fail(layer != NULL && layer->fit != NULL);
    uniq->active_layer = index;
    uniq->fit   = layer->fit;
    uniq->chans = (layer->fit->naxes[2] > 0) ? (int)layer->fit->naxes[2] : 1;
    gfit        = layer->fit;
}

/* ===================================================================== */
/* Public API — save_flis / load_flis                                    */
/* ===================================================================== */

/*
 * Comparator for g_slist_sort / g_slist_insert_sorted.
 * Sorts flis_layer_t* entries ascending by layer_order.
 * Used by both save_flis (for consistent HDU layout) and load_flis
 * (to maintain the sorted invariant of com.uniq->layers).
 */
static gint layer_order_cmp(gconstpointer a, gconstpointer b) {
    const flis_layer_t *la = (const flis_layer_t *)a;
    const flis_layer_t *lb = (const flis_layer_t *)b;
    return la->layer_order - lb->layer_order;
}

int save_flis(const gchar *filename) {
    if (!com.uniq || !com.uniq->layers) {
        siril_log_message(_("FLIS save: no layers to save\n"));
        return 1;
    }

    gchar *outpath = g_strdup(filename);

    if (g_unlink(outpath))
        siril_debug_print("g_unlink() failed (may not exist)\n");

    int status = 0;
    fitsfile *fptr = NULL;
    if (siril_fits_create_diskfile(&fptr, outpath, &status)) {
        report_fits_error(status);
        g_free(outpath);
        return 1;
    }

    /* Sort layers by layer_order (ascending) for consistent HDU layout */
    GSList *sorted = g_slist_copy(com.uniq->layers);
    sorted = g_slist_sort(sorted, (GCompareFunc)layer_order_cmp);

    /* Determine canvas dimensions from the base layer (lowest layer_order) */
    flis_layer_t *base = (flis_layer_t *)sorted->data;
    long canvas_w = base->fit ? base->fit->rx : 0;
    long canvas_h = base->fit ? base->fit->ry : 0;

    /* Determine whether an ICC profile should be written.  The profile
     * is image-level state on com.uniq, not per-layer. */
    cmsHPROFILE save_profile = (com.uniq) ? com.uniq->icc_profile : NULL;
    gboolean    save_managed = (com.uniq) ? com.uniq->color_managed : FALSE;
    gboolean icc_present = (save_profile != NULL && save_managed
                            && com.pref.fits_save_icc);

    /* ----------------------------------------------------------------
     * HDU 1 (index 0): primary thumbnail HDU
     * ---------------------------------------------------------------- */
    if (write_thumbnail_hdu(fptr, base->fit, canvas_w, canvas_h, icc_present)) {
        siril_log_error(_("FLIS: failed to write thumbnail HDU\n"));
        fits_close_file(fptr, &status);
        g_slist_free(sorted);
        g_free(outpath);
        return 1;
    }

    /* ----------------------------------------------------------------
     * HDU 2 (optional): ICC profile
     * ---------------------------------------------------------------- */
    if (icc_present) {
        if (write_icc_profile_to_fptr(fptr, save_profile)) {
            siril_log_warning(_("FLIS: warning — ICC profile write failed, continuing\n"));
            icc_present = FALSE;
        }
    }

    /* ----------------------------------------------------------------
     * HDU 2 or 3: FLIS_META binary table.
     * first_data_hdu = thumbnail(1) + optional_icc + meta_table + 1
     * All indices are 1-based (CFITSIO convention).
     * ---------------------------------------------------------------- */
    int first_data_hdu = icc_present ? 4 : 3;

    if (write_flis_meta_hdu(fptr, sorted, first_data_hdu)) {
        siril_log_error(_("FLIS: failed to write FLIS_META table\n"));
        fits_close_file(fptr, &status);
        g_slist_free(sorted);
        g_free(outpath);
        return 1;
    }

    /* ----------------------------------------------------------------
     * Data HDUs: layers and masks, in layer_order sequence.
     * Apply the same compression settings as standard FITS saves.
     * ---------------------------------------------------------------- */
    if (com.pref.comp.fits_enabled) {
        /* siril_fits_compress() needs a fits* with a valid fptr.
         * We use a temporary wrapper — it only reads fptr. */
        fits tmp_f = { .fptr = fptr };
        siril_fits_compress(&tmp_f);
    }

    int saved_ok = 1;
    for (GSList *l = sorted; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;

        /* --- Layer pixel data HDU --- */
        if (write_layer_hdu(fptr, lay)) {
            siril_log_error(_("FLIS: failed writing layer '%s'\n"),
                                    lay->layer_name ? lay->layer_name : "?");
            saved_ok = 0;
            break;
        }

        /* --- Layer mask HDU (LMASK) --- */
        if (lay->lmask) {
            gchar *lname = g_strdup_printf("%s Layer Mask",
                           lay->layer_name ? lay->layer_name : "Layer");
            int lmask_id = lay->item_id + 10000;
            if (write_mask_hdu(fptr, lay->lmask, lmask_id,
                               FLIS_TYPE_LMASK, lname,
                               BYTE_IMG, TBYTE, lay->lmask_active)) {
                siril_log_warning(_("FLIS: failed writing layer mask for '%s'\n"), lay->layer_name ? lay->layer_name : "?");
            }
            g_free(lname);
        }

        /* --- Processing mask HDU (MASK) --- */
        if (lay->fit && lay->fit->mask) {
            /* Processing masks are always written as float32 [0,1] in FLIS.
             * mask_t.bitpix follows FITS convention (8/16/32 = bits per pixel),
             * so convert uint8 or uint16 masks to float before writing. */
            mask_t *pmask = lay->fit->mask;
            gchar *mname = g_strdup_printf("%s Processing Mask",
                           lay->layer_name ? lay->layer_name : "Layer");
            int pmask_id = lay->item_id + 20000;
            size_t npix = lay->fit->rx * lay->fit->ry;

            float *float_data = NULL;
            void  *write_data = pmask->data;
            if (pmask->bitpix != 32) {
                float_data = malloc(npix * sizeof(float));
                if (float_data) {
                    if (pmask->bitpix == 8) {
                        const uint8_t *src = (const uint8_t *)pmask->data;
                        for (size_t i = 0; i < npix; i++)
                            float_data[i] = src[i] * (1.0f / 255.0f);
                    } else { /* 16 */
                        const uint16_t *src = (const uint16_t *)pmask->data;
                        for (size_t i = 0; i < npix; i++)
                            float_data[i] = src[i] * (1.0f / 65535.0f);
                    }
                    write_data = float_data;
                }
            }

            layermask_t tmp_lm = {
                .w      = lay->fit->rx,
                .h      = lay->fit->ry,
                .bitpix = 32,
                .data   = write_data
            };
            if (write_mask_hdu(fptr, &tmp_lm, pmask_id,
                               FLIS_TYPE_MASK, mname,
                               FLOAT_IMG, TFLOAT, lay->fit->mask_active)) {
                siril_log_warning(_("FLIS: failed writing processing mask for '%s'\n"), lay->layer_name ? lay->layer_name : "?");
            }
            free(float_data);
            g_free(mname);
        }
    }

    status = 0;
    fits_close_file(fptr, &status);
    g_slist_free(sorted);

    if (!status && saved_ok) {
        siril_log_message(_("FLIS: saved %d layer(s) to %s\n"),
                          g_slist_length(com.uniq->layers), outpath);
    } else if (status) {
        report_fits_error(status);
    }

    g_free(outpath);
    return (!saved_ok || status) ? 1 : 0;
}

/* ===================================================================== */
/* Public API — load_flis                                                */
/* ===================================================================== */

int load_flis(const gchar *filename) {
    int status = 0;
    fitsfile *fptr = NULL;

    /* Open the file without jumping to first image (we need HDU 1 = thumbnail) */
    fits_open_diskfile(&fptr, filename, READONLY, &status);
    if (status) {
        report_fits_error(status);
        return 1;
    }

    /* ----------------------------------------------------------------
     * HDU 1: validate FLIS identification
     * ---------------------------------------------------------------- */
    fits_movabs_hdu(fptr, 1, NULL, &status);
    if (status) { report_fits_error(status); fits_close_file(fptr, &status); return 1; }

    int is_flis = 0;
    fits_read_key(fptr, TLOGICAL, "FLIS", &is_flis, NULL, &status);
    if (status || !is_flis) {
        siril_log_error(_("FLIS: file %s is not a valid FLIS file (no FLIS=T keyword)\n"), filename);
        fits_close_file(fptr, &status);
        return 1;
    }
    status = 0;

    char ver[FLEN_VALUE] = { 0 };
    fits_read_key(fptr, TSTRING, "FLISVER", ver, NULL, &status);
    if (!status && strcmp(ver, FLIS_VERSION_STRING)) {
        siril_log_message(_("FLIS: file version %s differs from implementation version %s\n"),
                          ver, FLIS_VERSION_STRING);
    }
    status = 0;

    /* Read canvas dimensions */
    int canvas_w = 0, canvas_h = 0;
    fits_read_key(fptr, TINT, "FLISSIZX", &canvas_w, NULL, &status);
    fits_read_key(fptr, TINT, "FLISSIZY", &canvas_h, NULL, &status);
    status = 0;

    /* Check ICC flag */
    int icc_present = 0;
    fits_read_key(fptr, TLOGICAL, "FLISICC", &icc_present, NULL, &status);
    status = 0;

    /* ----------------------------------------------------------------
     * Locate the FLIS_META table.  It is at HDU 2 (no ICC) or HDU 3.
     * Verify by EXTNAME in case HDU layout shifted.
     * ---------------------------------------------------------------- */
    int meta_hdu = icc_present ? 3 : 2;
    int nhdus = 0;
    fits_get_num_hdus(fptr, &nhdus, &status);

    /* Verify expected position; fall back to scanning if wrong */
    fits_movabs_hdu(fptr, meta_hdu, NULL, &status);
    if (status) {
        siril_log_error(_("FLIS: cannot locate FLIS_META table\n"));
        fits_close_file(fptr, &status);
        return 1;
    }
    status = 0;

    char extname[FLEN_VALUE] = { 0 };
    fits_read_key(fptr, TSTRING, "EXTNAME", extname, NULL, &status);
    status = 0;
    if (g_ascii_strcasecmp(extname, "FLIS_META")) {
        /* Scan remaining HDUs for FLIS_META */
        gboolean found = FALSE;
        for (int h = 2; h <= nhdus; h++) {
            fits_movabs_hdu(fptr, h, NULL, &status); status = 0;
            fits_read_key(fptr, TSTRING, "EXTNAME", extname, NULL, &status); status = 0;
            if (!g_ascii_strcasecmp(extname, "FLIS_META")) {
                meta_hdu = h;
                found = TRUE;
                break;
            }
        }
        if (!found) {
            siril_log_error(_("FLIS: FLIS_META table not found\n"));
            fits_close_file(fptr, &status);
            return 1;
        }
        fits_movabs_hdu(fptr, meta_hdu, NULL, &status);
    }

    /* ----------------------------------------------------------------
     * Read the FLIS_META table
     * ---------------------------------------------------------------- */
    long nrows = 0;
    flis_meta_row_t *meta_rows = read_flis_meta(fptr, &nrows);
    if (!meta_rows) {
        fits_close_file(fptr, &status);
        return 1;
    }

    /* ----------------------------------------------------------------
     * Read ICC profile if present
     * ---------------------------------------------------------------- */
    cmsHPROFILE file_icc = NULL;
    if (icc_present) {
        fits_movabs_hdu(fptr, 2, NULL, &status); status = 0;
        file_icc = read_icc_profile_from_fptr(fptr);
    }

    /* ----------------------------------------------------------------
     * Build the layer list from the metadata table rows.
     *
     * Pass 1: create flis_layer_t for every LAYER row, loading pixel data.
     * Pass 2: attach LMASK and MASK items to their parent layers.
     * ---------------------------------------------------------------- */

    /* Map from item_id -> flis_layer_t* for parent lookup */
    GHashTable *id_map = g_hash_table_new(g_direct_hash, g_direct_equal);
    GSList *layers = NULL;
    GSList *groups = NULL;

    /* Pass 0: GROUP rows — build group objects first so layers can reference them */
    for (long r = 0; r < nrows; r++) {
        flis_meta_row_t *row = &meta_rows[r];
        if (g_ascii_strcasecmp(row->item_type, FLIS_TYPE_GROUP))
            continue;
        flis_group_t *grp = flis_group_new(row->layer_name);
        if (!grp) continue;
        grp->item_id    = row->item_id;
        grp->opacity    = row->opacity;
        grp->visible    = row->visible;
        grp->blend_mode = str_to_blend_mode(row->blend_mode);
        /* Parse COLLAPSED from metadata string */
        if (row->metadata && strstr(row->metadata, "COLLAPSED=1"))
            grp->collapsed = TRUE;
        /* Parse timestamps from metadata: CREATED=...; MODIFIED=... */
        {
            const char *p = row->metadata ? strstr(row->metadata, "CREATED=") : NULL;
            if (p) {
                p += 8;
                const char *end = strchr(p, ';');
                if (end && end > p)
                    grp->created = g_strndup(p, end - p);
                else if (*p)
                    grp->created = g_strdup(p);
            }
            const char *q = row->metadata ? strstr(row->metadata, "MODIFIED=") : NULL;
            if (q) {
                q += 9;
                const char *end = strchr(q, ';');
                if (end && end > q)
                    grp->modified = g_strndup(q, end - q);
                else if (*q)
                    grp->modified = g_strdup(q);
            }
        }
        groups = g_slist_append(groups, grp);
    }

    /* Pass 1: LAYER rows */
    for (long r = 0; r < nrows; r++) {
        flis_meta_row_t *row = &meta_rows[r];
        if (g_ascii_strcasecmp(row->item_type, FLIS_TYPE_LAYER))
            continue; /* not a LAYER row */

        fits *f = load_layer_from_hdu(fptr, row->hdu_index);
        if (!f) {
            siril_log_warning(
                _("FLIS: failed to load layer '%s' (HDU %d), skipping\n"), row->layer_name, row->hdu_index);
            continue;
        }

        flis_layer_t *layer = flis_layer_new(f, row->layer_name);
        if (!layer) { clearfits(f); free(f); continue; }

        layer->item_id     = row->item_id;
        layer->layer_order = row->layer_order;
        layer->blend_mode  = str_to_blend_mode(row->blend_mode);
        layer->opacity     = row->opacity;
        layer->visible     = row->visible;
        layer->position_x  = row->position_x;
        layer->position_y  = row->position_y;
        parse_metadata_string(row->metadata, layer);
        layer->group_id = row->group_id;

        layers = g_slist_insert_sorted(layers, layer,
                                       (GCompareFunc)layer_order_cmp);
        g_hash_table_insert(id_map,
                            GINT_TO_POINTER(row->item_id), layer);
    }

    /* Assign the file-level ICC profile to the base layer only.
     * The base layer is the first element of the sorted list (lowest
     * layer_order).  For any non-base layer whose HDU happened to carry
     * an embedded ICC profile, warn and discard it — ICC profiles are
     * an image-level concept in FLIS, not a per-layer one. */
    /* Per-fits ICC fields removed.  Any embedded profile that the
     * per-HDU loader picked up was already discarded by the load helper
     * (it only honours profiles when fit == gfit).  ICC state for the
     * FLIS lives on com.uniq and is set further down from file_icc. */

    /* Pass 2: LMASK and MASK rows */
    for (long r = 0; r < nrows; r++) {
        flis_meta_row_t *row = &meta_rows[r];
        gboolean is_lmask = !g_ascii_strcasecmp(row->item_type, FLIS_TYPE_LMASK);
        gboolean is_mask  = !g_ascii_strcasecmp(row->item_type, FLIS_TYPE_MASK);
        if (!is_lmask && !is_mask) continue;

        flis_layer_t *parent = g_hash_table_lookup(id_map,
                                    GINT_TO_POINTER(row->parent_id));
        if (!parent) {
            siril_log_warning(
                _("FLIS: mask '%s' references unknown parent ID %d, skipping\n"), row->layer_name, row->parent_id);
            continue;
        }

        if (is_lmask) {
            layermask_t *lm = load_mask_from_hdu(fptr, row->hdu_index, FALSE);
            if (lm) {
                layermask_free(parent->lmask); /* replace any existing */
                parent->lmask        = lm;
                parent->lmask_active = read_mask_act_from_hdu(fptr, row->hdu_index);
            }
        } else { /* MASK — processing mask, float */
            layermask_t *pm = load_mask_from_hdu(fptr, row->hdu_index, TRUE);
            if (pm && parent->fit) {
                /* Wrap into mask_t (Siril's existing processing mask struct) */
                if (parent->fit->mask) {
                    free_mask(parent->fit->mask);
                    parent->fit->mask = NULL;
                }
                mask_t *siril_mask = calloc(1, sizeof(mask_t));
                if (siril_mask) {
                    siril_mask->bitpix   = 32;
                    siril_mask->data     = pm->data;
                    pm->data             = NULL; /* transfer ownership */
                    parent->fit->mask    = siril_mask;
                }
                parent->fit->mask_active = read_mask_act_from_hdu(fptr, row->hdu_index);
                layermask_free(pm);
            } else {
                layermask_free(pm);
            }
        }
    }

    g_hash_table_destroy(id_map);
    for (long r = 0; r < nrows; r++)
        g_free(meta_rows[r].metadata);
    free(meta_rows);
    if (file_icc) cmsCloseProfile(file_icc);
    fits_close_file(fptr, &status);

    if (!layers) {
        siril_log_error(_("FLIS: no layers loaded from %s\n"), filename);
        return 1;
    }

    /* ----------------------------------------------------------------
     * Populate com.uniq
     * ---------------------------------------------------------------- */
    if (!com.uniq)
        com.uniq = g_new0(single, 1);
    else
        flis_free_layers(com.uniq);

    g_free(com.uniq->filename);
    com.uniq->filename  = g_strdup(filename);
    com.uniq->fileexist = TRUE;
    com.uniq->layers    = layers;
    com.uniq->groups    = groups;
    com.uniq->next_item_id = 1; /* will be set properly below */

    /* ICC profile (image-level — see flis_get_profiled_fit comment) lives
     * on com.uniq.  The base layer fit will continue to mirror it for the
     * duration of the migration via current_image_set_icc_profile so that
     * legacy fit->icc_profile readers still work. */
    if (com.uniq->icc_profile) {
        cmsCloseProfile(com.uniq->icc_profile);
        com.uniq->icc_profile = NULL;
    }
    if (file_icc) {
        com.uniq->icc_profile  = copyICCProfile(file_icc);
        com.uniq->color_managed = TRUE;
    } else {
        com.uniq->color_managed = FALSE;
    }

    /* Find the highest item_id to seed next_item_id */
    gint max_id = 0;
    for (GSList *l = layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (lay->item_id > max_id) max_id = lay->item_id;
    }
    for (GSList *g = groups; g; g = g->next) {
        flis_group_t *grp = (flis_group_t *)g->data;
        if (grp->item_id > max_id) max_id = grp->item_id;
    }
    com.uniq->next_item_id = max_id + 1;

    /* Set active layer to the base (index 0 = lowest layer_order after sort) */
    uniq_set_active_layer(com.uniq, 0);

    siril_log_message(_("FLIS: loaded %d layer(s) from %s (%dx%d canvas)\n"),
                      g_slist_length(layers), filename, canvas_w, canvas_h);
    return 0;
}

/* ===================================================================== */
/* Internal helpers shared by layer management functions                 */
/* ===================================================================== */

/*
 * Emit an ISO 8601 UTC timestamp string into a heap-allocated gchar*.
 * Returns a g_malloc'd string; caller must g_free().
 */
static gchar *flis_now_iso8601(void) {
    GDateTime *now = g_date_time_new_now_utc();
    gchar *s = g_date_time_format_iso8601(now);
    g_date_time_unref(now);
    return s;
}

/*
 * Return TRUE and log a warning if the layer is locked.
 * @op: short operation name for the log message.
 */
static gboolean flis_check_locked(const flis_layer_t *layer, const char *op) {
    if (layer->locked) {
        siril_log_warning(
            _("FLIS: cannot %s — layer '%s' is locked\n"),
            op, layer->layer_name ? layer->layer_name : "?");
        return TRUE;
    }
    return FALSE;
}

/*
 * Re-sort com.uniq->layers by layer_order and re-sync the active_layer
 * index so it still points at the same flis_layer_t* after sorting.
 * Must be called after any operation that modifies layer_order values.
 */
static void flis_resort_layers(flis_layer_t *keep_active) {
    com.uniq->layers = g_slist_sort(com.uniq->layers,
                                    (GCompareFunc)layer_order_cmp);
    if (keep_active) {
        gint idx = flis_layer_get_index(keep_active);
        if (idx >= 0)
            uniq_set_active_layer(com.uniq, idx);
    }
}

/* Public wrapper used by the undo mechanism to re-sort after restoring
 * layer_order values directly. */
void flis_sort_layer_stack(void) {
    if (!com.uniq || !com.uniq->layers) return;
    flis_layer_t *active = flis_active_layer();
    flis_resort_layers(active);
}

/*
 * Return the maximum layer_order value currently in the stack, or 0 if
 * the stack is empty.  Used when appending a new layer at the top.
 */
static gint flis_max_layer_order(void) {
    gint max = 0;
    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (lay->layer_order > max)
            max = lay->layer_order;
    }
    return max;
}

/* ===================================================================== */
/* Layer lookup helpers                                                  */
/* ===================================================================== */

flis_layer_t *flis_active_layer(void) {
    if (!com.uniq || !com.uniq->layers)
        return NULL;
    return (flis_layer_t *)g_slist_nth_data(com.uniq->layers,
                                             com.uniq->active_layer);
}

fits *flis_active_layer_fit(void) {
    flis_layer_t *lay = flis_active_layer();
    return lay ? lay->fit : NULL;
}

fits *flis_get_profiled_fit(void) {
    if (!is_current_image_flis() || !com.uniq || !com.uniq->layers)
        return gfit;
    flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
    return (base && base->fit) ? base->fit : gfit;
}

guint flis_composite_naxes2(void) {
    /* The FLIS display composite is always built as RGB (3 channels),
     * regardless of whether the individual layers are mono.  Return 3 for
     * any FLIS image so that ICC profile channel-count checks are performed
     * against the composite rather than the (potentially mono) base layer. */
    if (is_current_image_flis())
        return 3;
    return gfit->naxes[2];
}

gboolean flis_composite_is_chromatic(void) {
    /* "Chromatic" = the displayed composite is not equivalent to a mono
     * image broadcast to three channels.  An all-mono, no-tint (or only
     * greyscale-tint) FLIS produces R=G=B everywhere and so should still
     * present mono-style GUI controls (single-channel histogram, mono
     * vport tab) even though the composite buffer happens to be 3-channel.
     *
     * Any RGB layer brings chromaticity unconditionally.  A mono layer
     * with LAYER_COLOR tint contributes chroma only when the tint vector
     * is not greyscale (i.e. R != G or G != B) — a (0.5, 0.5, 0.5) tint
     * is just a brightness scalar. */
    if (!is_current_image_flis() || !com.uniq) return FALSE;
    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (!lay || !lay->fit) continue;
        if (lay->fit->naxes[2] >= 3)
            return TRUE;
        if (lay->has_tint &&
            (lay->layer_tint.r != lay->layer_tint.g ||
             lay->layer_tint.g != lay->layer_tint.b))
            return TRUE;
    }
    return FALSE;
}

/* Luminance-weighted collapse of three planar channels into one.
 * Uses Rec. 709 coefficients (0.2126 R + 0.7152 G + 0.0722 B). */
static void collapse_rgb_to_mono_float(float *r, float *g, float *b,
                                       float *out, size_t npixels) {
    for (size_t i = 0; i < npixels; i++)
        out[i] = 0.2126f * r[i] + 0.7152f * g[i] + 0.0722f * b[i];
}

static void collapse_rgb_to_mono_word(WORD *r, WORD *g, WORD *b,
                                      WORD *out, size_t npixels) {
    for (size_t i = 0; i < npixels; i++)
        out[i] = (WORD)(0.2126f * r[i] + 0.7152f * g[i] + 0.0722f * b[i] + 0.5f);
}

void flis_convert_layers_icc(cmsHPROFILE old_profile, cmsHPROFILE new_profile) {
    if (!is_current_image_flis() || !com.uniq || !com.uniq->layers) return;
    if (!old_profile || !new_profile) return;

    gboolean threaded = !processing_is_job_active();
    cmsUInt32Number intent = com.pref.icc.processing_intent;

    cmsUInt32Number rgb_float_type = get_planar_formatter_type(cmsSigRgbData, DATA_FLOAT, FALSE);
    cmsUInt32Number rgb_word_type  = get_planar_formatter_type(cmsSigRgbData, DATA_USHORT, FALSE);

    /* Build a single RGB→RGB transform, reused for every layer. */
    cmsHTRANSFORM xform_float = cmsCreateTransformTHR(
        threaded ? com.icc.context_threaded : com.icc.context_single,
        old_profile, rgb_float_type, new_profile, rgb_float_type,
        intent, com.icc.rendering_flags);
    cmsHTRANSFORM xform_word = cmsCreateTransformTHR(
        threaded ? com.icc.context_threaded : com.icc.context_single,
        old_profile, rgb_word_type, new_profile, rgb_word_type,
        intent, com.icc.rendering_flags);
    /* Packed interleaved transform for tint vectors (3 floats, single pixel).
     * The tint is defined as linear-light values, but we feed it through the
     * same ICC transform as the image data.  For values at 0.0 and 1.0 this
     * is exact; for intermediate values the TRC error is small and acceptable
     * for typical narrowband tint ratios. */
    cmsHTRANSFORM xform_tint = cmsCreateTransformTHR(
        threaded ? com.icc.context_threaded : com.icc.context_single,
        old_profile, TYPE_RGB_FLT, new_profile, TYPE_RGB_FLT,
        intent, com.icc.rendering_flags);

    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (!lay || !lay->fit) continue;
        fits *fit = lay->fit;
        size_t npixels = (size_t)fit->rx * fit->ry;
        cmsHTRANSFORM xform = (fit->type == DATA_FLOAT) ? xform_float : xform_word;
        if (!xform) continue;

        cmsUInt32Number dsize       = (fit->type == DATA_FLOAT) ? sizeof(float) : sizeof(WORD);
        cmsUInt32Number bytesperline  = fit->rx * dsize;
        cmsUInt32Number bytesperplane = npixels * dsize;

        if (fit->naxes[2] >= 3) {
            /* RGB layer: transform in-place. */
            void *data = (fit->type == DATA_FLOAT) ? (void *)fit->fdata : (void *)fit->data;
            cmsDoTransformLineStride(xform, data, data, fit->rx, fit->ry,
                                     bytesperline, bytesperline,
                                     bytesperplane, bytesperplane);
        } else {
            /* Mono layer: broadcast to 3 RGB planes, transform, collapse. */
            void *tmp = g_malloc(npixels * 3 * dsize);
            if (fit->type == DATA_FLOAT) {
                float *mono  = fit->fdata;
                float *tr    = (float *)tmp;
                float *tg    = tr + npixels;
                float *tb    = tg + npixels;
                for (size_t i = 0; i < npixels; i++)
                    tr[i] = tg[i] = tb[i] = mono[i];
                cmsDoTransformLineStride(xform, tmp, tmp, fit->rx, fit->ry,
                                         bytesperline, bytesperline,
                                         bytesperplane, bytesperplane);
                collapse_rgb_to_mono_float(tr, tg, tb, mono, npixels);
            } else {
                WORD *mono = fit->data;
                WORD *tr   = (WORD *)tmp;
                WORD *tg   = tr + npixels;
                WORD *tb   = tg + npixels;
                for (size_t i = 0; i < npixels; i++)
                    tr[i] = tg[i] = tb[i] = mono[i];
                cmsDoTransformLineStride(xform, tmp, tmp, fit->rx, fit->ry,
                                         bytesperline, bytesperline,
                                         bytesperplane, bytesperplane);
                collapse_rgb_to_mono_word(tr, tg, tb, mono, npixels);
            }
            g_free(tmp);

            /* Convert the tint to the new colorspace.  The tint encodes the
             * colour direction of the layer in the source space; leaving it
             * unchanged after an ICC conversion would produce the wrong hue
             * (e.g. an sRGB-green tint rendered as Rec2020-green is far
             * outside sRGB gamut and appears very wrong on a standard display).
             *
             * The relative-colorimetric transform can produce out-of-gamut
             * values when the source-space tint sits outside the destination
             * gamut (e.g. Rec2020 green (0,1,0) → sRGB ≈ (-0.59, 1.13, -0.10)).
             * Those values break the SCREEN/blend math downstream — the
             * (1 - source) factor inverts sign and the composite picks up a
             * cyan/teal cast instead of the green the user expects.  Clamp
             * into [0, 1]: a no-op for in-gamut tints, and for out-of-gamut
             * tints it lands on the same primary the monitor would have
             * displayed via the proofing transform anyway. */
            if (lay->has_tint && xform_tint) {
                float tint_in[3]  = { (float)lay->layer_tint.r,
                                      (float)lay->layer_tint.g,
                                      (float)lay->layer_tint.b };
                float tint_out[3] = { 0.f, 0.f, 0.f };
                cmsDoTransform(xform_tint, tint_in, tint_out, 1);
                lay->layer_tint.r = (double)CLAMP(tint_out[0], 0.0f, 1.0f);
                lay->layer_tint.g = (double)CLAMP(tint_out[1], 0.0f, 1.0f);
                lay->layer_tint.b = (double)CLAMP(tint_out[2], 0.0f, 1.0f);
            }
        }
    }

    if (xform_float) cmsDeleteTransform(xform_float);
    if (xform_word)  cmsDeleteTransform(xform_word);
    if (xform_tint)  cmsDeleteTransform(xform_tint);

    /* Profile bookkeeping is now on com.uniq, not per-layer.  The caller
     * (icc_convert_to_hook) calls siril_colorspace_transform after this
     * which lands the new profile via current_image_set_icc_profile. */
}

/* flis_update_layer_offset_after_crop:
 *
 * Must be called immediately after a crop operation has been applied to gfit
 * (i.e. after the pixel data has been trimmed and gfit->rx/ry updated) when
 * a FLIS image is loaded.
 *
 * @sel_x, @sel_y: top-left of the crop selection in canvas coordinates
 *                 (i.e. com.selection.x and com.selection.y before the crop)
 *
 * Two cases:
 *
 * 1. Base layer cropped: the canvas itself shrinks.  All non-base layers
 *    have their offsets reduced by (sel_x, sel_y) so they remain correctly
 *    positioned relative to the new canvas origin.
 *
 * 2. Non-base layer cropped: the layer's pixel data is now the selection
 *    area.  Its position on the canvas is the selection top-left, so set
 *    position_x = sel_x, position_y = sel_y.
 */
void flis_update_layer_offset_after_crop(gint sel_x, gint sel_y) {
    if (!is_current_image_flis() || !com.uniq || !com.uniq->layers) return;

    flis_layer_t *active = flis_active_layer();
    if (!active) return;

    flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;

    if (active == base) {
        /* Base layer cropped: adjust all non-base layer offsets */
        for (GSList *l = com.uniq->layers->next; l; l = l->next) {
            flis_layer_t *lay = (flis_layer_t *)l->data;
            if (!lay) continue;
            lay->position_x -= sel_x;
            lay->position_y -= sel_y;
        }
    } else {
        /* Non-base layer cropped: the data now sits at (sel_x, sel_y) */
        active->position_x = sel_x;
        active->position_y = sel_y;
    }

    gui_iface.flis_invalidate_composite();
    siril_debug_print("FLIS: layer offsets updated after crop (sel %d,%d)\n",
                      sel_x, sel_y);
}

/* flis_update_layer_offset_after_resize:
 *
 * Must be called immediately after a resize operation has been applied to gfit
 * (i.e. after the pixel data has been scaled and gfit->rx/ry updated) when
 * a FLIS image is loaded.
 *
 * @old_rx, @old_ry: canvas dimensions before the resize
 * @new_rx, @new_ry: canvas dimensions after the resize
 *
 * Two cases:
 *
 * 1. Base layer resized: the canvas itself scales.  All non-base layers have
 *    their centres scaled proportionally so they remain at the same relative
 *    position on the canvas.
 *
 * 2. Non-base layer resized: the layer's pixel data changes size but its
 *    centre on the canvas stays fixed.  position_x/y are adjusted by half
 *    the dimension delta so the centre is preserved.
 */
void flis_update_layer_offset_after_resize(gint old_rx, gint old_ry,
                                           gint new_rx, gint new_ry) {
    if (!is_current_image_flis() || !com.uniq || !com.uniq->layers) return;

    flis_layer_t *active = flis_active_layer();
    if (!active) return;

    flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;

    if (active == base) {
        /* Base layer resized: scale all non-base layer centres proportionally */
        for (GSList *l = com.uniq->layers->next; l; l = l->next) {
            flis_layer_t *lay = (flis_layer_t *)l->data;
            if (!lay) continue;
            double cx = lay->position_x + lay->fit->rx / 2.0;
            double cy = lay->position_y + lay->fit->ry / 2.0;
            cx *= (double)new_rx / old_rx;
            cy *= (double)new_ry / old_ry;
            lay->position_x = (gint)round(cx - lay->fit->rx / 2.0);
            lay->position_y = (gint)round(cy - lay->fit->ry / 2.0);
        }
    } else {
        /* Non-base layer resized: preserve its centre on the canvas */
        active->position_x += (old_rx - new_rx) / 2;
        active->position_y += (old_ry - new_ry) / 2;
    }

    gui_iface.flis_invalidate_composite();
    siril_debug_print("FLIS: layer offsets updated after resize (%dx%d -> %dx%d)\n",
                      old_rx, old_ry, new_rx, new_ry);
}

/* flis_update_layer_offset_after_rotate:
 *
 * Must be called immediately after a rotation operation has been applied to
 * gfit (i.e. after the pixel data has been transformed and gfit->rx/ry
 * updated) when a FLIS image is loaded.
 *
 * @old_rx, @old_ry: canvas dimensions before the rotation
 * @new_rx, @new_ry: canvas dimensions after the rotation
 * @angle:           rotation angle in degrees (positive = CCW in Siril convention)
 *
 * Two cases:
 *
 * 1. Base layer rotated: the canvas itself rotates.  All non-base layers have
 *    their centres rotated by the same angle around the old canvas centre and
 *    then mapped into the new canvas.  (Their pixel data is not rotated.)
 *
 * 2. Non-base layer rotated: the layer's pixel data is rotated but its centre
 *    on the canvas stays fixed.  position_x/y are adjusted by half the
 *    dimension delta so the centre is preserved.
 */
void flis_update_layer_offset_after_rotate(gint old_rx, gint old_ry,
                                           gint new_rx, gint new_ry,
                                           double angle) {
    if (!is_current_image_flis() || !com.uniq || !com.uniq->layers) return;

    flis_layer_t *active = flis_active_layer();
    if (!active) return;

    flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;

    if (active == base) {
        /* Base layer rotated: rotate each non-base layer's centre position
         * by the same angle around the old canvas centre, then map to the
         * new canvas.  The layer pixel data itself is not rotated. */
        double rad = angle * M_PI / 180.0;
        double cos_a = cos(rad);
        double sin_a = sin(rad);
        double ocx = old_rx / 2.0;
        double ocy = old_ry / 2.0;

        for (GSList *l = com.uniq->layers->next; l; l = l->next) {
            flis_layer_t *lay = (flis_layer_t *)l->data;
            if (!lay) continue;
            double cx = lay->position_x + lay->fit->rx / 2.0;
            double cy = lay->position_y + lay->fit->ry / 2.0;
            double dx = cx - ocx;
            double dy = cy - ocy;
            double new_cx = new_rx / 2.0 + dx * cos_a - dy * sin_a;
            double new_cy = new_ry / 2.0 + dx * sin_a + dy * cos_a;
            lay->position_x = (gint)round(new_cx - lay->fit->rx / 2.0);
            lay->position_y = (gint)round(new_cy - lay->fit->ry / 2.0);
        }
    } else {
        /* Non-base layer rotated: preserve its centre on the canvas */
        active->position_x += (old_rx - new_rx) / 2;
        active->position_y += (old_ry - new_ry) / 2;
    }

    gui_iface.flis_invalidate_composite();
    siril_debug_print("FLIS: layer offsets updated after rotate (%.2f deg, %dx%d -> %dx%d)\n",
                      angle, old_rx, old_ry, new_rx, new_ry);
}

/* flis_update_all_layer_offsets_after_rotate:
 *
 * Group-rotation variant: rotates the centre of EVERY non-base layer around
 * the old canvas centre and maps it onto the new canvas.  Used when an entire
 * layer group (which may include the base layer) has been rotated, so the
 * "active layer" reported by flis_active_layer() is the group node rather than
 * any individual layer.
 *
 * Non-group members that are NOT in the rotated set also need their canvas
 * positions updated here because the canvas itself has rotated.
 */
void flis_update_all_layer_offsets_after_rotate(gint old_rx, gint old_ry,
                                                gint new_rx, gint new_ry,
                                                double angle) {
    if (!is_current_image_flis() || !com.uniq || !com.uniq->layers) return;
    if (old_rx == new_rx && old_ry == new_ry && angle == 0.0) return;

    double rad   = angle * M_PI / 180.0;
    double cos_a = cos(rad);
    double sin_a = sin(rad);
    double ocx   = old_rx / 2.0;
    double ocy   = old_ry / 2.0;

    for (GSList *l = com.uniq->layers->next; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (!lay || !lay->fit) continue;
        double cx  = lay->position_x + lay->fit->rx / 2.0;
        double cy  = lay->position_y + lay->fit->ry / 2.0;
        double dx  = cx - ocx;
        double dy  = cy - ocy;
        double new_cx = new_rx / 2.0 + dx * cos_a - dy * sin_a;
        double new_cy = new_ry / 2.0 + dx * sin_a + dy * cos_a;
        lay->position_x = (gint)round(new_cx - lay->fit->rx / 2.0);
        lay->position_y = (gint)round(new_cy - lay->fit->ry / 2.0);
    }
    gui_iface.flis_invalidate_composite();
    siril_debug_print("FLIS: all layer offsets updated after group rotate (%.2f deg, %dx%d -> %dx%d)\n",
                      angle, old_rx, old_ry, new_rx, new_ry);
}

guint flis_canvas_rx(void) {
    if (!is_current_image_flis() || !com.uniq || !com.uniq->layers) return gfit->rx;
    flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
    return (base && base->fit) ? (guint)base->fit->rx : gfit->rx;
}

guint flis_canvas_ry(void) {
    if (!is_current_image_flis() || !com.uniq || !com.uniq->layers) return gfit->ry;
    flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;
    return (base && base->fit) ? (guint)base->fit->ry : gfit->ry;
}

gboolean flis_canvas_to_pixel_index(gint cx, gint cy_disp, guint canvas_ry,
                                    size_t *out_idx) {
    (void)canvas_ry;  /* unused — kept for API stability */
    gint pos_x = 0, pos_y = 0;

    if (is_current_image_flis() && com.uniq && com.uniq->layers) {
        flis_layer_t *active = flis_active_layer();
        flis_layer_t *base   = (flis_layer_t *)com.uniq->layers->data;
        if (active && active != base) {
            pos_x = active->position_x;
            pos_y = active->position_y;
        }
    }

    gint layer_rx = (gint)gfit->rx;
    gint layer_ry = (gint)gfit->ry;

    /* Coordinate convention (matches flis_compose.c, the display oracle):
     * a layer with position_y / lh occupies display top-down rows
     * [position_y, position_y + lh).  Within that range:
     *
     *   local_display_y = cy_disp - position_y       (0 = layer's top edge)
     *   local_fits_y    = lh - 1 - local_display_y   (flip to FITS bottom-up
     *                                                  within the layer)
     *
     * For the base layer (position_y == 0, lh == canvas_ry) this degenerates
     * to the standard `canvas_ry - 1 - cy_disp` formula.
     *
     * Earlier versions of this function treated position_y as FITS bottom-up
     * (matching the FLIS spec comment) but that disagreed with the compose
     * kernel and would have returned the wrong layer pixel for any sparse
     * upper layer.  The compose's interpretation wins because that's what
     * the user sees on screen. */
    gint local_x        = cx - pos_x;
    gint local_disp_y   = cy_disp - pos_y;
    gint local_fits_y   = layer_ry - 1 - local_disp_y;

    if (local_x < 0 || local_x >= layer_rx ||
        local_disp_y < 0 || local_disp_y >= layer_ry)
        return FALSE;

    *out_idx = (size_t)layer_rx * local_fits_y + local_x;
    return TRUE;
}

int flis_promote_from_gfit(const gchar *name) {
    if (!com.uniq) {
        siril_log_message(_("FLIS: flis_promote_from_gfit — no image loaded\n"));
        return 1;
    }
    if (is_current_image_flis()) {
        siril_log_message(_("FLIS: flis_promote_from_gfit — image is already a FLIS\n"));
        return 1;
    }
    if (!gfit || (!gfit->fdata && !gfit->data)) {
        siril_log_message(_("FLIS: flis_promote_from_gfit — gfit has no pixel data\n"));
        return 1;
    }

    /* Seed the item ID counter.  calloc in create_uniq_from_gfit() leaves
     * it at 0; FLIS item IDs must be >= 1. */
    com.uniq->next_item_id = 1;

    /* Wrap gfit as the base layer.  Ownership of the fits* transfers to
     * the layer; flis_layer_free() will clearfits()+free() it.
     * The caller must ensure gfit is heap-allocated. */
    flis_layer_t *base = flis_layer_new(gfit, name ? name : _("Background"));
    if (!base) return 1;

    base->item_id     = com.uniq->next_item_id++;
    base->layer_order = FLIS_ORDER_STEP;   /* first slot: 10 */
    base->created     = flis_now_iso8601();
    base->modified    = g_strdup(base->created);

    com.uniq->layers = g_slist_prepend(NULL, base);
    uniq_set_active_layer(com.uniq, 0);

    /* com.uniq->filename and fileexist are already correct from
     * create_uniq_from_gfit().  Mark the composite dirty. */
    gui_iface.flis_invalidate_composite();

    siril_debug_print("FLIS: promoted loaded image to single-layer FLIS\n");
    return 0;
}

flis_layer_t *flis_layer_get_by_id(gint item_id) {
    if (!com.uniq) return NULL;
    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (lay->item_id == item_id)
            return lay;
    }
    return NULL;
}

flis_layer_t *flis_layer_get_by_name(const gchar *name) {
    if (!com.uniq || !name) return NULL;
    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (lay->layer_name && !g_strcmp0(lay->layer_name, name))
            return lay;
    }
    return NULL;
}

gint flis_layer_get_index(const flis_layer_t *layer) {
    if (!com.uniq || !layer) return -1;
    gint i = 0;
    for (GSList *l = com.uniq->layers; l; l = l->next, i++) {
        if ((flis_layer_t *)l->data == layer)
            return i;
    }
    return -1;
}

flis_layer_t *flis_layer_get_by_fit(const fits *fit) {
    if (!com.uniq || !com.uniq->layers || !fit) return NULL;
    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (lay && lay->fit == fit) return lay;
    }
    return NULL;
}

gint flis_layer_count(void) {
    if (!com.uniq) return 0;
    return (gint)g_slist_length(com.uniq->layers);
}

/* ===================================================================== */
/* Group management                                                      */
/* ===================================================================== */

flis_group_t *flis_group_new(const gchar *name) {
    flis_group_t *grp = g_new0(flis_group_t, 1);
    if (!grp) { PRINT_ALLOC_ERR; return NULL; }
    grp->name       = g_strdup(name ? name : _("Group"));
    grp->visible    = TRUE;
    grp->opacity    = 1.0f;
    grp->blend_mode = FLIS_BLEND_PASS_THROUGH;
    grp->collapsed  = FALSE;
    return grp;
}

void flis_group_free(flis_group_t *group) {
    if (!group) return;
    g_free(group->name);
    g_free(group->created);
    g_free(group->modified);
    g_free(group);
}

void flis_free_groups(single *uniq) {
    if (!uniq || !uniq->groups) return;
    for (GSList *g = uniq->groups; g; g = g->next)
        flis_group_free((flis_group_t *)g->data);
    g_slist_free(uniq->groups);
    uniq->groups = NULL;
}

flis_group_t *flis_group_add(const gchar *name) {
    if (!com.uniq) return NULL;
    flis_group_t *grp = flis_group_new(name);
    if (!grp) return NULL;
    grp->item_id  = com.uniq->next_item_id++;
    grp->created  = flis_now_iso8601();
    grp->modified = g_strdup(grp->created);
    com.uniq->groups = g_slist_append(com.uniq->groups, grp);
    siril_log_message(_("FLIS: created group '%s' (id=%d)\n"),
                      grp->name, grp->item_id);
    return grp;
}

int flis_group_remove(flis_group_t *group) {
    if (!com.uniq || !group) return 1;
    /* Unassign all member layers */
    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (lay && lay->group_id == group->item_id)
            lay->group_id = 0;
    }
    com.uniq->groups = g_slist_remove(com.uniq->groups, group);
    siril_log_message(_("FLIS: removed group '%s' (layers unassigned)\n"),
                      group->name ? group->name : "?");
    flis_group_free(group);
    return 0;
}

int flis_group_delete_with_layers(flis_group_t *group) {
    if (!com.uniq || !group) return 1;
    /* Collect members first (flis_layer_remove modifies the list) */
    GSList *members = flis_group_get_layers(group);
    int err = 0;
    for (GSList *m = members; m; m = m->next) {
        flis_layer_t *lay = (flis_layer_t *)m->data;
        lay->group_id = 0;  /* unassign before removal so flis_layer_remove works */
        if (flis_layer_count() > 1) {
            flis_undo_purge_layer(lay->item_id);
            if (flis_layer_remove(lay)) err = 1;
        } else {
            siril_log_warning(
                _("FLIS: cannot delete last layer — keeping it\n"));
            err = 1;
        }
    }
    g_slist_free(members);
    /* Remove the group itself */
    com.uniq->groups = g_slist_remove(com.uniq->groups, group);
    flis_group_free(group);
    return err;
}

flis_group_t *flis_group_get_by_id(gint item_id) {
    if (!com.uniq || !com.uniq->groups) return NULL;
    for (GSList *g = com.uniq->groups; g; g = g->next) {
        flis_group_t *grp = (flis_group_t *)g->data;
        if (grp && grp->item_id == item_id) return grp;
    }
    return NULL;
}

GSList *flis_group_get_layers(const flis_group_t *group) {
    if (!com.uniq || !group) return NULL;
    GSList *result = NULL;
    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (lay && lay->group_id == group->item_id)
            result = g_slist_insert_sorted(result, lay,
                                           (GCompareFunc)layer_order_cmp);
    }
    return result;
}

gint flis_group_count(void) {
    if (!com.uniq) return 0;
    return (gint)g_slist_length(com.uniq->groups);
}

int flis_layer_set_group(flis_layer_t *layer, gint group_id) {
    if (!layer) return 1;
    layer->group_id = group_id;
    flis_layer_touch_modified(layer);
    return 0;
}

int flis_group_set_name(flis_group_t *group, const gchar *name) {
    if (!group) return 1;
    g_free(group->name);
    group->name = g_strdup(name ? name : "");
    flis_group_touch_modified(group);
    return 0;
}

int flis_group_set_visible(flis_group_t *group, gboolean visible) {
    if (!group) return 1;
    group->visible = visible;
    flis_group_touch_modified(group);
    return 0;
}

int flis_group_set_opacity(flis_group_t *group, gfloat opacity) {
    if (!group) return 1;
    group->opacity = CLAMP(opacity, 0.0f, 1.0f);
    flis_group_touch_modified(group);
    return 0;
}

int flis_group_set_blend_mode(flis_group_t *group, flis_blend_mode_t mode) {
    if (!group) return 1;
    group->blend_mode = mode;
    flis_group_touch_modified(group);
    return 0;
}

void flis_group_touch_modified(flis_group_t *group) {
    if (!group) return;
    g_free(group->modified);
    group->modified = flis_now_iso8601();
}

gboolean is_current_image_flis(void) {
    /* A FLIS file is loaded when com.uniq->layers is non-NULL.
     * This list is populated exclusively by load_flis() and cleared by
     * flis_free_layers().  A regular FITS image loaded via readfits()
     * never touches com.uniq->layers, so this test is unambiguous.
     * Layer count is intentionally not checked: a single-layer FLIS is
     * still a FLIS and must go through the composite path. */
    return com.uniq != NULL && com.uniq->layers != NULL;
}

/* ===================================================================== */
/* Timestamp                                                             */
/* ===================================================================== */

void flis_layer_touch_modified(flis_layer_t *layer) {
    g_free(layer->modified);
    layer->modified = flis_now_iso8601();
}

/* ===================================================================== */
/* Stack management                                                      */
/* ===================================================================== */

flis_layer_t *flis_layer_add(fits *fit, const gchar *name) {
    if (!com.uniq) {
        siril_log_message(_("FLIS: flis_layer_add — com.uniq is NULL\n"));
        return NULL;
    }

    flis_layer_t *layer = flis_layer_new(fit, name);
    if (!layer) return NULL;

    layer->item_id     = com.uniq->next_item_id++;
    layer->layer_order = flis_max_layer_order() + FLIS_ORDER_STEP;
    layer->created     = flis_now_iso8601();
    layer->modified    = g_strdup(layer->created);

    /* g_slist_insert_sorted keeps the list sorted by layer_order */
    com.uniq->layers = g_slist_insert_sorted(com.uniq->layers, layer,
                                              (GCompareFunc)layer_order_cmp);

    /* Per-fits ICC fields removed.  A new non-base layer's pixels are
     * assumed to already be in the FLIS's colour space (the on-disk
     * profile, if any, was discarded by the load helper since fit !=
     * gfit).  Colour-space conversion at flis_layer_add time is no
     * longer supported — convert before calling, or use Convert from
     * the Image → Color Management dialog after adding. */

    /* Activate the newly added layer */
    gint idx = flis_layer_get_index(layer);
    uniq_set_active_layer(com.uniq, idx);

    siril_log_message(_("FLIS: added layer '%s' (id=%d, order=%d)\n"),
                      layer->layer_name, layer->item_id, layer->layer_order);
    return layer;
}

/* -----------------------------------------------------------------------
 * flis_addlayer_hook — shared by the panel's Add Layer toolbar handler
 * and the flis_addlayer command (§4.3 parity contract).  Reads the FITS
 * file named in args->user, derives the layer name from the basename if
 * one wasn't supplied, transfers ownership to flis_layer_add.  Returns
 * 0 on success, non-zero on failure (file read / alloc / flis_layer_add).
 * ----------------------------------------------------------------------- */

void flis_addlayer_args_free(gpointer p) {
    struct flis_addlayer_args *a = p;
    if (!a) return;
    g_free(a->filename);
    g_free(a->name);
    g_free(a);
}

int flis_addlayer_hook(struct generic_layer_args *args) {
    struct flis_addlayer_args *a = (struct flis_addlayer_args *)args->user;
    if (!a || !a->filename || !*a->filename) {
        siril_log_error(_("flis_addlayer: no filename\n"));
        return 1;
    }
    fits *f = calloc(1, sizeof(fits));
    if (!f) { PRINT_ALLOC_ERR; return 1; }
    if (readfits(a->filename, f, NULL, com.pref.force_16bit ? FALSE : TRUE)) {
        siril_log_error(_("flis_addlayer: could not read '%s'\n"), a->filename);
        free(f);
        return 1;
    }
    gchar *name = a->name;
    gchar *derived = NULL;
    if (!name || !*name) {
        gchar *base = g_path_get_basename(a->filename);
        gchar *dot  = base ? strrchr(base, '.') : NULL;
        if (dot) *dot = '\0';
        derived = base;
        name = derived;
    }
    flis_layer_t *added = flis_layer_add(f, name);
    g_free(derived);
    if (!added) {
        /* flis_layer_add already logged the failure; clear the fits we
         * tried to transfer ownership of. */
        clearfits(f);
        free(f);
        return 1;
    }
    /* args->invalidate_item_id lets end_generic_layer route the cache
     * invalidation precisely.  Even though we used FLIS_INV_ALL here
     * (stack changed), capturing the id lets the panel select the new
     * layer after refresh, which is the natural UX. */
    args->invalidate_item_id = added->item_id;
    return 0;
}

/* -----------------------------------------------------------------------
 * Shared helper: build a layermask_t from a standalone FITS file.
 *
 * The file is read via readfits (handles all the format conversions).
 * The first channel is downsampled / converted into the requested
 * @bitpix (8 or 32).  Caller takes ownership of the returned struct.
 * Returns NULL on read failure or bad bitpix.
 *
 * For bitpix=8: pixel→byte via (WORD >> 8) or (float * 255).
 * For bitpix=32: pixel→float in [0,1].
 * ----------------------------------------------------------------------- */
static layermask_t *flis_layermask_from_fits_file(const char *path, int bitpix) {
    if (bitpix != 8 && bitpix != 32) return NULL;
    fits *f = calloc(1, sizeof(fits));
    if (!f) { PRINT_ALLOC_ERR; return NULL; }
    if (readfits(path, f, NULL, TRUE)) {
        siril_log_error(_("flis: could not read mask file '%s'\n"), path);
        free(f);
        return NULL;
    }
    const size_t w = (size_t)f->rx;
    const size_t h = (size_t)f->ry;
    const size_t n = w * h;
    layermask_t *m = calloc(1, sizeof(layermask_t));
    if (!m) { PRINT_ALLOC_ERR; clearfits(f); free(f); return NULL; }
    m->w = w; m->h = h; m->bitpix = (guint8)bitpix;
    m->data = malloc(n * (bitpix / 8));
    if (!m->data) { PRINT_ALLOC_ERR; free(m); clearfits(f); free(f); return NULL; }

    const gboolean is_float = (f->type == DATA_FLOAT);
    if (bitpix == 8) {
        uint8_t *out = m->data;
        if (is_float) {
            const float *src = f->fpdata[0];
            for (size_t i = 0; i < n; i++) {
                float v = src[i] * 255.f;
                out[i] = (uint8_t)(v < 0.f ? 0 : v > 255.f ? 255 : (int)v);
            }
        } else {
            const WORD *src = f->pdata[0];
            for (size_t i = 0; i < n; i++) out[i] = (uint8_t)(src[i] >> 8);
        }
    } else {  /* bitpix == 32 */
        float *out = m->data;
        if (is_float) {
            memcpy(out, f->fpdata[0], n * sizeof(float));
        } else {
            const WORD *src = f->pdata[0];
            for (size_t i = 0; i < n; i++)
                out[i] = (float)src[i] * INV_USHRT_MAX_SINGLE;
        }
    }
    clearfits(f);
    free(f);
    return m;
}

/* -----------------------------------------------------------------------
 * flis_setmask_hook / flis_clearmask_hook — shared by the panel's mask
 * Add… / Remove… buttons and the matching flis_setmask / flis_clearmask
 * commands.  The target layer's id is in args->invalidate_item_id;
 * setmask additionally reads filename + bitpix from a flis_setmask_args
 * payload in args->user.  Both wrap flis_layer_set_lmask /
 * flis_layer_remove_lmask, which honour the layer-lock check.
 * ----------------------------------------------------------------------- */

void flis_setmask_args_free(gpointer p) {
    struct flis_setmask_args *a = p;
    if (!a) return;
    g_free(a->filename);
    g_free(a);
}

int flis_setmask_hook(struct generic_layer_args *args) {
    struct flis_setmask_args *a = (struct flis_setmask_args *)args->user;
    if (!a || !a->filename || !*a->filename) {
        siril_log_error(_("flis_setmask: no mask filename\n"));
        return 1;
    }
    flis_layer_t *lay = flis_layer_get_by_id(args->invalidate_item_id);
    if (!lay) {
        siril_log_error(_("flis_setmask: no layer with id %d\n"),
                        args->invalidate_item_id);
        return 1;
    }
    int bitpix = a->bitpix > 0 ? a->bitpix : 8;
    if (bitpix != 8 && bitpix != 32) {
        siril_log_error(_("flis_setmask: bitpix must be 8 or 32 (got %d)\n"), bitpix);
        return 1;
    }
    layermask_t *m = flis_layermask_from_fits_file(a->filename, bitpix);
    if (!m) return 1;
    /* Reject if mask dimensions don't match the layer. */
    if ((size_t)lay->fit->rx != m->w || (size_t)lay->fit->ry != m->h) {
        siril_log_error(_("flis_setmask: mask size %zux%zu does not match layer "
                          "'%s' size %ux%u\n"),
                        m->w, m->h,
                        lay->layer_name ? lay->layer_name : "?",
                        lay->fit->rx, lay->fit->ry);
        layermask_free(m);
        return 1;
    }
    /* flis_layer_set_lmask transfers ownership of the mask. */
    if (flis_layer_set_lmask(lay, m)) {
        layermask_free(m);
        return 1;
    }
    siril_log_message(_("flis_setmask: layer '%s' mask loaded from '%s'\n"),
                      lay->layer_name ? lay->layer_name : "?",
                      a->filename);
    return 0;
}

/* -----------------------------------------------------------------------
 * flis_addgroup_hook / flis_setgroup_hook (§4.3 slice 3).  Shared by the
 * panel's Create Group / Move-to-group entries and the matching
 * flis_addgroup / flis_setgroup commands.  Both produce identical state
 * changes via the existing flis_group_add / flis_layer_set_group
 * primitives.
 * ----------------------------------------------------------------------- */

void flis_addgroup_args_free(gpointer p) {
    struct flis_addgroup_args *a = p;
    if (!a) return;
    g_free(a->name);
    g_free(a);
}

int flis_addgroup_hook(struct generic_layer_args *args) {
    struct flis_addgroup_args *a = (struct flis_addgroup_args *)args->user;
    /* If no name was supplied, pick "Group N" for the smallest N not
     * currently in use.  Cheap O(N²) on group count which is tiny. */
    gchar *auto_name = NULL;
    const char *name = (a && a->name && *a->name) ? a->name : NULL;
    if (!name) {
        for (int n = 1; ; n++) {
            auto_name = g_strdup_printf(_("Group %d"), n);
            gboolean taken = FALSE;
            for (GSList *g = com.uniq ? com.uniq->groups : NULL; g; g = g->next) {
                flis_group_t *gg = (flis_group_t *)g->data;
                if (gg && gg->name && !g_strcmp0(gg->name, auto_name)) {
                    taken = TRUE; break;
                }
            }
            if (!taken) break;
            g_free(auto_name);
        }
        name = auto_name;
    }
    flis_group_t *grp = flis_group_add(name);
    g_free(auto_name);
    if (!grp) return 1;
    args->invalidate_item_id = grp->item_id;   /* for downstream selection */
    return 0;
}

void flis_setgroup_args_free(gpointer p) {
    struct flis_setgroup_args *a = p;
    if (!a) return;
    g_free(a);
}

int flis_setgroup_hook(struct generic_layer_args *args) {
    struct flis_setgroup_args *a = (struct flis_setgroup_args *)args->user;
    if (!a) return 1;
    flis_layer_t *lay = flis_layer_get_by_id(args->invalidate_item_id);
    if (!lay) {
        siril_log_error(_("flis_setgroup: no layer with id %d\n"),
                        args->invalidate_item_id);
        return 1;
    }
    if (a->group_id != 0) {
        flis_group_t *grp = flis_group_get_by_id(a->group_id);
        if (!grp) {
            siril_log_error(_("flis_setgroup: no group with id %d\n"), a->group_id);
            return 1;
        }
    }
    return flis_layer_set_group(lay, a->group_id);
}

int flis_clearmask_hook(struct generic_layer_args *args) {
    flis_layer_t *lay = flis_layer_get_by_id(args->invalidate_item_id);
    if (!lay) {
        siril_log_error(_("flis_clearmask: no layer with id %d\n"),
                        args->invalidate_item_id);
        return 1;
    }
    if (!lay->lmask) {
        siril_log_message(_("flis_clearmask: layer '%s' has no mask\n"),
                          lay->layer_name ? lay->layer_name : "?");
        return 0;
    }
    return flis_layer_remove_lmask(lay);
}

int flis_layer_remove(flis_layer_t *layer) {
    if (!com.uniq || !layer) return 1;
    if (flis_check_locked(layer, "remove layer")) return 1;

    gint count = flis_layer_count();
    if (count <= 1) {
        siril_log_warning(
            _("FLIS: cannot remove the last remaining layer\n"));
        return 1;
    }

    gint idx = flis_layer_get_index(layer);
    if (idx < 0) {
        siril_log_message(_("FLIS: flis_layer_remove — layer not found in stack\n"));
        return 1;
    }

    com.uniq->layers = g_slist_remove(com.uniq->layers, layer);

    /* Choose a sensible new active layer: prefer the one that was below;
     * if the removed layer was the base, use the new base. */
    gint new_idx = (idx > 0) ? idx - 1 : 0;
    uniq_set_active_layer(com.uniq, new_idx);

    siril_log_message(_("FLIS: removed layer '%s'\n"),
                      layer->layer_name ? layer->layer_name : "?");
    flis_layer_free(layer);
    return 0;
}

/* -------------------------------------------------------------------------
 * flis_layer_install_render
 *
 * Transfers ownership of all pixel data from @rendered (a float RGB fits*
 * returned by flis_render_layers) into @lay->fit.  The layer's existing
 * pixel buffers are freed.  ICC profile, keywords, mask, and all other
 * metadata fields already on the layer's fit are preserved.  The @rendered
 * shell is freed after the transfer.
 * ------------------------------------------------------------------------- */
static void flis_layer_install_render(flis_layer_t *lay, fits *rendered) {
    fits *f = lay->fit;

    /* Free existing pixel buffers */
    if (f->fdata) { free(f->fdata); f->fdata = NULL; }
    if (f->data)  { free(f->data);  f->data  = NULL; }
    for (int i = 0; i < 3; i++) {
        f->fpdata[i] = NULL;
        f->pdata[i]  = NULL;
    }

    /* Transfer float RGB data ownership */
    f->fdata       = rendered->fdata;
    f->fpdata[0]   = rendered->fpdata[0];
    f->fpdata[1]   = rendered->fpdata[1];
    f->fpdata[2]   = rendered->fpdata[2];
    f->type        = DATA_FLOAT;
    f->bitpix      = FLOAT_IMG;
    f->orig_bitpix = FLOAT_IMG;
    f->naxis       = rendered->naxis;
    f->naxes[0]    = rendered->naxes[0];
    f->naxes[1]    = rendered->naxes[1];
    f->naxes[2]    = rendered->naxes[2];
    f->rx          = rendered->rx;
    f->ry          = rendered->ry;

    /* ICC profile transfer.  rendered's profile was copied from com.uniq
     * (by flis_render_layers) — after flatten the FLIS is collapsing to
     * a plain FITS, so we keep the profile on com.uniq (it's still the
     * current image's profile) and just mirror onto f via the legacy
     * path until the post-flatten sequencing clears the layer stack. */
    /* Per-fits ICC fields removed.  The composite rendered by
     * flis_render_layers also no longer carries an icc_profile.  The
     * post-flatten image's ICC state continues to live on com.uniq
     * unchanged. */

    invalidate_stats_from_fit(f);

    /* Null rendered's pointers before freeing to prevent double-free */
    rendered->fdata = NULL;
    rendered->fpdata[0] = rendered->fpdata[1] = rendered->fpdata[2] = NULL;
    free(rendered);
}

/* -------------------------------------------------------------------------
 * flis_merge_down_layer
 *
 * Merges @top onto the layer directly below it in the stack.  The top
 * layer's compositing parameters (blend mode, opacity, layer mask) are
 * applied during rendering and then the layer is deleted.  The target
 * layer retains its own blend mode, opacity, and masks.
 *
 * This operation is destructive and purges the undo history for both layers.
 * Returns 0 on success, 1 on failure.
 * ------------------------------------------------------------------------- */
int flis_merge_down_layer(flis_layer_t *top) {
    if (!com.uniq || !top) return 1;
    if (flis_check_locked(top, _("merge down"))) return 1;

    /* Find the layer with the highest layer_order below top */
    flis_layer_t *bottom = NULL;
    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (!lay || lay == top) continue;
        if (lay->layer_order < top->layer_order) {
            if (!bottom || lay->layer_order > bottom->layer_order)
                bottom = lay;
        }
    }
    if (!bottom) {
        siril_log_warning(_("Merge Down: no layer below the current layer.\n"));
        return 1;
    }

    /* Render [bottom, top] as a 2-layer composite */
    GSList *tmp = g_slist_append(NULL, bottom);
    tmp = g_slist_append(tmp, top);
    fits *merged = flis_render_layers(tmp);
    g_slist_free(tmp);
    if (!merged) {
        siril_log_error(_("Merge Down: rendering failed.\n"));
        return 1;
    }

    /* Capture names for the closing log line BEFORE flis_layer_free(top)
     * dangling-pointers top->layer_name (use-after-free otherwise). */
    gchar *top_name_copy = g_strdup(top->layer_name ? top->layer_name : "?");
    const char *bottom_name = bottom->layer_name ? bottom->layer_name : "?";

    /* Destructive operation: purge undo history for both layers */
    flis_undo_purge_layer(top->item_id);
    flis_undo_purge_layer(bottom->item_id);

    /* Clear the bottom layer's masks (consumed by the merge) */
    if (bottom->fit->mask) { free_mask(bottom->fit->mask); bottom->fit->mask = NULL; }
    layermask_free(bottom->lmask);
    bottom->lmask = NULL;

    /* Install merged pixels into the bottom layer */
    flis_layer_install_render(bottom, merged);

    /* Remove the top layer from the stack (flis_layer_free handles its masks) */
    com.uniq->layers = g_slist_remove(com.uniq->layers, top);
    flis_layer_free(top);

    /* Make bottom the active layer */
    gint idx = flis_layer_get_index(bottom);
    if (idx >= 0)
        uniq_set_active_layer(com.uniq, idx);

    gui_iface.flis_invalidate_composite();

    siril_log_message(_("FLIS: merged '%s' down into '%s'\n"),
                      top_name_copy, bottom_name);
    g_free(top_name_copy);
    return 0;
}

/* -------------------------------------------------------------------------
 * flis_flatten_all
 *
 * Composites all visible layers into a single base layer.  All non-base
 * layers are deleted; the base layer's blend mode is reset to Normal and
 * all masks are cleared.
 *
 * This operation is destructive and purges the entire undo history.
 * Returns 0 on success, 1 on failure.
 * ------------------------------------------------------------------------- */
int flis_flatten_all(void) {
    if (!com.uniq || !com.uniq->layers) return 1;
    if (flis_layer_count() <= 1) return 0; /* already flat */

    /* Render all visible layers */
    fits *flat = flis_render_layers(com.uniq->layers);
    if (!flat) {
        siril_log_error(_("Flatten: rendering failed.\n"));
        return 1;
    }

    /* Base layer is the first in the sorted list (lowest layer_order) */
    flis_layer_t *base = (flis_layer_t *)com.uniq->layers->data;

    /* Remove and free all layers except the base, purging their undo entries */
    GSList *rest = g_slist_copy(com.uniq->layers->next);
    for (GSList *l = rest; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        flis_undo_purge_layer(lay->item_id);
        com.uniq->layers = g_slist_remove(com.uniq->layers, lay);
        flis_layer_free(lay);
    }
    g_slist_free(rest);

    /* Purge base's undo history too */
    flis_undo_purge_layer(base->item_id);

    /* Clear the base layer's masks */
    if (base->fit->mask) { free_mask(base->fit->mask); base->fit->mask = NULL; }
    layermask_free(base->lmask);
    base->lmask = NULL;

    /* Reset compositing parameters: base is now the sole layer */
    base->blend_mode = FLIS_BLEND_NORMAL;
    base->opacity    = 1.0f;
    base->has_tint   = FALSE;
    base->visible    = TRUE;

    /* Install the flattened pixels */
    flis_layer_install_render(base, flat);

    uniq_set_active_layer(com.uniq, 0);
    gui_iface.flis_invalidate_composite();

    siril_log_message(_("FLIS: image flattened to single layer '%s'\n"),
                      base->layer_name ? base->layer_name : "?");
    return 0;
}

flis_layer_t *flis_layer_duplicate(const flis_layer_t *src) {
    if (!com.uniq || !src) return NULL;

    /* Deep-copy the pixel data */
    fits *new_fit = calloc(1, sizeof(fits));
    if (!new_fit) { PRINT_ALLOC_ERR; return NULL; }

    if (copyfits((fits *)src->fit, new_fit,
                 CP_ALLOC | CP_COPYA | CP_FORMAT | CP_COPYMASK, 0)) {
        free(new_fit);
        return NULL;
    }

    /* copyfits(CP_FORMAT) NULLs out pointer-typed header fields (date_obs,
     * wcslib, unknown_keys) to avoid double-free.  Re-populate them with
     * deep copies so the duplicate carries complete metadata, including the
     * plate-solve solution needed for future co-registration. */
    copy_fits_metadata(src->fit, new_fit);

    gchar *new_name = g_strdup_printf("%s copy",
                      src->layer_name ? src->layer_name : "Layer");
    flis_layer_t *dup = flis_layer_new(new_fit, new_name);
    g_free(new_name);
    if (!dup) { clearfits(new_fit); free(new_fit); return NULL; }

    /* Copy compositing and metadata fields */
    dup->blend_mode  = src->blend_mode;
    dup->opacity     = src->opacity;
    dup->visible     = src->visible;
    dup->has_tint    = src->has_tint;
    dup->layer_tint  = src->layer_tint;
    dup->position_x  = src->position_x;
    dup->position_y  = src->position_y;
    dup->locked      = FALSE; /* duplicates start unlocked */
    dup->group_id    = src->group_id;

    /* Deep-copy the layer mask if present */
    if (src->lmask) {
        size_t nbytes = src->lmask->w * src->lmask->h
                        * (src->lmask->bitpix / 8);
        dup->lmask = calloc(1, sizeof(layermask_t));
        if (dup->lmask) {
            dup->lmask->w      = src->lmask->w;
            dup->lmask->h      = src->lmask->h;
            dup->lmask->bitpix = src->lmask->bitpix;
            dup->lmask->data   = malloc(nbytes);
            if (dup->lmask->data) {
                memcpy(dup->lmask->data, src->lmask->data, nbytes);
            } else {
                PRINT_ALLOC_ERR;
                layermask_free(dup->lmask);
                dup->lmask = NULL;
            }
        }
    }

    /* Assign identity and position above the source layer */
    dup->item_id     = com.uniq->next_item_id++;
    dup->layer_order = src->layer_order + FLIS_ORDER_STEP;
    dup->created     = flis_now_iso8601();
    dup->modified    = g_strdup(dup->created);

    com.uniq->layers = g_slist_insert_sorted(com.uniq->layers, dup,
                                              (GCompareFunc)layer_order_cmp);

    gint idx = flis_layer_get_index(dup);
    uniq_set_active_layer(com.uniq, idx);

    siril_log_message(_("FLIS: duplicated layer '%s' → '%s' (id=%d)\n"),
                      src->layer_name ? src->layer_name : "?",
                      dup->layer_name, dup->item_id);
    return dup;
}

/* No-op: profile state lives on com.uniq and is invariant under reorders
 * that change which layer is the base.  Kept as a stub so call sites
 * documented by the previous architecture still compile and have a place
 * to be removed in a follow-up sweep. */
static void flis_transfer_icc_to_new_base(flis_layer_t *old_base,
                                           flis_layer_t *new_base) {
    (void)old_base;
    (void)new_base;
}

int flis_layer_move_up(flis_layer_t *layer) {
    if (!com.uniq || !layer) return -1;
    if (flis_check_locked(layer, "move layer up")) return -1;

    /* The list is sorted ascending by layer_order; "up" = towards end */
    GSList *node = g_slist_find(com.uniq->layers, layer);
    if (!node) return -1;
    GSList *next = node->next;
    if (!next) return 1;  /* already at top */

    flis_layer_t *neighbour = (flis_layer_t *)next->data;
    if (!neighbour) return -1;

    /* If the neighbour belongs to a group, jump past the entire group:
     * find the topmost (highest layer_order) member of that group. */
    if (neighbour->group_id != 0) {
        flis_group_t *_grp = flis_group_get_by_id(neighbour->group_id);
        if (_grp) {
            GSList *_members = flis_group_get_layers(_grp);  /* sorted ascending */
            if (_members) {
                GSList *_last = g_slist_last(_members);
                if (_last && _last->data)
                    neighbour = (flis_layer_t *)_last->data;
                g_slist_free(_members);
            }
        }
    }

    if (flis_check_locked(neighbour, "swap with locked layer above")) return -1;

    flis_layer_t *old_base = (flis_layer_t *)com.uniq->layers->data;

    /* Swap the two layer_order values, then re-sort */
    gint tmp              = layer->layer_order;
    layer->layer_order    = neighbour->layer_order;
    neighbour->layer_order = tmp;

    flis_resort_layers(layer);

    flis_layer_t *new_base = (flis_layer_t *)com.uniq->layers->data;
    if (new_base != old_base)
        flis_transfer_icc_to_new_base(old_base, new_base);

    flis_layer_touch_modified(layer);
    return 0;
}

int flis_layer_move_down(flis_layer_t *layer) {
    if (!com.uniq || !layer) return -1;
    if (flis_check_locked(layer, "move layer down")) return -1;

    gint idx = flis_layer_get_index(layer);
    if (idx < 0) return -1;
    if (idx == 0) return 1;  /* already at bottom */

    flis_layer_t *neighbour = (flis_layer_t *)g_slist_nth_data(
                                  com.uniq->layers, idx - 1);
    if (!neighbour) return -1;

    /* If the neighbour belongs to a group, jump past the entire group:
     * find the bottommost (lowest layer_order) member of that group. */
    if (neighbour->group_id != 0) {
        flis_group_t *_grp = flis_group_get_by_id(neighbour->group_id);
        if (_grp) {
            GSList *_members = flis_group_get_layers(_grp);  /* sorted ascending */
            if (_members) {
                flis_layer_t *_first = (flis_layer_t *)_members->data;
                if (_first)
                    neighbour = _first;
                g_slist_free(_members);
            }
        }
    }

    if (flis_check_locked(neighbour, "swap with locked layer below")) return -1;

    flis_layer_t *old_base = (flis_layer_t *)com.uniq->layers->data;

    gint tmp               = layer->layer_order;
    layer->layer_order     = neighbour->layer_order;
    neighbour->layer_order = tmp;

    flis_resort_layers(layer);

    flis_layer_t *new_base = (flis_layer_t *)com.uniq->layers->data;
    if (new_base != old_base)
        flis_transfer_icc_to_new_base(old_base, new_base);

    flis_layer_touch_modified(layer);
    return 0;
}

/* ===================================================================== */
/* Compositing properties                                                */
/* ===================================================================== */

int flis_layer_set_blend_mode(flis_layer_t *layer, flis_blend_mode_t mode) {
    if (!layer) return 1;
    if (flis_check_locked(layer, "set blend mode")) return 1;
    layer->blend_mode = mode;
    flis_layer_touch_modified(layer);
    return 0;
}

int flis_layer_set_opacity(flis_layer_t *layer, gfloat opacity) {
    if (!layer) return 1;
    if (flis_check_locked(layer, "set opacity")) return 1;
    if (opacity < 0.0f) opacity = 0.0f;
    if (opacity > 1.0f) opacity = 1.0f;
    layer->opacity = opacity;
    flis_layer_touch_modified(layer);
    return 0;
}

int flis_layer_set_visible(flis_layer_t *layer, gboolean visible) {
    if (!layer) return 1;
    /* Visibility is a display toggle; allowed even on locked layers */
    layer->visible = visible;
    return 0;
}

/* ===================================================================== */
/* Layer name and lock                                                   */
/* ===================================================================== */

int flis_layer_set_name(flis_layer_t *layer, const gchar *name) {
    if (!layer) return 1;
    if (flis_check_locked(layer, "rename layer")) return 1;
    g_free(layer->layer_name);
    layer->layer_name = g_strdup(name ? name : "Layer");
    flis_layer_touch_modified(layer);
    return 0;
}

int flis_layer_set_locked(flis_layer_t *layer, gboolean locked) {
    if (!layer) return 1;
    layer->locked = locked;
    /* Deliberately does not update the modified timestamp: locking is
     * an editorial state, not a data change. */
    return 0;
}

/* ===================================================================== */
/* Layer mask (LMASK) operations                                         */
/* ===================================================================== */

int flis_layer_set_lmask(flis_layer_t *layer, layermask_t *lmask) {
    if (!layer) return 1;
    if (flis_check_locked(layer, "set layer mask")) return 1;

    if (lmask) {
        /* Dimensions must match the layer's pixel data */
        if ((size_t)layer->fit->rx != lmask->w ||
            (size_t)layer->fit->ry != lmask->h) {
            siril_log_warning(
                _("FLIS: layer mask size %zux%zu does not match layer "
                  "'%s' size %ux%u\n"),
                lmask->w, lmask->h,
                layer->layer_name ? layer->layer_name : "?",
                layer->fit->rx, layer->fit->ry);
            return 1;
        }
    }

    layermask_free(layer->lmask);
    layer->lmask = lmask;   /* NULL is valid: removes the mask */
    flis_layer_touch_modified(layer);
    return 0;
}

int flis_layer_remove_lmask(flis_layer_t *layer) {
    return flis_layer_set_lmask(layer, NULL);
}

int flis_layer_move_lmask(flis_layer_t *from, flis_layer_t *to) {
    if (!from || !to) return 1;

    if (!from->lmask) {
        siril_log_warning(
            _("FLIS: flis_layer_move_lmask — source layer '%s' has no "
              "layer mask\n"),
            from->layer_name ? from->layer_name : "?");
        return 1;
    }
    if (flis_check_locked(from, "detach layer mask")) return 1;
    if (flis_check_locked(to,   "attach layer mask")) return 1;

    /* Size check against destination */
    if ((size_t)to->fit->rx != from->lmask->w ||
        (size_t)to->fit->ry != from->lmask->h) {
        siril_log_warning(
            _("FLIS: cannot move layer mask — mask size %zux%zu does not "
              "match destination layer '%s' size %ux%u\n"),
            from->lmask->w, from->lmask->h,
            to->layer_name ? to->layer_name : "?",
            to->fit->rx, to->fit->ry);
        return 1;
    }

    /* Transfer ownership */
    layermask_free(to->lmask);
    to->lmask   = from->lmask;
    from->lmask = NULL;

    flis_layer_touch_modified(from);
    flis_layer_touch_modified(to);
    return 0;
}

/* ===================================================================== */
/* Mono layer tinting                                                    */
/* ===================================================================== */

int flis_layer_set_tint(flis_layer_t *layer, double r, double g, double b) {
    if (!layer) return 1;
    if (flis_check_locked(layer, "set tint")) return 1;

    if (!layer->fit || layer->fit->naxes[2] != 1) {
        siril_log_warning(
            _("FLIS: flis_layer_set_tint — layer '%s' is not a mono layer\n"), layer->layer_name ? layer->layer_name : "?");
        return 1;
    }

    /* Clamp each component to [0.0, 1.0] */
    layer->layer_tint.r = CLAMP(r, 0.0, 1.0);
    layer->layer_tint.g = CLAMP(g, 0.0, 1.0);
    layer->layer_tint.b = CLAMP(b, 0.0, 1.0);
    layer->has_tint     = TRUE;
    flis_layer_touch_modified(layer);
    return 0;
}

int flis_layer_clear_tint(flis_layer_t *layer) {
    if (!layer) return 1;
    if (flis_check_locked(layer, "clear tint")) return 1;
    layer->has_tint   = FALSE;
    layer->layer_tint = (flis_tint_t){ 1.0, 1.0, 1.0 };
    flis_layer_touch_modified(layer);
    return 0;
}

/* ===================================================================== */
/* Background neutralisation across tinted mono layers                   */
/* ===================================================================== */

/*
 * Solve T * a = b via the Moore-Penrose pseudoinverse.
 *
 * T is (3 × N) stored row-major: T_data[channel * N + layer].
 * b is a 3-vector (the target, typically (1,1,1)).
 * out_a receives the N-vector solution.
 *
 * Case N <= 3: SVD of T directly (3×N, M >= N).
 * Case N >  3: SVD of T^T (N×3) — gives minimum-norm solution.
 *
 * Returns 0 on success, non-zero if the matrix is rank-deficient.
 */
static int pseudoinverse_solve(const double *T_data, int N,
                                const double b[3], double *out_a) {
    const double EPS_REL = 1e-9;
    int ret = 0;

    if (N <= 3) {
        /* --- SVD of T (3×N) --- */
        gsl_matrix *A    = gsl_matrix_alloc(3, N);
        gsl_matrix *V    = gsl_matrix_alloc(N, N);
        gsl_vector *S    = gsl_vector_alloc(N);
        gsl_vector *work = gsl_vector_alloc(N);
        gsl_vector *bv   = gsl_vector_alloc(3);
        gsl_vector *tmp  = gsl_vector_alloc(N);
        gsl_vector *av   = gsl_vector_alloc(N);

        for (int r = 0; r < 3; r++)
            for (int c = 0; c < N; c++)
                gsl_matrix_set(A, r, c, T_data[r * N + c]);
        for (int r = 0; r < 3; r++)
            gsl_vector_set(bv, r, b[r]);

        gsl_linalg_SV_decomp(A, V, S, work);   /* A becomes U */

        double s_max = gsl_vector_get(S, 0);
        if (s_max < 1e-15) { ret = 1; goto cleanup_le3; }
        double eps = EPS_REL * s_max;

        /* tmp = U^T * b */
        gsl_blas_dgemv(CblasTrans, 1.0, A, bv, 0.0, tmp);
        /* apply S^+ */
        for (int k = 0; k < N; k++) {
            double sv = gsl_vector_get(S, k);
            double t  = gsl_vector_get(tmp, k);
            gsl_vector_set(tmp, k, (sv > eps) ? t / sv : 0.0);
        }
        /* a = V * tmp */
        gsl_blas_dgemv(CblasNoTrans, 1.0, V, tmp, 0.0, av);

        for (int i = 0; i < N; i++)
            out_a[i] = gsl_vector_get(av, i);

cleanup_le3:
        gsl_matrix_free(A); gsl_matrix_free(V);
        gsl_vector_free(S); gsl_vector_free(work);
        gsl_vector_free(bv); gsl_vector_free(tmp); gsl_vector_free(av);

    } else {
        /* --- SVD of T^T (N×3) for minimum-norm solution --- */
        gsl_matrix *A    = gsl_matrix_alloc(N, 3);
        gsl_matrix *V    = gsl_matrix_alloc(3, 3);
        gsl_vector *S    = gsl_vector_alloc(3);
        gsl_vector *work = gsl_vector_alloc(3);
        gsl_vector *bv   = gsl_vector_alloc(3);
        gsl_vector *tmp  = gsl_vector_alloc(3);
        gsl_vector *av   = gsl_vector_alloc(N);

        /* A = T^T: A[i][c] = T[c][i] */
        for (int r = 0; r < N; r++)
            for (int c = 0; c < 3; c++)
                gsl_matrix_set(A, r, c, T_data[c * N + r]);
        for (int r = 0; r < 3; r++)
            gsl_vector_set(bv, r, b[r]);

        gsl_linalg_SV_decomp(A, V, S, work);   /* A becomes U_t */

        double s_max = gsl_vector_get(S, 0);
        if (s_max < 1e-15) { ret = 1; goto cleanup_gt3; }
        double eps = EPS_REL * s_max;

        /* tmp = V^T * b */
        gsl_blas_dgemv(CblasTrans, 1.0, V, bv, 0.0, tmp);
        /* apply S^+ */
        for (int k = 0; k < 3; k++) {
            double sv = gsl_vector_get(S, k);
            double t  = gsl_vector_get(tmp, k);
            gsl_vector_set(tmp, k, (sv > eps) ? t / sv : 0.0);
        }
        /* a = U_t * tmp */
        gsl_blas_dgemv(CblasNoTrans, 1.0, A, tmp, 0.0, av);

        for (int i = 0; i < N; i++)
            out_a[i] = gsl_vector_get(av, i);

cleanup_gt3:
        gsl_matrix_free(A); gsl_matrix_free(V);
        gsl_vector_free(S); gsl_vector_free(work);
        gsl_vector_free(bv); gsl_vector_free(tmp); gsl_vector_free(av);
    }

    return ret;
}

/*
 * Scale every pixel in a mono layer by the given factor, clamping to [0, 1]
 * for float data or [0, USHRT_MAX] for uint16.
 */
static void scale_layer_pixels(fits *fit, double scale) {
    size_t n = fit->rx * fit->ry;
    if (fit->type == DATA_FLOAT && fit->fdata) {
        float s = (float)scale;
        float *p = fit->fpdata[0];
        for (size_t i = 0; i < n; i++)
            p[i] = CLAMP(p[i] * s, 0.0f, 1.0f);
    } else if (fit->type == DATA_USHORT && fit->data) {
        WORD *p = fit->pdata[0];
        for (size_t i = 0; i < n; i++) {
            double v = p[i] * scale;
            p[i] = (WORD)(v >= USHRT_MAX ? USHRT_MAX : (v < 0.0 ? 0 : (WORD)v));
        }
    }
}

/*
 * flis_background_neutralise — scale each mono layer so that the screen-blend
 * composite has a neutral (equal-RGB) background.
 *
 * The algorithm (from the linear approximation for small background values):
 *   M^c = Σ_i  s_i · m_i · T_i^c   (composite background in channel c)
 *
 * Setting M^R = M^G = M^B = M_target gives the 3×N linear system T·a = 1
 * where a_i = s_i·m_i and T_i^c is the normalised tint of layer i.
 * Solved via SVD pseudoinverse; normalised so that the overall composite
 * brightness is preserved.
 *
 * Only mono (single-channel) layers participate; RGB layers are skipped.
 *
 * Returns 0 on success, non-zero on failure.
 */
int flis_background_neutralise(void) {
    if (!is_current_image_flis() || !com.uniq || !com.uniq->layers) return 1;

    /* --- Collect eligible layers (mono only) --- */
    int total = g_slist_length(com.uniq->layers);
    flis_layer_t **layers_arr = calloc(total, sizeof(flis_layer_t *));
    double *medians = calloc(total, sizeof(double));
    double *T_data  = calloc(3 * total, sizeof(double)); /* row-major [ch][layer] */
    if (!layers_arr || !medians || !T_data) {
        PRINT_ALLOC_ERR;
        free(layers_arr); free(medians); free(T_data);
        return 1;
    }

    int N = 0;
    for (GSList *l = com.uniq->layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (!lay || !lay->fit) continue;
        if (lay->fit->naxes[2] != 1) continue;   /* skip RGB layers */
        if (!lay->fit->fdata && !lay->fit->data)  continue;

        /* Tint (already normalised [0,1]; default {1,1,1} for untinted) */
        T_data[0 * total + N] = lay->layer_tint.r;
        T_data[1 * total + N] = lay->layer_tint.g;
        T_data[2 * total + N] = lay->layer_tint.b;

        /* Background median (normalised to [0,1]) */
        imstats *st = statistics(NULL, -1, lay->fit, 0, NULL, STATS_BASIC, SINGLE_THREADED);
        if (!st) { free(layers_arr); free(medians); free(T_data); return 1; }
        double med = st->median;
        free_stats(st);
        if (lay->fit->type == DATA_USHORT)
            med /= USHRT_MAX_DOUBLE;
        medians[N] = med;

        layers_arr[N] = lay;
        N++;
    }

    if (N == 0) {
        siril_log_warning(_("FLIS: background neutralise — no eligible mono layers found\n"));
        free(layers_arr); free(medians); free(T_data);
        return 1;
    }

    /* Compact T to N columns */
    double *T = calloc(3 * N, sizeof(double));
    if (!T) { PRINT_ALLOC_ERR; free(layers_arr); free(medians); free(T_data); return 1; }
    for (int c = 0; c < 3; c++)
        for (int i = 0; i < N; i++)
            T[c * N + i] = T_data[c * total + i];
    free(T_data);

    /* Compute old_bg: average composite background brightness before any scaling.
     * This is the per-channel target we want to preserve after neutralisation. */
    double old_bg = 0.0;
    for (int i = 0; i < N; i++)
        old_bg += medians[i] * (T[0*N+i] + T[1*N+i] + T[2*N+i]);
    old_bg /= 3.0;

    if (old_bg <= 0.0) {
        siril_log_warning(_("FLIS: background neutralise — background level is zero, nothing to do\n"));
        free(layers_arr); free(medians); free(T);
        return 1;
    }

    /* --- Solve T · a = (1,1,1)^T for ALL layers simultaneously ---
     *
     * a_i = s_i · m_i is the target composite contribution of layer i.
     * The actual scale factor is then s_i = old_bg · a_i / m_i.
     *
     * Key insight: even when a layer requires a_i < 0 (infeasible with positive
     * scaling), the OTHER layers' coefficients from this global solve still correctly
     * balance the channels where the infeasible layer has zero tint.  Re-solving
     * after removing the infeasible layer destroys that balance.  Therefore we
     * always apply the coefficients from this single global solve, and merely
     * leave any infeasible layer unscaled.
     */
    double *a = calloc(N, sizeof(double));
    if (!a) { PRINT_ALLOC_ERR; free(layers_arr); free(medians); free(T); return 1; }

    const double b_unit[3] = { 1.0, 1.0, 1.0 };
    if (pseudoinverse_solve(T, N, b_unit, a)) {
        siril_log_warning(_("FLIS: background neutralise — tint matrix is rank-deficient, cannot solve\n"));
        free(layers_arr); free(medians); free(T); free(a);
        return 1;
    }

    /* Flag layers whose tint geometry requires a negative scale factor.
     * These are left unscaled (s_i = 1).  The channels where such a layer has
     * non-zero tint will not be fully neutralised; channels where it has zero
     * tint are unaffected and will be correctly balanced by the other layers.
     *
     * Feasibility condition: the vector (1,1,1) must lie in the positive span of
     * the tint vectors.  An infeasible layer is one whose tint "over-contributes"
     * to some channels relative to what the other layers can compensate for.
     * The fix is always to make the infeasible layer's tint more balanced across
     * all three channels (ideally equal R=G=B, or at least no single channel
     * significantly dominant). */
    for (int i = 0; i < N; i++) {
        if (a[i] <= 0.0 || medians[i] <= 0.0) {
            const char *name = layers_arr[i]->layer_name ? layers_arr[i]->layer_name : "?";

            /* Find the dominant channel in this layer's tint */
            double tr = T[0*N+i], tg = T[1*N+i], tb = T[2*N+i];
            const char *dom = (tr >= tg && tr >= tb) ? "R"
                            : (tg >= tr && tg >= tb) ? "G" : "B";
            const char *low = (tr <= tg && tr <= tb) ? "R"
                            : (tg <= tr && tg <= tb) ? "G" : "B";

            siril_log_warning(
                _("FLIS: background neutralise — layer '%s' is incompatible with "
                  "neutral balance (coefficient %.4f < 0).\n"
                  "  Its tint (%s-dominant) cannot be compensated by the other layers.\n"
                  "  To fix: change this layer's tint so that the %s component is "
                  "increased (aim for equal R=G=B, or at least no dominant primary).\n"
                  "  Example: for an SHO palette, assign each layer to a distinct "
                  "primary (e.g. SII=#FF0000, Ha=#00FF00, OIII=#0000FF).\n"), name, a[i], dom, low);
            a[i] = -1.0;  /* sentinel: leave unscaled */
        }
    }

    /* --- Apply scale factors and report predicted composite --- */
    siril_log_info(_("FLIS: background neutralise scale factors:\n"));
    double new_bg[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < N; i++) {
        const char *name = layers_arr[i]->layer_name ? layers_arr[i]->layer_name : "?";
        double s;
        if (a[i] < 0.0) {
            s = 1.0;
            siril_log_warning("  %-24s  ×1.0000  (left unscaled — tint infeasible)\n", name);
        } else {
            s = (old_bg * a[i]) / medians[i];
            siril_log_info("  %-24s  ×%.4f  (median %.5f → %.5f)\n",
                name, s, medians[i], medians[i] * s);
            scale_layer_pixels(layers_arr[i]->fit, s);
            invalidate_stats_from_fit(layers_arr[i]->fit);
        }
        for (int c = 0; c < 3; c++)
            new_bg[c] += s * medians[i] * T[c * N + i];
    }
    siril_log_info(
        _("FLIS: predicted composite background  R:%.5f  G:%.5f  B:%.5f  (target %.5f each)\n"), new_bg[0], new_bg[1], new_bg[2], old_bg);

    free(layers_arr); free(medians); free(T); free(a);
    gui_iface.flis_invalidate_composite();
    return 0;
}

/*
 * flis_background_neutralise_layers — same as flis_background_neutralise but
 * operates only on the specified layer subset (a GSList of flis_layer_t*).
 * If @layer_subset is NULL, falls back to operating on all layers.
 */
int flis_background_neutralise_layers(GSList *layer_subset) {
    if (!is_current_image_flis() || !com.uniq || !com.uniq->layers) return 1;
    GSList *target = layer_subset ? layer_subset : com.uniq->layers;

    int total = g_slist_length(target);
    flis_layer_t **layers_arr = calloc(total, sizeof(flis_layer_t *));
    double *medians = calloc(total, sizeof(double));
    double *T_data  = calloc(3 * total, sizeof(double));
    if (!layers_arr || !medians || !T_data) {
        PRINT_ALLOC_ERR;
        free(layers_arr); free(medians); free(T_data);
        return 1;
    }

    int N = 0;
    for (GSList *l = target; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (!lay || !lay->fit) continue;
        if (lay->fit->naxes[2] != 1) continue;
        if (!lay->fit->fdata && !lay->fit->data) continue;

        T_data[0 * total + N] = lay->has_tint ? lay->layer_tint.r : 1.0;
        T_data[1 * total + N] = lay->has_tint ? lay->layer_tint.g : 1.0;
        T_data[2 * total + N] = lay->has_tint ? lay->layer_tint.b : 1.0;

        imstats *st = statistics(NULL, -1, lay->fit, 0, NULL, STATS_BASIC, SINGLE_THREADED);
        if (!st) { free(layers_arr); free(medians); free(T_data); return 1; }
        double med = st->median;
        free_stats(st);
        if (lay->fit->type == DATA_USHORT)
            med /= USHRT_MAX_DOUBLE;
        medians[N] = med;
        layers_arr[N] = lay;
        N++;
    }

    if (N == 0) {
        siril_log_warning(_("FLIS: background neutralise — no eligible mono layers found\n"));
        free(layers_arr); free(medians); free(T_data);
        return 1;
    }

    double *T = calloc(3 * N, sizeof(double));
    if (!T) { PRINT_ALLOC_ERR; free(layers_arr); free(medians); free(T_data); return 1; }
    for (int c = 0; c < 3; c++)
        for (int i = 0; i < N; i++)
            T[c * N + i] = T_data[c * total + i];
    free(T_data);

    double old_bg = 0.0;
    for (int i = 0; i < N; i++)
        old_bg += medians[i] * (T[0*N+i] + T[1*N+i] + T[2*N+i]);
    old_bg /= 3.0;

    if (old_bg <= 0.0) {
        siril_log_warning(_("FLIS: background neutralise — background level is zero\n"));
        free(layers_arr); free(medians); free(T);
        return 1;
    }

    double *a = calloc(N, sizeof(double));
    if (!a) { PRINT_ALLOC_ERR; free(layers_arr); free(medians); free(T); return 1; }

    const double b_unit[3] = { 1.0, 1.0, 1.0 };
    if (pseudoinverse_solve(T, N, b_unit, a)) {
        siril_log_warning(_("FLIS: background neutralise — tint matrix rank-deficient\n"));
        free(layers_arr); free(medians); free(T); free(a);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        if (a[i] <= 0.0 || medians[i] <= 0.0) {
            const char *name = layers_arr[i]->layer_name ? layers_arr[i]->layer_name : "?";
            double tr = T[0*N+i], tg = T[1*N+i], tb = T[2*N+i];
            const char *dom = (tr >= tg && tr >= tb) ? "R" : (tg >= tr && tg >= tb) ? "G" : "B";
            const char *low = (tr <= tg && tr <= tb) ? "R" : (tg <= tr && tg <= tb) ? "G" : "B";
            siril_log_warning(
                _("FLIS: background neutralise — layer '%s' incompatible "
                  "(coeff %.4f < 0; %s-dominant tint, increase %s component)\n"), name, a[i], dom, low);
            a[i] = -1.0;
        }
    }

    siril_log_info(_("FLIS: background neutralise scale factors:\n"));
    double new_bg[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < N; i++) {
        const char *name = layers_arr[i]->layer_name ? layers_arr[i]->layer_name : "?";
        double s;
        if (a[i] < 0.0) {
            s = 1.0;
            siril_log_warning("  %-24s  x1.0000  (left unscaled — tint infeasible)\n", name);
        } else {
            s = (old_bg * a[i]) / medians[i];
            siril_log_info("  %-24s  x%.4f  (median %.5f -> %.5f)\n",
                name, s, medians[i], medians[i] * s);
            scale_layer_pixels(layers_arr[i]->fit, s);
            invalidate_stats_from_fit(layers_arr[i]->fit);
        }
        for (int c = 0; c < 3; c++)
            new_bg[c] += s * medians[i] * T[c * N + i];
    }
    siril_log_info(
        _("FLIS: predicted composite background  R:%.5f  G:%.5f  B:%.5f  (target %.5f each)\n"), new_bg[0], new_bg[1], new_bg[2], old_bg);

    free(layers_arr); free(medians); free(T); free(a);
    gui_iface.flis_invalidate_composite();
    return 0;
}

