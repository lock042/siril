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

#include "core/siril.h"
#include "core/proto.h"
#include "core/masks.h"
#include "core/siril_log.h"
#include "core/icc_profile.h"
#include "image_format_fits.h"
#include "image_format_flis.h"

/* -----------------------------------------------------------------------
 * FLIS_META binary table column definitions.
 *
 * The METADATA column uses a fixed 256A string rather than the 1PA
 * variable-length form described in the spec.  This avoids CFITSIO heap
 * complexity; 256 characters is sufficient for all current metadata keys.
 * A future revision may switch to 1PA when more metadata keys are defined.
 * ----------------------------------------------------------------------- */
#define FLIS_META_NCOLS      13

static const char *FLIS_COL_NAMES[FLIS_META_NCOLS] = {
    "ITEM_ID", "ITEM_TYPE", "HDU_INDEX", "PARENT_ID", "LAYER_ORDER",
    "LAYER_NAME", "COLOR_MDL", "BLEND_MODE", "OPACITY", "VISIBLE",
    "POSITION_X", "POSITION_Y", "METADATA"
};
static const char *FLIS_COL_FMTS[FLIS_META_NCOLS] = {
    "1J", "8A", "1J", "1J", "1J",
    "32A", "4A", "16A", "1E", "1L",
    "1J", "1J", "256A"
};
static const char *FLIS_COL_UNITS[FLIS_META_NCOLS] = {
    "", "", "", "", "",
    "", "", "", "", "",
    "pix", "pix", ""
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

/* ITEM_TYPE string constants */
#define FLIS_TYPE_LAYER  "LAYER"
#define FLIS_TYPE_MASK   "MASK"
#define FLIS_TYPE_LMASK  "LMASK"

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
    char    metadata[257];   /* 256A + NUL */
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
        siril_log_color_message(_("FLIS: thumbnail generation failed, writing empty HDU\n"), "salmon");
        /* Fall through to write a degenerate 1x1 placeholder */
        tw = 1; th = 1;
        thumb = calloc(1, 1);
        if (!thumb) return 1;
    }

    long nchans = (base_fit && base_fit->naxes[2] >= 3) ? 3 : 1;
    int naxis = (nchans == 3) ? 3 : 2;
    long naxes[3] = { tw, th, nchans };

    /* Resize the primary HDU to hold the thumbnail */
    if (fits_resize_img(fptr, BYTE_IMG, naxis, naxes, &status)) {
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
            /* Truncate metadata to 256 chars */
            if (meta && strlen(meta) > 256) meta[256] = '\0';
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
            g_free(mname);
            row++;
        }

        if (status) {
            report_fits_error(status);
            return 1;
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
                          int bitpix_fits, int cfitsio_type) {
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
        strptr = row->metadata;
        fits_read_col(fptr, TSTRING, COL_METADATA,    r, 1, 1, &anull,     &strptr,            NULL, &status);

        if (status) {
            siril_log_color_message(_("FLIS: error reading metadata row %ld: %d\n"), "salmon", r, status);
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

    layer->fit        = fit;
    layer->layer_name = g_strdup(name ? name : "Layer");
    layer->blend_mode = FLIS_BLEND_NORMAL;
    layer->opacity    = 1.0f;
    layer->visible    = TRUE;
    layer->has_tint   = FALSE;
    layer->layer_tint = (flis_tint_t){ 1.0, 1.0, 1.0 };
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

    /* Ensure the filename has the .flis extension */
    gchar *outpath;
    if (g_str_has_suffix(filename, ".flis")) {
        outpath = g_strdup(filename);
    } else {
        outpath = g_strdup_printf("%s.flis", filename);
    }

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

    /* Determine whether an ICC profile should be written.
     * Use the base layer's ICC profile for the whole file. */
    gboolean icc_present = (base->fit &&
                            base->fit->icc_profile != NULL &&
                            base->fit->color_managed &&
                            com.pref.fits_save_icc);

    /* ----------------------------------------------------------------
     * HDU 1 (index 0): primary thumbnail HDU
     * ---------------------------------------------------------------- */
    if (write_thumbnail_hdu(fptr, base->fit, canvas_w, canvas_h, icc_present)) {
        siril_log_color_message(_("FLIS: failed to write thumbnail HDU\n"), "red");
        fits_close_file(fptr, &status);
        g_slist_free(sorted);
        g_free(outpath);
        return 1;
    }

    /* ----------------------------------------------------------------
     * HDU 2 (optional): ICC profile
     * ---------------------------------------------------------------- */
    if (icc_present) {
        if (write_icc_profile_to_fptr(fptr, base->fit->icc_profile)) {
            siril_log_color_message(_("FLIS: warning — ICC profile write failed, continuing\n"), "salmon");
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
        siril_log_color_message(_("FLIS: failed to write FLIS_META table\n"), "red");
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
            siril_log_color_message(_("FLIS: failed writing layer '%s'\n"), "red",
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
                               BYTE_IMG, TBYTE)) {
                siril_log_color_message(_("FLIS: failed writing layer mask for '%s'\n"),
                                        "salmon", lay->layer_name ? lay->layer_name : "?");
            }
            g_free(lname);
        }

        /* --- Processing mask HDU (MASK) --- */
        if (lay->fit && lay->fit->mask) {
            /* Processing masks in Siril are 32-bit float */
            mask_t *pmask = lay->fit->mask;
            gchar *mname = g_strdup_printf("%s Processing Mask",
                           lay->layer_name ? lay->layer_name : "Layer");
            int pmask_id = lay->item_id + 20000;

            /* Wrap in a layermask_t for the generic writer */
            layermask_t tmp_lm = {
                .w      = lay->fit->rx,
                .h      = lay->fit->ry,
                .bitpix = 32,
                .data   = pmask->data
            };
            if (write_mask_hdu(fptr, &tmp_lm, pmask_id,
                               FLIS_TYPE_MASK, mname,
                               FLOAT_IMG, TFLOAT)) {
                siril_log_color_message(_("FLIS: failed writing processing mask for '%s'\n"),
                                        "salmon", lay->layer_name ? lay->layer_name : "?");
            }
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
        siril_log_color_message(_("FLIS: file %s is not a valid FLIS file (no FLIS=T keyword)\n"),
                                "red", filename);
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
        siril_log_color_message(_("FLIS: cannot locate FLIS_META table\n"), "red");
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
            siril_log_color_message(_("FLIS: FLIS_META table not found\n"), "red");
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

    /* Pass 1: LAYER rows */
    for (long r = 0; r < nrows; r++) {
        flis_meta_row_t *row = &meta_rows[r];
        if (g_ascii_strcasecmp(row->item_type, FLIS_TYPE_LAYER))
            continue; /* not a LAYER row */

        fits *f = load_layer_from_hdu(fptr, row->hdu_index);
        if (!f) {
            siril_log_color_message(
                _("FLIS: failed to load layer '%s' (HDU %d), skipping\n"),
                "salmon", row->layer_name, row->hdu_index);
            continue;
        }

        /* Assign the file-level ICC profile to each layer */
        if (file_icc) {
            f->icc_profile = copyICCProfile(file_icc);
            color_manage(f, TRUE);
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

        layers = g_slist_insert_sorted(layers, layer,
                                       (GCompareFunc)layer_order_cmp);
        g_hash_table_insert(id_map,
                            GINT_TO_POINTER(row->item_id), layer);
    }

    /* Pass 2: LMASK and MASK rows */
    for (long r = 0; r < nrows; r++) {
        flis_meta_row_t *row = &meta_rows[r];
        gboolean is_lmask = !g_ascii_strcasecmp(row->item_type, FLIS_TYPE_LMASK);
        gboolean is_mask  = !g_ascii_strcasecmp(row->item_type, FLIS_TYPE_MASK);
        if (!is_lmask && !is_mask) continue;

        flis_layer_t *parent = g_hash_table_lookup(id_map,
                                    GINT_TO_POINTER(row->parent_id));
        if (!parent) {
            siril_log_color_message(
                _("FLIS: mask '%s' references unknown parent ID %d, skipping\n"),
                "salmon", row->layer_name, row->parent_id);
            continue;
        }

        if (is_lmask) {
            layermask_t *lm = load_mask_from_hdu(fptr, row->hdu_index, FALSE);
            if (lm) {
                layermask_free(parent->lmask); /* replace any existing */
                parent->lmask = lm;
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
                    siril_mask->bitpix = 32;
                    siril_mask->data   = pm->data;
                    pm->data           = NULL; /* transfer ownership */
                    parent->fit->mask  = siril_mask;
                }
                layermask_free(pm);
            } else {
                layermask_free(pm);
            }
        }
    }

    g_hash_table_destroy(id_map);
    free(meta_rows);
    if (file_icc) cmsCloseProfile(file_icc);
    fits_close_file(fptr, &status);

    if (!layers) {
        siril_log_color_message(_("FLIS: no layers loaded from %s\n"), "red", filename);
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
    com.uniq->next_item_id = 1; /* will be set properly below */

    /* Find the highest item_id to seed next_item_id */
    gint max_id = 0;
    for (GSList *l = layers; l; l = l->next) {
        flis_layer_t *lay = (flis_layer_t *)l->data;
        if (lay->item_id > max_id) max_id = lay->item_id;
    }
    com.uniq->next_item_id = max_id + 1;

    /* Set active layer to the base (index 0 = lowest layer_order after sort) */
    uniq_set_active_layer(com.uniq, 0);

    siril_log_message(_("FLIS: loaded %d layer(s) from %s (%dx%d canvas)\n"),
                      g_slist_length(layers), filename, canvas_w, canvas_h);
    return 0;
}
