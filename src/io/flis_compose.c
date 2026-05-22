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
 * flis_compose.c — FLIS layer compositing kernel.
 *
 * Extracted from the original flis-branch src/gui/image_display.c so that
 * the pure-algorithm compositing kernel lives in src/io/ alongside the
 * FLIS container format itself.  The kernel takes a GSList of
 * flis_layer_t* and returns a newly allocated fits* containing the float
 * RGB composite.  It is callable from non-GUI code paths:
 *   - merge-down  (flis_merge_down_layer in image_format_flis.c)
 *   - flatten     (flis_flatten_all)
 *   - thumbnail   (save-time HDU 0 generation)
 *   - tests       (per-blend-mode pixel-correctness suites in stage 1.7)
 *
 * The display layer (src/gui-gtk4/image_display.c, stage 3) keeps the
 * composite cache and the texture upload pipeline; it calls this kernel
 * when the per-layer GdkTexture path is unavailable (fallback) and uses
 * it as the correctness oracle for the GPU compose path.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "core/icc_profile.h"
#include "io/image_format_flis.h"
#include "io/image_format_fits.h"
#include "flis_compose.h"

static fits *flis_render_layers_internal(GSList *layers, gboolean sub_composite);

fits *flis_render_layers(GSList *layers) {
    return flis_render_layers_internal(layers, FALSE);
}

/* =========================================================================
 * FLIS compositing inner-loop macros
 *
 * Option B: the blend mode switch is hoisted above all pixel loops.
 *   Each case of FLIS_DISPATCH expands to a dedicated inner loop whose
 *   body is a pure arithmetic expression with no branches.  The compiler
 *   can auto-vectorise each instantiation independently.
 *
 * Option C: masked and unmasked loops are separate macro expansions.
 *   Neither contains a per-pixel branch for the mask check.
 *
 * Variable conventions inside macro bodies (underscore prefix avoids
 * clashing with outer-scope names):
 *   _y, _x, _i, _row  -- loop indices and pixel offset
 *   _alpha             -- per-pixel effective alpha (masked variant only)
 *   _c, _k             -- used in the USHORT conversion pass
 *
 * Blend expression variables (declared inside each inner loop):
 *   sr, sg, sb  -- source pixel RGB [0,1] after tint applied
 *   dr, dg, db  -- destination pixel RGB [0,1] from previous composite
 *
 * Tint scalars TR/TG/TB are loop-invariant float constants.  For RGB or
 * untinted mono layers they are 1.f, which the compiler folds so that
 * pre_r[_i]*1.f compiles to a plain load.
 * ========================================================================= */

/* OMP pragma helpers -- emit nothing when OpenMP is not compiled in */
#ifdef _OPENMP
#  define FLIS_OMP_PAR_FOR \
     _Pragma("omp parallel for num_threads(com.max_thread) schedule(static)")
#  define FLIS_OMP_SIMD          _Pragma("omp simd")
#  define FLIS_OMP_PAR_FOR_SIMD \
     _Pragma("omp parallel for simd num_threads(com.max_thread) schedule(static)")
#else
#  define FLIS_OMP_PAR_FOR
#  define FLIS_OMP_SIMD
#  define FLIS_OMP_PAR_FOR_SIMD
#endif

/* Clamp to [0,1] via fminf/fmaxf -- these map to VMINPS/VMAXPS on AVX */
#define FLIS_C01(v)  fminf(fmaxf((v), 0.f), 1.f)

/* Porter-Duff 'over': blended*alpha + dst*(1-alpha) */
#define FLIS_OVER(b, d, a)  FLIS_C01((b) * (a) + (d) * (1.f - (a)))

/* Soft-light sub-expressions (W3C formula, two nested ternaries).
 * Both sqrtf branches are computed by the SIMD unit and blended with
 * a select; this is vectorisable on all modern ISAs. */
#define FLIS_SL_D(d) \
    ((d) <= 0.25f ? (((16.f*(d) - 12.f)*(d) + 4.f)*(d)) : sqrtf(d))
#define FLIS_SL(s, d) \
    ((s) <= 0.5f \
        ? (d) - (1.f - 2.f*(s))*(d)*(1.f-(d)) \
        : (d) + (2.f*(s) - 1.f)*(FLIS_SL_D(d) - (d)))

/* Color-dodge and color-burn with division-by-zero guard */
#define FLIS_DODGE(s, d) \
    ((s) >= 1.f ? 1.f : fminf((d) / (1.f - (s)), 1.f))
#define FLIS_BURN(s, d) \
    ((s) <= 0.f ? 0.f : 1.f - fminf((1.f - (d)) / (s), 1.f))

/* -------------------------------------------------------------------------
 * Unified blend loop macros — work for both full-canvas and sparse layers.
 *
 * Variables from the surrounding scope of flis_composite_build:
 *   y0, y1  — intersection row range  (0..H for full-canvas)
 *   x0, x1  — intersection col range  (0..W for full-canvas)
 *   oy, ox  — layer origin on canvas  (0,0 for full-canvas)
 *   lW      — layer pixel width       (= W for full-canvas)
 *
 * Zero-pixel convention: source pixels where all channels are 0 are "no
 * data" (rotation/undistortion borders).  _valid zeroes the alpha
 * branchlessly for such pixels, leaving the destination unchanged.
 * FLIS_BASE_LOOP is exempt — zero base pixels initialise the canvas.
 * ------------------------------------------------------------------------- */
#define FLIS_BASE_LOOP(TR, TG, TB) \
do { \
    FLIS_OMP_PAR_FOR \
    for (gint _y = y0; _y < y1; _y++) { \
        const size_t _lrow = (size_t)(_y - oy) * lW + (size_t)(x0 - ox); \
        const size_t _crow = (size_t)_y * W    + (size_t)x0; \
        FLIS_OMP_SIMD \
        for (gint _dx = 0; _dx < (x1 - x0); _dx++) { \
            const size_t _li = _lrow + (size_t)_dx; \
            const size_t _i  = _crow + (size_t)_dx; \
            out_r[_i] = FLIS_C01(pre_r[_li] * (TR)); \
            out_g[_i] = FLIS_C01(pre_g[_li] * (TG)); \
            out_b[_i] = FLIS_C01(pre_b[_li] * (TB)); \
        } \
    } \
} while (0)

#define FLIS_BLEND_LOOP(TR, TG, TB, BR, BG, BB) \
do { \
    FLIS_OMP_PAR_FOR \
    for (gint _y = y0; _y < y1; _y++) { \
        const size_t _lrow = (size_t)(_y - oy) * lW + (size_t)(x0 - ox); \
        const size_t _crow = (size_t)_y * W    + (size_t)x0; \
        FLIS_OMP_SIMD \
        for (gint _dx = 0; _dx < (x1 - x0); _dx++) { \
            const size_t _li = _lrow + (size_t)_dx; \
            const size_t _i  = _crow + (size_t)_dx; \
            const float sr = pre_r[_li] * (TR); \
            const float sg = pre_g[_li] * (TG); \
            const float sb = pre_b[_li] * (TB); \
            const float _valid = (float)((sr + sg + sb) > 0.f); \
            const float dr = out_r[_i]; \
            const float dg = out_g[_i]; \
            const float db = out_b[_i]; \
            out_r[_i] = FLIS_OVER((BR), dr, global_opacity * _valid); \
            out_g[_i] = FLIS_OVER((BG), dg, global_opacity * _valid); \
            out_b[_i] = FLIS_OVER((BB), db, global_opacity * _valid); \
        } \
    } \
} while (0)

#define FLIS_BLEND_LOOP_MASKED(TR, TG, TB, BR, BG, BB) \
do { \
    FLIS_OMP_PAR_FOR \
    for (gint _y = y0; _y < y1; _y++) { \
        const size_t _lrow = (size_t)(_y - oy) * lW + (size_t)(x0 - ox); \
        const size_t _crow = (size_t)_y * W    + (size_t)x0; \
        FLIS_OMP_SIMD \
        for (gint _dx = 0; _dx < (x1 - x0); _dx++) { \
            const size_t _li   = _lrow + (size_t)_dx; \
            const size_t _i    = _crow + (size_t)_dx; \
            const float sr = pre_r[_li] * (TR); \
            const float sg = pre_g[_li] * (TG); \
            const float sb = pre_b[_li] * (TB); \
            const float _valid = (float)((sr + sg + sb) > 0.f); \
            const float _alpha = global_opacity \
                                 * (mask_data[_li] * INV_UCHAR_MAX_SINGLE) \
                                 * _valid; \
            const float dr = out_r[_i]; \
            const float dg = out_g[_i]; \
            const float db = out_b[_i]; \
            out_r[_i] = FLIS_OVER((BR), dr, _alpha); \
            out_g[_i] = FLIS_OVER((BG), dg, _alpha); \
            out_b[_i] = FLIS_OVER((BB), db, _alpha); \
        } \
    } \
} while (0)

/* HSL loop: opacity is passed into flis_blend_hsl_pixel and consumed there
 * (component-space interpolation).  No outer FLIS_OVER is applied. */
#define FLIS_HSL_LOOP(TR, TG, TB) \
do { \
    FLIS_OMP_PAR_FOR \
    for (gint _y = y0; _y < y1; _y++) { \
        const size_t _lrow = (size_t)(_y - oy) * lW + (size_t)(x0 - ox); \
        const size_t _crow = (size_t)_y * W    + (size_t)x0; \
        for (gint _dx = 0; _dx < (x1 - x0); _dx++) { \
            const size_t _li = _lrow + (size_t)_dx; \
            const size_t _i  = _crow + (size_t)_dx; \
            const float sr = pre_r[_li] * (TR); \
            const float sg = pre_g[_li] * (TG); \
            const float sb = pre_b[_li] * (TB); \
            const float _valid = (float)((sr + sg + sb) > 0.f); \
            const float dr = out_r[_i]; \
            const float dg = out_g[_i]; \
            const float db = out_b[_i]; \
            flis_blend_hsl_pixel(mode, global_opacity * _valid, \
                                 sr, sg, sb, dr, dg, db, \
                                 &out_r[_i], &out_g[_i], &out_b[_i]); \
        } \
    } \
} while (0)

#define FLIS_HSL_LOOP_MASKED(TR, TG, TB) \
do { \
    FLIS_OMP_PAR_FOR \
    for (gint _y = y0; _y < y1; _y++) { \
        const size_t _lrow = (size_t)(_y - oy) * lW + (size_t)(x0 - ox); \
        const size_t _crow = (size_t)_y * W    + (size_t)x0; \
        for (gint _dx = 0; _dx < (x1 - x0); _dx++) { \
            const size_t _li   = _lrow + (size_t)_dx; \
            const size_t _i    = _crow + (size_t)_dx; \
            const float sr = pre_r[_li] * (TR); \
            const float sg = pre_g[_li] * (TG); \
            const float sb = pre_b[_li] * (TB); \
            const float _valid = (float)((sr + sg + sb) > 0.f); \
            const float _alpha = global_opacity \
                                 * (mask_data[_li] * INV_UCHAR_MAX_SINGLE) \
                                 * _valid; \
            const float dr = out_r[_i]; \
            const float dg = out_g[_i]; \
            const float db = out_b[_i]; \
            flis_blend_hsl_pixel(mode, _alpha, \
                                 sr, sg, sb, dr, dg, db, \
                                 &out_r[_i], &out_g[_i], &out_b[_i]); \
        } \
    } \
} while (0)

/* CHROMA loop: calls flis_blend_chroma_pixel which handles opacity
 * via FLIS_OVER internally. */
#define FLIS_CHROMA_LOOP(TR, TG, TB) \
do { \
    FLIS_OMP_PAR_FOR \
    for (gint _y = y0; _y < y1; _y++) { \
        const size_t _lrow = (size_t)(_y - oy) * lW + (size_t)(x0 - ox); \
        const size_t _crow = (size_t)_y * W    + (size_t)x0; \
        for (gint _dx = 0; _dx < (x1 - x0); _dx++) { \
            const size_t _li = _lrow + (size_t)_dx; \
            const size_t _i  = _crow + (size_t)_dx; \
            const float sr = pre_r[_li] * (TR); \
            const float sg = pre_g[_li] * (TG); \
            const float sb = pre_b[_li] * (TB); \
            const float _valid = (float)((sr + sg + sb) > 0.f); \
            const float dr = out_r[_i]; \
            const float dg = out_g[_i]; \
            const float db = out_b[_i]; \
            flis_blend_chroma_pixel(global_opacity * _valid, \
                                    sr, sg, sb, dr, dg, db, \
                                    &out_r[_i], &out_g[_i], &out_b[_i]); \
        } \
    } \
} while (0)

#define FLIS_CHROMA_LOOP_MASKED(TR, TG, TB) \
do { \
    FLIS_OMP_PAR_FOR \
    for (gint _y = y0; _y < y1; _y++) { \
        const size_t _lrow = (size_t)(_y - oy) * lW + (size_t)(x0 - ox); \
        const size_t _crow = (size_t)_y * W    + (size_t)x0; \
        for (gint _dx = 0; _dx < (x1 - x0); _dx++) { \
            const size_t _li   = _lrow + (size_t)_dx; \
            const size_t _i    = _crow + (size_t)_dx; \
            const float sr = pre_r[_li] * (TR); \
            const float sg = pre_g[_li] * (TG); \
            const float sb = pre_b[_li] * (TB); \
            const float _valid = (float)((sr + sg + sb) > 0.f); \
            const float _alpha = global_opacity \
                                 * (mask_data[_li] * INV_UCHAR_MAX_SINGLE) \
                                 * _valid; \
            const float dr = out_r[_i]; \
            const float dg = out_g[_i]; \
            const float db = out_b[_i]; \
            flis_blend_chroma_pixel(_alpha, \
                                    sr, sg, sb, dr, dg, db, \
                                    &out_r[_i], &out_g[_i], &out_b[_i]); \
        } \
    } \
} while (0)

/* -------------------------------------------------------------------------
 * FLIS_DISPATCH -- hoisted blend mode switch (Option B).
 *
 * Called once per layer.  Selects the masked or unmasked loop variant
 * (Option C) and then the blend-mode-specific inner loop.  Each case
 * expands to a fully inlined, branchless, SIMD-vectorisable loop body.
 *
 * TR, TG, TB: tint scalars passed through to the inner loops.
 * ------------------------------------------------------------------------- */
#define FLIS_DISPATCH(TR, TG, TB) \
do { \
    if (mask_data) { \
        switch (mode) { \
            case FLIS_BLEND_NORMAL: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, sr,sg,sb); break; \
            case FLIS_BLEND_MULTIPLY: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, \
                    sr*dr, sg*dg, sb*db); break; \
            case FLIS_BLEND_SCREEN: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, \
                    1.f-(1.f-sr)*(1.f-dr), \
                    1.f-(1.f-sg)*(1.f-dg), \
                    1.f-(1.f-sb)*(1.f-db)); break; \
            case FLIS_BLEND_OVERLAY: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, \
                    dr<0.5f?2.f*sr*dr:1.f-2.f*(1.f-sr)*(1.f-dr), \
                    dg<0.5f?2.f*sg*dg:1.f-2.f*(1.f-sg)*(1.f-dg), \
                    db<0.5f?2.f*sb*db:1.f-2.f*(1.f-sb)*(1.f-db)); break; \
            case FLIS_BLEND_SOFT_LIGHT: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, \
                    FLIS_SL(sr,dr), FLIS_SL(sg,dg), FLIS_SL(sb,db)); break; \
            case FLIS_BLEND_HARD_LIGHT: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, \
                    sr<0.5f?2.f*sr*dr:1.f-2.f*(1.f-sr)*(1.f-dr), \
                    sg<0.5f?2.f*sg*dg:1.f-2.f*(1.f-sg)*(1.f-dg), \
                    sb<0.5f?2.f*sb*db:1.f-2.f*(1.f-sb)*(1.f-db)); break; \
            case FLIS_BLEND_COLOR_DODGE: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, \
                    FLIS_DODGE(sr,dr), \
                    FLIS_DODGE(sg,dg), \
                    FLIS_DODGE(sb,db)); break; \
            case FLIS_BLEND_COLOR_BURN: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, \
                    FLIS_BURN(sr,dr), \
                    FLIS_BURN(sg,dg), \
                    FLIS_BURN(sb,db)); break; \
            case FLIS_BLEND_DARKEN: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, \
                    fminf(sr,dr), fminf(sg,dg), fminf(sb,db)); break; \
            case FLIS_BLEND_LIGHTEN: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, \
                    fmaxf(sr,dr), fmaxf(sg,dg), fmaxf(sb,db)); break; \
            case FLIS_BLEND_DIFFERENCE: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, \
                    fabsf(sr-dr), fabsf(sg-dg), fabsf(sb-db)); break; \
            case FLIS_BLEND_EXCLUSION: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, \
                    sr+dr-2.f*sr*dr, \
                    sg+dg-2.f*sg*dg, \
                    sb+db-2.f*sb*db); break; \
            default: \
                FLIS_BLEND_LOOP_MASKED(TR,TG,TB, sr,sg,sb); break; \
        } \
    } else { \
        switch (mode) { \
            case FLIS_BLEND_NORMAL: \
                FLIS_BLEND_LOOP(TR,TG,TB, sr,sg,sb); break; \
            case FLIS_BLEND_MULTIPLY: \
                FLIS_BLEND_LOOP(TR,TG,TB, \
                    sr*dr, sg*dg, sb*db); break; \
            case FLIS_BLEND_SCREEN: \
                FLIS_BLEND_LOOP(TR,TG,TB, \
                    1.f-(1.f-sr)*(1.f-dr), \
                    1.f-(1.f-sg)*(1.f-dg), \
                    1.f-(1.f-sb)*(1.f-db)); break; \
            case FLIS_BLEND_OVERLAY: \
                FLIS_BLEND_LOOP(TR,TG,TB, \
                    dr<0.5f?2.f*sr*dr:1.f-2.f*(1.f-sr)*(1.f-dr), \
                    dg<0.5f?2.f*sg*dg:1.f-2.f*(1.f-sg)*(1.f-dg), \
                    db<0.5f?2.f*sb*db:1.f-2.f*(1.f-sb)*(1.f-db)); break; \
            case FLIS_BLEND_SOFT_LIGHT: \
                FLIS_BLEND_LOOP(TR,TG,TB, \
                    FLIS_SL(sr,dr), FLIS_SL(sg,dg), FLIS_SL(sb,db)); break; \
            case FLIS_BLEND_HARD_LIGHT: \
                FLIS_BLEND_LOOP(TR,TG,TB, \
                    sr<0.5f?2.f*sr*dr:1.f-2.f*(1.f-sr)*(1.f-dr), \
                    sg<0.5f?2.f*sg*dg:1.f-2.f*(1.f-sg)*(1.f-dg), \
                    sb<0.5f?2.f*sb*db:1.f-2.f*(1.f-sb)*(1.f-db)); break; \
            case FLIS_BLEND_COLOR_DODGE: \
                FLIS_BLEND_LOOP(TR,TG,TB, \
                    FLIS_DODGE(sr,dr), \
                    FLIS_DODGE(sg,dg), \
                    FLIS_DODGE(sb,db)); break; \
            case FLIS_BLEND_COLOR_BURN: \
                FLIS_BLEND_LOOP(TR,TG,TB, \
                    FLIS_BURN(sr,dr), \
                    FLIS_BURN(sg,dg), \
                    FLIS_BURN(sb,db)); break; \
            case FLIS_BLEND_DARKEN: \
                FLIS_BLEND_LOOP(TR,TG,TB, \
                    fminf(sr,dr), fminf(sg,dg), fminf(sb,db)); break; \
            case FLIS_BLEND_LIGHTEN: \
                FLIS_BLEND_LOOP(TR,TG,TB, \
                    fmaxf(sr,dr), fmaxf(sg,dg), fmaxf(sb,db)); break; \
            case FLIS_BLEND_DIFFERENCE: \
                FLIS_BLEND_LOOP(TR,TG,TB, \
                    fabsf(sr-dr), fabsf(sg-dg), fabsf(sb-db)); break; \
            case FLIS_BLEND_EXCLUSION: \
                FLIS_BLEND_LOOP(TR,TG,TB, \
                    sr+dr-2.f*sr*dr, \
                    sg+dg-2.f*sg*dg, \
                    sb+db-2.f*sb*db); break; \
            default: \
                FLIS_BLEND_LOOP(TR,TG,TB, sr,sg,sb); break; \
        } \
    } \
} while (0)

/* =========================================================================
 * HSL helper functions -- kept as regular functions because HSL conversions
 * are not vectorisable per-channel and the function-call overhead is
 * negligible relative to the conversion cost.
 * ========================================================================= */

static void flis_rgb_to_hsl(float r, float g, float b,
                             float *h, float *s, float *l) {
    float maxc = fmaxf(r, fmaxf(g, b));
    float minc = fminf(r, fminf(g, b));
    float delta = maxc - minc;
    *l = (maxc + minc) * 0.5f;
    if (delta < 1e-6f) { *h = *s = 0.f; return; }
    *s = delta / (1.f - fabsf(2.f * *l - 1.f));
    if      (maxc == r) *h = fmodf((g - b) / delta, 6.f) / 6.f;
    else if (maxc == g) *h = ((b - r) / delta + 2.f)     / 6.f;
    else                *h = ((r - g) / delta + 4.f)     / 6.f;
    if (*h < 0.f) *h += 1.f;
}

static void flis_hsl_to_rgb(float h, float s, float l,
                             float *r, float *g, float *b) {
    float c = (1.f - fabsf(2.f * l - 1.f)) * s;
    float x = c * (1.f - fabsf(fmodf(h * 6.f, 2.f) - 1.f));
    float m = l - c * 0.5f;
    int   seg = (int)(h * 6.f) % 6;
    float r1, g1, b1;
    switch (seg) {
        case 0: r1=c; g1=x; b1=0.f; break;
        case 1: r1=x; g1=c; b1=0.f; break;
        case 2: r1=0.f; g1=c; b1=x; break;
        case 3: r1=0.f; g1=x; b1=c; break;
        case 4: r1=x; g1=0.f; b1=c; break;
        default:r1=c; g1=0.f; b1=x; break;
    }
    *r = r1 + m; *g = g1 + m; *b = b1 + m;
}

static inline float flis_lum(float r, float g, float b) {
    return 0.2126f * r + 0.7152f * g + 0.0722f * b;
}

static void flis_clip_color(float *r, float *g, float *b) {
    float l = flis_lum(*r, *g, *b);
    float n = fminf(*r, fminf(*g, *b));
    float x = fmaxf(*r, fmaxf(*g, *b));
    if (n < 0.f) {
        float d = l - n;
        *r = l + (*r - l) * l / d;
        *g = l + (*g - l) * l / d;
        *b = l + (*b - l) * l / d;
    }
    if (x > 1.f) {
        float d = x - l;
        *r = l + (*r - l) * (1.f - l) / d;
        *g = l + (*g - l) * (1.f - l) / d;
        *b = l + (*b - l) * (1.f - l) / d;
    }
}

/* Shortest-path hue interpolation on the colour circle [0, 1). */
static inline float flis_lerp_hue(float h0, float h1, float t) {
    float d = h1 - h0;
    if (d >  0.5f) d -= 1.f;
    if (d < -0.5f) d += 1.f;
    float h = h0 + d * t;
    if (h <  0.f) h += 1.f;
    if (h >= 1.f) h -= 1.f;
    return h;
}

/* HSL blend modes — opacity is interpolated in HSL component space so that
 * partial opacity gives a genuine mix of the HSL component rather than an
 * RGB-space blend of the fully-replaced result.  At opacity == 1.0 the
 * output is identical to the classical "replace one component" definition.
 *
 * The function writes the final composited pixel directly; the caller must
 * NOT apply a further FLIS_OVER afterwards. */
static void flis_blend_hsl_pixel(flis_blend_mode_t mode, float opacity,
                                 float sr, float sg, float sb,
                                 float dr, float dg, float db,
                                 float *or, float *og, float *ob) {
    float sh, ss, sl, dh, ds, dl;
    flis_rgb_to_hsl(sr, sg, sb, &sh, &ss, &sl);
    flis_rgb_to_hsl(dr, dg, db, &dh, &ds, &dl);
    float rh, rs, rl;
    switch (mode) {
        case FLIS_BLEND_HUE:
            rh = flis_lerp_hue(dh, sh, opacity); rs = ds;                    rl = dl; break;
        case FLIS_BLEND_SATURATION:
            rh = dh;                              rs = ds + (ss - ds) * opacity; rl = dl; break;
        case FLIS_BLEND_COLOR:
            rh = flis_lerp_hue(dh, sh, opacity); rs = ds + (ss - ds) * opacity; rl = dl; break;
        case FLIS_BLEND_LUMINOSITY:
            rh = dh;                              rs = ds;                    rl = dl + (sl - dl) * opacity; break;
        default:
            *or = dr; *og = dg; *ob = db; return;
    }
    flis_hsl_to_rgb(rh, rs, rl, or, og, ob);
    flis_clip_color(or, og, ob);
}

/* LRGB chroma-transfer blend mode (FLIS_BLEND_CHROMA).
 *
 * Takes the colour direction from the source (top) layer and scales it to
 * match the luminance of the destination canvas.  This is the classical
 * luma/chroma ratio method for LRGB composition in astrophotography:
 *   result = source_rgb × (dest_luma / source_luma)
 *
 * BT.709 luma coefficients are used on the perceptually-stretched data
 * as-is; formal linearisation is inappropriate here because the ICC profile
 * on stretched astrophotography data is a consistency marker only.
 *
 * Opacity is applied via FLIS_OVER in RGB space after the transfer. */
static void flis_blend_chroma_pixel(float opacity,
                                    float sr, float sg, float sb,
                                    float dr, float dg, float db,
                                    float *or, float *og, float *ob) {
    const float src_lum = flis_lum(sr, sg, sb);
    const float dst_lum = flis_lum(dr, dg, db);
    float br, bg, bb;
    if (src_lum < 1e-6f) {
        br = dr; bg = dg; bb = db;   /* no colour data — pass through dest */
    } else {
        const float scale = dst_lum / src_lum;
        br = FLIS_C01(sr * scale);
        bg = FLIS_C01(sg * scale);
        bb = FLIS_C01(sb * scale);
    }
    *or = FLIS_OVER(br, dr, opacity);
    *og = FLIS_OVER(bg, dg, opacity);
    *ob = FLIS_OVER(bb, db, opacity);
}

/* -------------------------------------------------------------------------
 * flis_render_layers  (public) / flis_composite_build  (private wrapper)
 *
 * flis_render_layers: composites a caller-supplied list of FLIS layers into
 * a newly allocated float-RGB fits*.  This is the shared rendering kernel
 * used by the display pipeline, merge-down, flatten, and the flat-FITS save
 * path.  The caller owns the returned fits and must clearfits()+free() it.
 *
 * flis_composite_build: thin wrapper that calls flis_render_layers on the
 * full com.uniq->layers list and installs the result as flis_display_composite.
 *
 * Output: float RGB [0,1], same dimensions as the base layer (canvas).
 * Layer data may be float or USHORT, mono or RGB.
 * Layer masks are 8-bit greyscale (LMASK type).
 * All blending is performed in linear [0,1] float.
 *
 * The blend mode switch is hoisted above pixel loops (Option B).
 * Masked and unmasked paths are separately compiled (Option C).
 * USHORT sources are converted to float once per layer, not per pixel.
 * The inner loops are tagged with omp simd for auto-vectorisation.
 * ------------------------------------------------------------------------- */
static fits *flis_render_layers_internal(GSList *layers, gboolean sub_composite) {
    if (!layers) return NULL;

    flis_layer_t *base = (flis_layer_t *)layers->data;
    if (!base || !base->fit) return NULL;

    const guint W  = base->fit->rx;
    const guint H  = base->fit->ry;
    const size_t N = (size_t)W * H;

    fits *out = calloc(1, sizeof(fits));
    if (!out) { PRINT_ALLOC_ERR; return NULL; }
    out->fdata = malloc(N * 3 * sizeof(float));
    if (!out->fdata) { PRINT_ALLOC_ERR; free(out); return NULL; }

    memset(out->fdata, 0, N * 3 * sizeof(float));

    out->fpdata[RLAYER] = out->fdata;
    out->fpdata[GLAYER] = out->fdata + N;
    out->fpdata[BLAYER] = out->fdata + N * 2;
    out->type        = DATA_FLOAT;
    out->bitpix      = FLOAT_IMG;
    out->orig_bitpix = FLOAT_IMG;
    out->naxis       = 3;
    out->naxes[0]    = W;
    out->naxes[1]    = H;
    out->naxes[2]    = 3;
    out->rx          = W;
    out->ry          = H;

    /* The composite is always 3-channel RGB float, even when every input
     * layer is mono.  Copying a *Gray* ICC profile (which any mono base
     * layer naturally has) onto an RGB buffer is wrong: the display
     * pipeline would then color-transform the 3-plane composite as if it
     * were a single-channel gray source, collapsing R/G/B to luminance
     * and erasing any per-channel content contributed by tinted or RGB
     * layers above the base.  Match the profile's channel count to the
     * composite's, falling back to the working RGB standard when the
     * base supplies a gray profile. */
    if (base->fit->icc_profile && base->fit->color_managed) {
        gboolean base_is_gray = (base->fit->naxes[2] == 1);
        if (!base_is_gray) {
            out->icc_profile = copyICCProfile(base->fit->icc_profile);
            color_manage(out, TRUE);
        } else if (com.icc.working_standard) {
            out->icc_profile = copyICCProfile(com.icc.working_standard);
            color_manage(out, TRUE);
        }
        /* If no working RGB profile is available, leave the composite
         * unmanaged — display path treats it as raw pixels, which is
         * correct for a freshly assembled RGB float buffer. */
    }

    float *out_r = out->fpdata[RLAYER];
    float *out_g = out->fpdata[GLAYER];
    float *out_b = out->fpdata[BLAYER];

    gboolean first_layer = TRUE;

    /* Pre-pass: build sub-composites for non-PASS_THROUGH groups.
     * These are composited separately and then blended using the group blend mode. */
    GHashTable *grp_composites = NULL;
    GHashTable *grp_trigger    = NULL;
    GHashTable *skip_in_main   = NULL;

    if (!sub_composite && com.uniq && com.uniq->groups) {
        for (GSList *_gg = com.uniq->groups; _gg; _gg = _gg->next) {
            flis_group_t *_grp = (flis_group_t *)_gg->data;
            if (!_grp || _grp->blend_mode == FLIS_BLEND_PASS_THROUGH) continue;
            if (!_grp->visible) {
                if (!skip_in_main)
                    skip_in_main = g_hash_table_new(g_direct_hash, g_direct_equal);
                GSList *_ms = flis_group_get_layers(_grp);
                for (GSList *_m = _ms; _m; _m = _m->next)
                    g_hash_table_insert(skip_in_main, _m->data, GINT_TO_POINTER(1));
                g_slist_free(_ms);
                continue;
            }

            GSList *_ms = flis_group_get_layers(_grp);
            fits *_sub = flis_render_layers_internal(_ms, TRUE);
            g_slist_free(_ms);
            if (!_sub) continue;

            if (!grp_composites) grp_composites = g_hash_table_new(g_direct_hash, g_direct_equal);
            if (!grp_trigger)    grp_trigger    = g_hash_table_new(g_direct_hash, g_direct_equal);
            if (!skip_in_main)   skip_in_main   = g_hash_table_new(g_direct_hash, g_direct_equal);

            g_hash_table_insert(grp_composites, GINT_TO_POINTER(_grp->item_id), _sub);

            /* Find first member in the sorted layers list to use as trigger */
            for (GSList *_nl = layers; _nl; _nl = _nl->next) {
                flis_layer_t *_tl = (flis_layer_t *)_nl->data;
                if (_tl && _tl->group_id == _grp->item_id) {
                    g_hash_table_insert(grp_trigger, _tl, GINT_TO_POINTER(_grp->item_id));
                    break;
                }
            }

            /* Mark all members for skipping */
            GSList *_ms2 = flis_group_get_layers(_grp);
            for (GSList *_m = _ms2; _m; _m = _m->next)
                g_hash_table_insert(skip_in_main, _m->data, GINT_TO_POINTER(1));
            g_slist_free(_ms2);
        }
    }

    for (GSList *node = layers; node; node = node->next) {
        flis_layer_t *lay = (flis_layer_t *)node->data;
        if (!lay || !lay->fit || !lay->visible) continue;

        fits              *lfit          = lay->fit;
        layermask_t       *lmask         = lay->lmask;
        const gboolean     mono          = (lfit->naxes[2] == 1);
        const gboolean     is_float      = (lfit->type == DATA_FLOAT);
        const gboolean     hsl_mode      = (lay->blend_mode == FLIS_BLEND_HUE        ||
                                            lay->blend_mode == FLIS_BLEND_SATURATION  ||
                                            lay->blend_mode == FLIS_BLEND_COLOR       ||
                                            lay->blend_mode == FLIS_BLEND_LUMINOSITY);
        const gboolean     chroma_mode   = (lay->blend_mode == FLIS_BLEND_CHROMA);
        flis_blend_mode_t  mode          = lay->blend_mode;
        /* Group compositing: skip layers that belong to a non-PASS_THROUGH group
         * (they are handled via the group sub-composite). */
        if (!sub_composite && skip_in_main &&
                g_hash_table_contains(skip_in_main, lay)) {
            /* Check if this layer is the trigger point for a group composite */
            if (grp_trigger) {
                gpointer _gid_ptr = g_hash_table_lookup(grp_trigger, lay);
                if (_gid_ptr) {
                    gint _gid = GPOINTER_TO_INT(_gid_ptr);
                    fits *_gsub = grp_composites
                                  ? (fits *)g_hash_table_lookup(grp_composites, GINT_TO_POINTER(_gid))
                                  : NULL;
                    flis_group_t *_gobj = flis_group_get_by_id(_gid);
                    if (_gsub && _gobj) {
                        /* Blend group sub-composite using the group's blend mode and opacity.
                         * Re-use the standard dispatch macros by shadowing local variables. */
                        const float *pre_r        = _gsub->fpdata[RLAYER];
                        const float *pre_g        = _gsub->fpdata[GLAYER];
                        const float *pre_b        = _gsub->fpdata[BLAYER];
                        const float  global_opacity = _gobj->opacity;
                        const flis_blend_mode_t mode = _gobj->blend_mode;
                        const gboolean hsl_mode2   = (mode == FLIS_BLEND_HUE ||
                                                      mode == FLIS_BLEND_SATURATION ||
                                                      mode == FLIS_BLEND_COLOR ||
                                                      mode == FLIS_BLEND_LUMINOSITY);
                        const gboolean chroma_mode2 = (mode == FLIS_BLEND_CHROMA);
                        const uint8_t *mask_data   = NULL;
                        const gint    ox = 0, oy = 0;
                        const gint    x0 = 0, x1 = (gint)W;
                        const gint    y0 = 0, y1 = (gint)H;
                        const guint   lW = W, lH = H;
                        (void)lH;  /* kept for symmetry with sparse-layer paths */
                        const float   tr = 1.f, tg = 1.f, tb = 1.f;
                        if (first_layer) {
                            FLIS_BASE_LOOP(tr, tg, tb);
                            first_layer = FALSE;
                        } else if (hsl_mode2) {
                            FLIS_HSL_LOOP(tr, tg, tb);
                        } else if (chroma_mode2) {
                            FLIS_CHROMA_LOOP(tr, tg, tb);
                        } else {
                            FLIS_DISPATCH(tr, tg, tb);
                        }
                    }
                }
            }
            continue;
        }

        /* Group modifiers: if this layer belongs to a hidden group, skip it;
         * if the group has reduced opacity, scale the effective opacity. */
        gfloat effective_opacity = lay->opacity;
        if (!sub_composite && lay->group_id != 0) {
            flis_group_t *grp = flis_group_get_by_id(lay->group_id);
            if (grp) {
                if (!grp->visible) continue;
                effective_opacity *= grp->opacity;
            }
        }
        const float        global_opacity = effective_opacity;

        /* Pre-compute source float pointers.
         *
         * For float source: use fpdata directly (zero allocation, zero copy).
         * For USHORT source: convert once to a temporary float buffer using a
         * vectorised parallel pass, then blend from that buffer.
         *
         * Mono sources use a single plane; all three pre_* pointers are set to
         * the same array.  Tint (tr/tg/tb) then differentiates the channels
         * during the blend loop via a scalar broadcast multiply.              */
        float        *float_buf = NULL;
        const float  *pre_r, *pre_g, *pre_b;

        if (is_float) {
            pre_r = lfit->fpdata[RLAYER];
            pre_g = mono ? lfit->fpdata[RLAYER] : lfit->fpdata[GLAYER];
            pre_b = mono ? lfit->fpdata[RLAYER] : lfit->fpdata[BLAYER];
        } else {
            const size_t n_planes = mono ? 1 : 3;
            float_buf = malloc(N * n_planes * sizeof(float));
            if (!float_buf) { PRINT_ALLOC_ERR; continue; }

            const WORD *sw[3] = {
                lfit->pdata[RLAYER],
                mono ? lfit->pdata[RLAYER] : lfit->pdata[GLAYER],
                mono ? lfit->pdata[RLAYER] : lfit->pdata[BLAYER]
            };
            for (size_t _c = 0; _c < n_planes; _c++) {
                float      *fw  = float_buf + _c * N;
                const WORD *src = sw[_c];
                FLIS_OMP_PAR_FOR_SIMD
                for (size_t _k = 0; _k < N; _k++)
                    fw[_k] = src[_k] * INV_USHRT_MAX_SINGLE;
            }

            pre_r = float_buf;
            pre_g = mono ? float_buf : float_buf + N;
            pre_b = mono ? float_buf : float_buf + N * 2;
        }

        const float tr = (mono && lay->has_tint) ? (float)lay->layer_tint.r : 1.f;
        const float tg = (mono && lay->has_tint) ? (float)lay->layer_tint.g : 1.f;
        const float tb = (mono && lay->has_tint) ? (float)lay->layer_tint.b : 1.f;

        const uint8_t *mask_data = (lmask && lay->lmask_active) ? (const uint8_t *)lmask->data : NULL;

        /* Intersection of the layer rectangle with the canvas.
         * position_y is in FITS convention (origin bottom-left), convert to
         * top-down canvas row: oy = H - position_y - lH.
         * For full-canvas: ox=0, oy=0, lW=W, lH=H → x0=0,x1=W,y0=0,y1=H. */
        const gint  ox  = lay->position_x;
        const guint lW  = (guint)lfit->rx;
        const guint lH  = (guint)lfit->ry;
        const gint  oy  = (gint)H - lay->position_y - (gint)lH;
        const gint  x0  = MAX(0,      ox);
        const gint  x1  = MIN((gint)W, ox + (gint)lW);
        const gint  y0  = MAX(0,      oy);
        const gint  y1  = MIN((gint)H, oy + (gint)lH);

        if (x0 >= x1 || y0 >= y1) {
            /* Layer entirely outside the canvas */
            if (first_layer) first_layer = FALSE;
            free(float_buf);
            continue;
        }

        if (first_layer) {
            FLIS_BASE_LOOP(tr, tg, tb);
            first_layer = FALSE;
        } else if (hsl_mode) {
            if (mask_data) {
                FLIS_HSL_LOOP_MASKED(tr, tg, tb);
            } else {
                FLIS_HSL_LOOP(tr, tg, tb);
            }
        } else if (chroma_mode) {
            if (mask_data) {
                FLIS_CHROMA_LOOP_MASKED(tr, tg, tb);
            } else {
                FLIS_CHROMA_LOOP(tr, tg, tb);
            }
        } else {
            FLIS_DISPATCH(tr, tg, tb);
        }

        free(float_buf);  /* no-op for float sources (float_buf == NULL) */
    }

    /* Cleanup group compositing tables */
    if (grp_composites) {
        GHashTableIter _it;
        gpointer _k, _v;
        g_hash_table_iter_init(&_it, grp_composites);
        while (g_hash_table_iter_next(&_it, &_k, &_v)) {
            fits *_gsub = (fits *)_v;
            if (_gsub) { free(_gsub->fdata); free(_gsub); }
        }
        g_hash_table_destroy(grp_composites);
    }
    if (grp_trigger)  g_hash_table_destroy(grp_trigger);
    if (skip_in_main) g_hash_table_destroy(skip_in_main);

    return out;
}
