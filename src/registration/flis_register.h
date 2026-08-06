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
 */

#ifndef SRC_REGISTRATION_FLIS_REGISTER_H_
#define SRC_REGISTRATION_FLIS_REGISTER_H_

#include "core/siril.h"
#include "opencv/opencv.h"
#include "registration/registration.h"

/* Symbolic identifiers for the registration methods that the FLIS
 * layer-registration dialog and command expose.  Maps to the existing
 * registration_function family in registration.h. */
typedef enum {
	FLIS_REG_GLOBAL = 0,   /* register_star_alignment — single-pass global star alignment */
	FLIS_REG_2PASS,        /* register_multi_step_global — two-pass global */
	FLIS_REG_DFT,          /* register_shift_dft — DFT cross-correlation, REQUIRES_SQUARED_SELECTION */
	FLIS_REG_KOMBAT,       /* register_kombat — pattern matching, REQUIRES_ANY_SELECTION */
	FLIS_REG_N_METHODS
} flis_reg_method_id;

/**
 * flis_register_resolve_method:
 *
 * Returns the registration_function pointer for the given id, along
 * with its selection requirement and the appropriate transformation
 * type for the method (SHIFT_TRANSFORMATION for DFT/KOMBAT, HOMOGRAPHY
 * for the global methods).  Returns NULL if @id is out of range.
 */
registration_function flis_register_resolve_method(flis_reg_method_id id,
                                                    selection_type *out_sel,
                                                    transformation_type *out_tx);

/**
 * flis_register_layers:
 * @ref_lay:        reference layer (NULL → first layer of @target_layers)
 * @target_layers:  layers to register (NULL → every layer in com.uniq->layers)
 * @method:         registration_function from registration.h.  When NULL
 *                   defaults to register_star_alignment (single-pass
 *                   global star alignment) — the safest method to run
 *                   on a freshly-loaded FLIS without any selection.
 * @sel_req:        selection requirement for @method (see registration.h's
 *                   selection_type).  The function validates com.selection
 *                   against this requirement and refuses with an error
 *                   when the requirement isn't met (e.g. DFT shift needs
 *                   a square selection).
 * @tx_type:        transformation type to record on the registration
 *                   args (SHIFT_TRANSFORMATION for DFT/KOMBAT;
 *                   HOMOGRAPHY_TRANSFORMATION for the global methods).
 *                   When @method is NULL, defaults to HOMOGRAPHY.
 * @interpolation:  resampling interpolation used by register_apply_reg.
 * @clamp:          if TRUE, clamp the cubic / Lanczos passes to suppress
 *                   ringing on bright pixels.
 * @method_id:      which flis_reg_method_id @method came from.  The
 *                   function pointer alone is ambiguous (GLOBAL and 2PASS
 *                   both resolve to register_multi_step_global), and the
 *                   NDE record has to name the SETTING the user chose so a
 *                   replay that must re-solve can reproduce it.  Pass
 *                   FLIS_REG_N_METHODS when @method is NULL / unknown.
 *
 * Builds an internal sequence from @target_layers and runs @method
 * followed by register_apply_reg with FRAMING_MAX.  Layer pixel data
 * is resampled in place; canvas offsets are written back to each
 * layer's position_x / position_y.  Captures ONE joint NDE record
 * (flis.register, nde_joint.h) spanning every participant.  Synchronous —
 * wrap in a worker thread for async dispatch.  Returns 0 on success,
 * non-zero on failure (logged).
 */
int flis_register_layers(flis_layer_t *ref_lay,
                         GSList       *target_layers,
                         registration_function method,
                         selection_type        sel_req,
                         transformation_type   tx_type,
                         opencv_interpolation interpolation,
                         gboolean       clamp,
                         flis_reg_method_id method_id);

/* One participant's registration RESULT: the transform the warp needs, the
 * size it warps into and where the layer lands on the canvas.  Produced by
 * both the live capture and the replay-time re-solve, so the two paths
 * cannot drift in what they mean by "the answer". */
typedef struct {
	double H[9];      /* row-major, post-framing */
	int    out_rx, out_ry;
	int    pos_x, pos_y;
} flis_reg_solution;

/**
 * flis_register_solve:
 *
 * Re-solve a registration OFFLINE, against caller-supplied pixels rather
 * than the live layer stack — the L2 half of the flis.register NDE replay
 * (nde_joint.h).  @fits_in are @n resolved participant images in the
 * record's participant order and @ref_index says which is the reference;
 * they are NOT modified.  @out receives @n solutions in the same order.
 *
 * Runs the same two steps the live path runs (the registration method, then
 * register_apply_reg with FRAMING_MAX) on an internal sequence built over
 * COPIES of @fits_in, then reads the transforms back out.
 *
 * @canvas_index is the participant index of the document's base layer, or
 * -1 when the base layer is not a participant; together with @ref_pos_x /
 * @ref_pos_y (where the reference layer sat before the registration) it
 * reproduces the live path's choice of canvas origin exactly.
 *
 * Returns 0 on success.  Non-zero means the solve did not run — the caller
 * is expected to fall back to whatever it recorded (L3).
 */
int flis_register_solve(fits **fits_in, int n, int ref_index,
                        registration_function method,
                        selection_type sel_req,
                        transformation_type tx_type,
                        rectangle selection,
                        int canvas_index,
                        int ref_pos_x, int ref_pos_y,
                        flis_reg_solution *out);

#endif /* SRC_REGISTRATION_FLIS_REGISTER_H_ */
