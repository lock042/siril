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
 *
 * Builds an internal sequence from @target_layers and runs @method
 * followed by register_apply_reg with FRAMING_MAX.  Layer pixel data
 * is resampled in place; canvas offsets are written back to each
 * layer's position_x / position_y.  Synchronous — wrap in a worker
 * thread for async dispatch.  Returns 0 on success, non-zero on
 * failure (logged).
 */
int flis_register_layers(flis_layer_t *ref_lay,
                         GSList       *target_layers,
                         registration_function method,
                         selection_type        sel_req,
                         transformation_type   tx_type,
                         opencv_interpolation interpolation,
                         gboolean       clamp);

#endif /* SRC_REGISTRATION_FLIS_REGISTER_H_ */
