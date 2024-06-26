/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "algos/statistics.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"

#include "arithm.h"

/*****************************************************************************
 *       S I R I L      A R I T H M E T I C      O P E R A T I O N S         *
 ****************************************************************************/

static int soper_ushort_to_ushort(fits *a, float scalar, image_operator oper) {
	if (!a) return 1;
	if (!a->data) return 1;
	WORD *data = NULL;
	size_t i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	if (!n) return 1;
	data = a->data;
	float invnorm = (a->bitpix == BYTE_IMG) ? INV_UCHAR_MAX_SINGLE : 1.f;
	if (oper == OPER_DIV) {
		scalar = 1.0f / scalar;
		oper = OPER_MUL;
	}
	if (oper == OPER_MUL) scalar *= invnorm;

	switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				float pixel = ushort_to_float_bitpix(a, data[i]);
				if (invnorm == 1.f) {
					data[i] = float_to_ushort_range(pixel + scalar);
				} else {
					data[i] = float_to_ushort_range(invnorm * (pixel + scalar));
				}
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				float pixel = ushort_to_float_bitpix(a, data[i]);
				if (invnorm == 1.f) {
					data[i] = float_to_ushort_range(pixel - scalar);
				} else {
					data[i] = float_to_ushort_range(invnorm * (pixel - scalar));
				}
			}
			break;
		case OPER_MUL:
		default:
			for (i = 0; i < n; ++i) {
				data[i] = roundf_to_WORD((float)data[i] * scalar);
			}
			break;
	}
	invalidate_stats_from_fit(a);
	return 0;
}

static int soper_ushort_to_float(fits *a, float scalar, image_operator oper) {
	if (!a) return 1;
	if (!a->data) return 1;
	WORD *data = NULL;
	float *result= NULL;
	size_t i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	if (!n) return 1;
	data = a->data;
	result = malloc(n * sizeof(float));
	if (!result) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	if (oper == OPER_DIV) {
		scalar = 1.0f / scalar;
		oper = OPER_MUL;
	}

	switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				float pixel = ushort_to_float_bitpix(a, data[i]);
				result[i] = pixel + scalar;
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				float pixel = ushort_to_float_bitpix(a, data[i]);
				result[i] = pixel - scalar;
			}
			break;
		case OPER_MUL:
		default:
			for (i = 0; i < n; ++i) {
				float pixel = ushort_to_float_bitpix(a, data[i]);
				result[i] = pixel * scalar;
			}
			break;
	}
	fit_replace_buffer(a, result, DATA_FLOAT);
	return 0;
}

int soper_unscaled_div_ushort_to_float(fits *a, int scalar) {
	if (!a) return 1;
	if (!a->data) return 1;
	size_t i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	if (!n) return 1;
	WORD *data = a->data;
	float *result = malloc(n * sizeof(float));
	if (!result) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	float operand = 1.0f / (float)scalar;
	for (i = 0; i < n; ++i) {
		result[i] = data[i] * operand;
	}
	fit_replace_buffer(a, result, DATA_FLOAT);
	return 0;
}

static int soper_float(fits *a, float scalar, image_operator oper) {
	if (!a) return 1;
	if (!a->fdata) return 1;
	float *data = a->fdata;
	size_t i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	if (!n) return 1;
	if (oper == OPER_DIV) {
		scalar = 1.0f / scalar;
		oper = OPER_MUL;
	}

	switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				data[i] = data[i] + scalar;
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				data[i] = data[i] - scalar;
			}
			break;
		case OPER_MUL:
		default:
			for (i = 0; i < n; ++i) {
				data[i] = data[i] * scalar;
			}
			break;
	}
	invalidate_stats_from_fit(a);
	return 0;
}

/* equivalent to (map simple_operation a), with simple_operation being
 * (lambda (pixel) (oper pixel scalar))
 * scalar is applied on images with data in [0, 1]
 */
int soper(fits *a, float scalar, image_operator oper, gboolean conv_to_float) {
	if (oper == OPER_DIV && scalar == 0.f) {
		siril_log_message(_("Cannot divide by zero, aborting."));
		return 1;
	}
	if (a->type == DATA_USHORT) {
		if (conv_to_float)
			return soper_ushort_to_float(a, scalar, oper);
		return soper_ushort_to_ushort(a, scalar, oper);
	}
	if (a->type == DATA_FLOAT)
		return soper_float(a, scalar, oper);
	return 1;
}

static int imoper_to_ushort(fits *a, fits *b, image_operator oper, float factor) {
	if (!a) return 1;
	if (!a->data) return 1;
	if (!b) return 1;
	size_t i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	if (!n) return 1;
	if (memcmp(a->naxes, b->naxes, sizeof a->naxes)) {
		siril_log_color_message(_("Images must have same dimensions.\n"), "red");
		return 1;
	}

	if (b->type == DATA_USHORT) {
		if (!b->data) return 1;
		WORD *abuf = a->data, *bbuf = b->data;
		if (oper == OPER_DIV) {
			for (i = 0; i < n; ++i) {
				if (bbuf[i] == 0)
					abuf[i] = 0;
				else {
					float aval = (float) abuf[i];
					float bval = (float) bbuf[i];
					if (factor != 1.0f)
						abuf[i] = roundf_to_WORD(factor * (aval / bval));
					else
						abuf[i] = roundf_to_WORD(aval / bval);
				}
			}
		} else if (oper == OPER_MUL) {
			for (i = 0; i < n; ++i) {
				if (bbuf[i] == 0)
					abuf[i] = 0;
				else {
					float aval = (float) abuf[i];
					float bval = (float) bbuf[i];
					if (factor != 1.0f)
						abuf[i] = roundf_to_WORD(factor * (aval * bval));
					else
						abuf[i] = roundf_to_WORD(aval * bval);
				}
			}
		} else {
			for (i = 0; i < n; ++i) {
				int aval = (int) abuf[i];
				int bval = (int) bbuf[i];
				switch (oper) {
				case OPER_ADD:
					abuf[i] = truncate_to_WORD(aval + bval);
					break;
				case OPER_SUB:
					abuf[i] = truncate_to_WORD(aval - bval);
					break;
				default:	// OPER_MUL, OPER_DIV handled above
					break;
				}

				if (factor != 1.0f)
					abuf[i] = roundf_to_WORD(factor * (float) abuf[i]);
			}
		}
	} else if (b->type == DATA_FLOAT) {
		if (!b->fdata) return 1;
		WORD *abuf = a->data;
		float *bbuf = b->fdata;
		float norm = (a->bitpix == BYTE_IMG) ? UCHAR_MAX_SINGLE : USHRT_MAX_SINGLE;

		if (oper == OPER_DIV) {
			for (i = 0; i < n; ++i) {
				if (bbuf[i] == 0.f)
					abuf[i] = 0;
				else {
					float aval = (float) abuf[i];
					float bval = bbuf[i] * norm;
					if (factor != 1.0f)
						abuf[i] = roundf_to_WORD(factor * (aval / bval));
					else
						abuf[i] = roundf_to_WORD(aval / bval);
				}
			}
		} else if (oper == OPER_MUL) {
			for (i = 0; i < n; ++i) {
				if (bbuf[i] == 0.f)
					abuf[i] = 0;
				else {
					float aval = (float) abuf[i];
					float bval = bbuf[i] * norm;
					if (factor != 1.0f)
						abuf[i] = roundf_to_WORD(factor * (aval * bval));
					else
						abuf[i] = roundf_to_WORD(aval * bval);
				}
			}
		} else {
			for (i = 0; i < n; ++i) {
				int aval = (int) abuf[i];
				int bval = roundf_to_int(bbuf[i] * norm);
				switch (oper) {
				case OPER_ADD:
					abuf[i] = truncate_to_WORD(aval + bval);
					break;
				case OPER_SUB:
					abuf[i] = truncate_to_WORD(aval - bval);
					break;
				default:	// OPER_MUL and OPER_DIV handled above
					break;
				}
				if (factor != 1.0f)
					abuf[i] = roundf_to_WORD(factor * (float) abuf[i]);
			}
		}
	}
	invalidate_stats_from_fit(a);
	a->neg_ratio = 0.0f;
	return 0;
}

int imoper_to_float(fits *a, fits *b, image_operator oper, float factor) {
	if (!a) return 1;
	if (!b) return 1;
	size_t n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	if (!n) return 1;
	float *result = NULL;

	if (memcmp(a->naxes, b->naxes, sizeof a->naxes)) {
		siril_log_color_message(_("Images must have same dimensions.\n"), "red");
		return 1;
	}

	if (a->type == DATA_FLOAT) {
		if (!a->fdata) return 1;
		result = a->fdata;
	}
	else if (a->type == DATA_USHORT) {
		result = malloc(n * sizeof(float));
		if (!result) {
			PRINT_ALLOC_ERR;
			return 1;
		}
	}
	else return 1;

	if (b->type == DATA_USHORT && (!b->data)) { free(result); return 1; }
	if (b->type == DATA_FLOAT && (!b->fdata)) { free(result); return 1; }
	if (!(b->type == DATA_FLOAT || b->type == DATA_USHORT)) { free(result); return 1; }

	size_t nb_negative = 0;
	for (size_t i = 0; i < n; ++i) {
		float aval = a->type == DATA_USHORT ? ushort_to_float_bitpix(a, a->data[i]) : a->fdata[i];
		float bval = b->type == DATA_USHORT ? ushort_to_float_bitpix(b, b->data[i]) : b->fdata[i];
		switch (oper) {
			case OPER_ADD:
				result[i] = aval + bval;
				break;
			case OPER_SUB:
				result[i] = aval - bval;
				break;
			case OPER_MUL:
				result[i] = aval * bval;
				break;
			case OPER_DIV:
				if (bval == 0.0f)
					result[i] = 0.0f;
				else result[i] = aval / bval;
		}
		if (factor != 1.0f)
			result[i] *= factor;
		if (result[i] > 1.0f)	// should we truncate by default?
			result[i] = 1.0f;
		if (result[i] < -1.0f)	// Here we want to clip all garbages pixels.
			result[i] = 0.0f;
		if (result[i] < 0.0f)
			nb_negative++;
	}
	a->neg_ratio = (float)((double)nb_negative / n);
	if (a->type == DATA_USHORT)
		fit_replace_buffer(a, result, DATA_FLOAT);
	else invalidate_stats_from_fit(a);
	return 0;
}

/* applies operation of image a with image b, for all their layers:
 * a = factor * a oper b
 * returns 0 on success */
static int imoper_with_factor(fits *a, fits *b, image_operator oper, float factor, gboolean allow_32bits) {
	// ushort result can only be forced when both input images are ushort
	if (allow_32bits && !com.pref.force_16bit)
		return imoper_to_float(a, b, oper, factor);
	else {
		if (a->type == DATA_USHORT)
			return imoper_to_ushort(a, b, oper, factor);
		siril_log_color_message(_("Image operations can only be kept 16 bits if first input images are 16 bits. Aborting.\n"), "red");
	}
	return 1;
}

int imoper(fits *a, fits *b, image_operator oper, gboolean allow_32bits) {
	return imoper_with_factor(a, b, oper, 1.0f, allow_32bits);
}

/* a = coef * a / b
 * a is expected as USHORT, returned as FLOAT. b is expected as USHORT.
 * If overflow, siril_fdiv returns 1*/
int siril_fdiv(fits *a, fits *b, float coef, gboolean allow_32bits) {
	return imoper_with_factor(a, b, OPER_DIV, coef, allow_32bits);
}

// a = max(a, b)
int addmax(fits *a, fits *b) {
	size_t i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];

	if (memcmp(a->naxes, b->naxes, sizeof a->naxes)) {
		siril_log_color_message(_("Images must have same dimensions.\n"), "red");
		return 1;
	}
	if (a->type != b->type) {
		siril_log_color_message(_("Images must have same data type.\n"), "red");
		return 1;
	}
	g_assert(a->naxes[2] == 1 || a->naxes[2] == 3);

	if (a->type == DATA_USHORT) {
		WORD *abuf = a->data, *bbuf = b->data;
		for (i = 0; i < n; ++i) {
			if (bbuf[i] > abuf[i])
				abuf[i] = bbuf[i];
		}
	} else {
		float *abuf = a->fdata, *bbuf = b->fdata;
		for (i = 0; i < n; ++i) {
			if (bbuf[i] > abuf[i])
				abuf[i] = bbuf[i];
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

void rgbblend(blend_data *data, float* r, float* g, float* b, float m_CB) {
	float *out[3] = { r, g, b };

	// Compute the maximum values using fmaxf for better vectorization
	float sfmax = fmaxf(fmaxf(data->sf[0], data->sf[1]), data->sf[2]);
	float tfmax = fmaxf(fmaxf(data->tf[0], data->tf[1]), data->tf[2]);

	// Compute the difference
	float d = sfmax - tfmax;

	// Calculate k based on conditions
	float k = (tfmax + m_CB * d > 1.f) ? ((d != 0.f) ? fminf(m_CB, (1.f - tfmax) / d) : m_CB) : m_CB;

	// Loop over channels
	for (size_t chan = 0; chan < 3; chan++) {
		// Compute the output using the precomputed k
		*out[chan] = data->do_channel[chan] ? (1.f - k) * data->tf[chan] + k * data->sf[chan] : *out[chan];
	}
}
