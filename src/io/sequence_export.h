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
#ifndef SRC_IO_SEQUENCE_EXPORT_H_
#define SRC_IO_SEQUENCE_EXPORT_H_

#include "core/siril.h"
#include "io/sequence.h"

/* same order as in the combo box 'combo_export_preset' */
typedef enum {
	PRESET_RATIO,
	PRESET_FREE,
	PRESET_FULLHD,
	PRESET_4K,
	PRESET_8K
} export_presets;

struct exportseq_args {
	sequence *seq;
	seq_image_filter filtering_criterion;
	double filtering_parameter;

	char *basename;
	export_format output;
	gboolean normalize;

	gboolean tiff_compression; /* compression of TIFF files */

	int film_fps;		/* has to be int for avi or ffmpeg */
	int film_quality;	/* [1, 5], for mp4 and webm */

	gboolean resample;
	int32_t dest_width, dest_height;

	gboolean crop;
	rectangle crop_area;
};

/* Start the sequence export in a new thread. Returns TRUE on success. */
gboolean sequence_export_start(struct exportseq_args *args);

#endif /* SRC_IO_SEQUENCE_EXPORT_H_ */
