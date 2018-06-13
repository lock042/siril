/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2018 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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

/* This is the main file for the multi-point planetary processing */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <glib.h>
#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/processing.h"
#include "core/proto.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "stacking/sum.h"
#include "io/sequence.h"
#include "gui/progress_and_log.h"

//static fits refimage;

static gpointer sequence_analysis_thread_func(gpointer p);
static gboolean end_reference_image_stacking(gpointer p);

/* First step: running the global registration and building the reference image */
void on_planetary_analysis_clicked(GtkButton *button, gpointer user_data) {
	if (!sequence_is_loaded()) return;
	GtkComboBox *cbbt_layers = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "comboboxreglayer"));

	struct registration_args *reg_args;
	reg_args = calloc(1, sizeof(struct registration_args));
	reg_args->func = register_cog;
	reg_args->seq = &com.seq;
	reg_args->reference_image = sequence_find_refimage(&com.seq);
	reg_args->process_all_frames = FALSE;
	reg_args->follow_star = FALSE;
	reg_args->matchSelection = FALSE;
	reg_args->translation_only = TRUE;
	reg_args->x2upscale = FALSE;
	reg_args->layer = gtk_combo_box_get_active(cbbt_layers);
	reg_args->run_in_thread = TRUE;
	reg_args->load_new_sequence = FALSE;

	char *msg = siril_log_color_message(_("Starting sequence analysis\n"), "red");
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);
	set_cursor_waiting(TRUE);

	start_in_new_thread(sequence_analysis_thread_func, reg_args);
}

static gpointer sequence_analysis_thread_func(gpointer p) {
	struct registration_args *reg_args = (struct registration_args *) p;

	reg_args->retval = reg_args->func(reg_args);

	if (reg_args->retval)
		return GINT_TO_POINTER(-1);

	if (reg_args->seq->reference_image == -1) {
		// set new reference image: should we do it all the time?
		// also done in generated sequence in global.c
		reg_args->seq->reference_image = sequence_find_refimage(reg_args->seq);
	}
	writeseqfile(reg_args->seq);

	/* sequence quality evaluation is over, now we generate the reference
	 * image by taking a few of the best images */

	struct stacking_args *stack_args = calloc(1, sizeof(struct stacking_args));
	stack_args->method = stack_summing_generic;
	stack_args->seq = reg_args->seq;
	stack_args->filtering_criterion = stack_filter_quality;
	stack_args->filtering_parameter = 0.9;
	stack_args->nb_images_to_stack = compute_nb_filtered_images(stack_args);
	stack_args->image_indices = malloc(stack_args->nb_images_to_stack * sizeof(int));
	stack_fill_list_of_unfiltered_images(stack_args);
	snprintf(stack_args->description, sizeof stack_args->description,
			_("Creating the reference image for the multi-point planetary processing"));
	char *output_filename = malloc(strlen(reg_args->seq->seqname) + 20);
	sprintf(output_filename, "%s%sreference%s", reg_args->seq->seqname,
			ends_with(com.seq.seqname, "_") ?  "" :
			(ends_with(com.seq.seqname, "-") ? "" : "_"), com.ext);
	stack_args->output_filename = (const char*)output_filename;
	stack_args->output_overwrite = TRUE;
	stack_args->reglayer = reg_args->layer;
	gettimeofday(&stack_args->t_start, NULL);
	free(reg_args);

	stack_args->retval = stack_summing_generic(stack_args);
	if (stack_args->retval)
		return GINT_TO_POINTER(-1);
	/* the reference image is now in gfit */

	siril_add_idle(end_stacking, stack_args);
	siril_add_idle(end_reference_image_stacking, stack_args);
	
	return NULL;
}

static gboolean end_reference_image_stacking(gpointer p) {
	struct stacking_args *args = (struct stacking_args *)p;

	free(args);
	return FALSE;
}
