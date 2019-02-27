
/* about FFTW: http://www.fftw.org/fftw3_doc/Introduction.html
 * The standard FFTW distribution works most efficiently for arrays whose size
 * can be factored into small primes (2, 3, 5, and 7), and otherwise it uses a
 * slower general-purpose routine.
 * Maybe we need to constrain the zone size to these numbers, that's not too
 * many that do not register (11 13 17 19 23 31 37 41 43 47 53 59 61 67 71 ...)
 * wait, do they need to decompose down to only 2 3 5 and 7? If so the list is
 * longer.
 *
 * Applying the phase correlation method to a pair of images produces a third
 * image which contains a single peak. The location of this peak corresponds to
 * the relative translation between the images.
 * The Fourier-Mellin transform extends phase correlation to handle images
 * transformed by both translation and rotation.
 */
int the_multipoint_dft_registration(struct mpr_args *args) {
	int zone_idx, nb_zones, frame;
	int abort = 0;
	int retval = 0;
	stacking_zone *zone;
	fftw_complex **ref, **in, **out, **convol;
	fftw_plan *fplan, *bplan;	// forward and backward plans
	regdata *regparam = args->seq->regparam[args->layer];
	if (!regparam) return -1;

	nb_zones = get_number_of_zones();
	if (nb_zones < 1) {
		fprintf(stderr, "cannot do the multi-point registration if no zone is defined\n");
		return -1;
	}
	ref = malloc(nb_zones * sizeof(fftw_complex*));
	in = malloc(nb_zones * sizeof(fftw_complex*));
	out = malloc(nb_zones * sizeof(fftw_complex*));
	convol = malloc(nb_zones * sizeof(fftw_complex*));
	fplan = malloc(nb_zones * sizeof(fftw_plan));
	bplan = malloc(nb_zones * sizeof(fftw_plan));

	/* reading zones in the reference image, single threaded init */
	fprintf(stdout, "loading reference zones\n");
	//gchar *wisdom_file = get_configdir_file_path(FFTW_WISDOM_FILE);
	// TO BE FIXED IF USED AGAIN SOME DAY
	gchar *wisdom_file = "/tmp/"FFTW_WISDOM_FILE;
	if (wisdom_file) {
		if (fftw_import_wisdom_from_filename(wisdom_file))
			fprintf(stdout, "FFTW wisdom restored\n");
	}

#ifdef USE_SEQUENCE_REF_IMAGE
	fits reffits = { 0 };
	if (seq_read_frame(args->seq, args->seq->reference_image, &reffits))
		return -1;
	fprintf(stdout, "Using the reference image from sequence (%d) instead "
			"of the stacked reference\n", args->seq->reference_image);
#endif

	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		zone = &com.stacking_zones[zone_idx];
		int side = get_side(zone);
		int nb_pixels = side * side;
		// allocate aligned DFT buffers with fftw_malloc
		ref[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
		in[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
		out[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
		convol[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);

		start_timer();
		fplan[zone_idx] = fftw_plan_dft_2d(side, side, ref[zone_idx], out[zone_idx], FFTW_FORWARD, FFTW_PATIENT);
		bplan[zone_idx] = fftw_plan_dft_2d(side, side, convol[zone_idx], out[zone_idx], FFTW_BACKWARD, FFTW_PATIENT);
		fprintf(stdout, "plan %d creation time: %ld microsec\n", zone_idx,
				stop_timer_elapsed_mus());
		// the plan creation time should be long only the first time a zone size is used

#ifdef USE_SEQUENCE_REF_IMAGE
		copy_image_zone_to_fftw(&reffits, zone, ref[zone_idx], args->layer);
#else
		// use the refimage stacked from the best of sequence
		copy_image_zone_to_fftw(&refimage, zone, ref[zone_idx], args->layer);
#endif
		fftw_execute_dft(fplan[zone_idx], ref[zone_idx], in[zone_idx]);
	}
#ifdef USE_SEQUENCE_REF_IMAGE
	clearfits(&reffits);
#endif

	// in out and convol can probably be freed here, or even allocated once and reused
	// for the above loop
	if (wisdom_file) {
		if (fftw_export_wisdom_to_filename(wisdom_file))
			fprintf(stdout, "FFTW wisdom saved\n");
		g_free(wisdom_file);
	}

	fftw_complex **zones = calloc(nb_zones, sizeof(fftw_complex*));
	fftw_complex **out2 = calloc(nb_zones, sizeof(fftw_complex*));
	fftw_complex **convol2 = calloc(nb_zones, sizeof(fftw_complex*));

	/* for each image, we read the zones with global shift and register them */
	/* for sequences that require demosaicing, seq_read_frame is usually the longest
	 * operation in this loop, so parallelizing the zone alignment does not help
	 * much, either with fftw3_threads or with OpenMP in the zone loop below. The
	 * best way to speed things up is probably to execute this loop in parallel.
	 * The fftw plan execution is thread-safe so it should not be a problem. */
	for (frame = 0; frame < args->seq->number; frame++) {
		fits fit = { 0 };
		int i;
		if (abort) continue;
		/* TODO: this filtering is required for global mode, not for local */
		/*if (!args->filtering_criterion(args->seq, args->layer,
					frame, args->filtering_parameter) || abort)
			continue;*/

		if (seq_read_frame(args->seq, frame, &fit)) {
			abort = 1;
			continue;
		}
		fprintf(stdout, "aligning zones for image %d\n", frame);

		/* reading zones in the reference image */
#ifdef _OPENMP
//#pragma omp parallel for num_threads(com.max_thread) private(zone_idx, zone) schedule(static) if((args->seq->type == SEQ_REGULAR && fits_is_reentrant()) || args->seq->type == SEQ_SER)
#endif
		for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
			if (abort) continue;
			zone = &com.stacking_zones[zone_idx];
			int side = get_side(zone);
			int nb_pixels = side * side;
			stacking_zone shifted_zone = { .centre =
				{ .x = zone->centre.x - regparam[frame].shiftx,
					.y = zone->centre.y + regparam[frame].shifty },
				.half_side = zone->half_side };

			if (!zones[zone_idx]) {	// keep across images
			       zones[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
			       out2[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
			       convol2[zone_idx] = fftw_malloc(sizeof(fftw_complex) * nb_pixels);
			}

			copy_image_zone_to_fftw(&fit, &shifted_zone, zones[zone_idx],
					args->layer);

			// forward transformation zone -> out2
			fftw_execute_dft(fplan[zone_idx], zones[zone_idx], out2[zone_idx]);

			// compute a new fourier domain image defined as the fourier
			// representation of the reference image on this zone * conj(out2)
			for (i = 0; i < nb_pixels; i++) {
				convol2[zone_idx][i] = in[zone_idx][i] * conj(out2[zone_idx][i]);
			}

			// backward transformation of this new image to out2
			fftw_execute_dft(bplan[zone_idx], convol2[zone_idx], out2[zone_idx]);

			// searching for the real part peak in out2, which is the shift
			// between the reference image and this image in this zone
			int shift = 0;
			for (i = 1; i < nb_pixels; ++i) {
				if (creal(out2[zone_idx][i]) > creal(out2[zone_idx][shift])) {
					shift = i;
				}
			}
			int shifty = shift / side;
			int shiftx = shift % side;
			if (shifty > zone->half_side) {
				shifty -= side;
			}
			if (shiftx > zone->half_side) {
				shiftx -= side;
			}

			/* shitfs for this image and this zone is (shiftx, shifty) + the global shifts */
			fprintf(stdout, "frame %d, zone %d adjustment shifts: %d,%d\n", frame, zone_idx, shiftx, shifty);
			zone->mpregparam[frame].x = (double)shiftx + regparam[frame].shiftx;
			zone->mpregparam[frame].y = (double)shifty - regparam[frame].shifty;
		}

		clearfits(&fit);

		// TODO: this will require a fix similar to what's done in the generic
		// function, especially if it's executed in parallel
		set_progress_bar_data(NULL, frame/(double)args->seq->number);
	}

cleaning_all:
	for (zone_idx = 0; zone_idx < nb_zones; zone_idx++) {
		fftw_free(zones[zone_idx]);
		fftw_free(out2[zone_idx]);
		fftw_free(convol2[zone_idx]);

		fftw_destroy_plan(fplan[zone_idx]);
		fftw_destroy_plan(bplan[zone_idx]);
		fftw_free(in[zone_idx]);
		fftw_free(out[zone_idx]);
		fftw_free(ref[zone_idx]);
		fftw_free(convol[zone_idx]);
	}
	free(zones); free(out2); free(convol2);
	free(fplan); free(bplan);
	free(in); free(out); free(ref); free(convol);

	return retval;
}


// copy the image zone into a double buffer, it remains upside-down
static int copy_image_zone_to_fftw(fits *fit, const stacking_zone *zone, fftw_complex *dest,
		int layer) {
	int side = get_side(zone);
	// start coordinates on the displayed image, but images are read upside-down
	int startx = round_to_int(zone->centre.x - zone->half_side);
	int starty = round_to_int(zone->centre.y - zone->half_side);
	if (startx < 0 || startx >= fit->rx - side || starty < 0 || starty >= fit->ry - side) {
		/* this zone is partly outside the image, I don't think there's
		 * much we can do for it, it just has to be ignored for this
		 * image for the stacking. */
		return -1;
	}

	WORD *from = fit->pdata[layer] + (fit->ry - starty - side - 1) * fit->rx + startx;
	int stridefrom = fit->rx - side;
	int i, j;

	for (i = 0; i < side; ++i) {
		for (j = 0; j < side; ++j) {
			*dest++ = (double)*from++;
		}
		from += stridefrom;
	}
	return 0;
}


