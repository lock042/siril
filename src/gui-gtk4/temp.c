/* run the SPCC using the existing star list of the image from the provided file */
gpointer spectrophotometric_cc_standalone(gpointer p) {
	struct photometric_cc_data *args = (struct photometric_cc_data *)p;
	float kw[3];
	coeff bg[3];
	float mins[3];
	float maxs[3];
	int norm_channel;

	/* run peaker to measure FWHM of the image to adjust photometry settings */
	args->fwhm = measure_image_FWHM(args->fit, -1);
	if (args->fwhm <= 0.0f) {
		siril_log_message(_("Error computing FWHM for photometry settings adjustment\n"));
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}

	double ra, dec;
	center2wcs(args->fit, &ra, &dec);
	double resolution = get_wcs_image_resolution(args->fit);
	uint64_t sqr_radius = ((uint64_t) gfit->rx * (uint64_t) gfit->rx + (uint64_t) gfit->ry * (uint64_t) gfit->ry) / 4;
	double radius = resolution * sqrt((double)sqr_radius);	// in degrees
	double mag = args->mag_mode == LIMIT_MAG_ABSOLUTE ?
		args->magnitude_arg : compute_mag_limit_from_fov(radius * 2.0);
	if (args->mag_mode == LIMIT_MAG_AUTO_WITH_OFFSET)
		mag += args->magnitude_arg;

	// Initialize filters if required
	if (!spcc_filters_initialized) {
		init_spcc_filters();
		spcc_filters_initialized = TRUE;
	}

	if (!isrgb(args->fit)) {
		siril_log_message(_("Photometric color correction will do nothing for monochrome images\n"));
		return 0;
	}

	rectangle *bkg_sel = NULL;
	if (!args->bg_auto)
		bkg_sel = &(args->bg_area);

	/* we use the median of each channel to sort them by level and select
	 * the reference channel expressed in terms of order of middle median value */
	if (get_stats_coefficients(args->fit, bkg_sel, bg, mins, maxs, &norm_channel)) {
		siril_log_message(_("failed to compute statistics on image, aborting\n"));
		free(args);
		return GINT_TO_POINTER(1);
	}
	siril_log_message(_("Normalizing on %s channel.\n"), (norm_channel == 0) ? _("red") : ((norm_channel == 1) ? _("green") : _("blue")));

	/* set photometry parameters to values adapted to the image */
	struct phot_config backup = com.pref.phot_set;
	com.pref.phot_set.force_radius = FALSE;
	com.pref.phot_set.inner = max(7.0, 3.0 * args->fwhm);
	com.pref.phot_set.outer = com.pref.phot_set.inner + 10;
	siril_log_message(_("Photometry radii set to %.1f for inner and %.1f for outer\n"),
			com.pref.phot_set.inner, com.pref.phot_set.outer);

	set_progress_bar_data(_("Spectrophotometric color calibration in progress..."), PROGRESS_RESET);

	//	Calculate filter responses at 2nm spacing from 380nm to 700nm
	//	Either look up the responses (for OSC) or calculate as a product of filter_resp * QE
	xpsampled response[3] = { { xpsampled_wl, { 0.0 } }, { xpsampled_wl, { 0.0 } }, { xpsampled_wl, { 0.0 } } };
	for (int chan = 0 ; chan < 3 ; chan++) {
		get_spectrum_from_ui(&response[chan], chan);
	}

	// Calculate effective primaries (and get the desired white point while we're
	// at it) for later
	cmsCIExyYTRIPLE primaries;
	primaries.Red = xpsampled_to_xyY(&response[RED], CMF_1931);
	primaries.Green = xpsampled_to_xyY(&response[GREEN], CMF_1931);
	primaries.Blue = xpsampled_to_xyY(&response[BLUE], CMF_1931);
	cmsCIExyY spcc_whitepoint;
//	spcc_whitepoint = get_spcc_whitepoint_from_UI(); TODO
	memcpy(&spcc_whitepoint, &Whitepoint_D58, sizeof(cmsCIExyY));

	//	Do a Gaia DR3 cone search and retrieve the reference stars CSV and the
	//	xp_sampled FITS
	gchar *datalink_path = NULL;
	if (!args->ref_stars) {
		args->ref_stars = calloc(1, sizeof(siril_catalogue));
		args->ref_stars->cat_index = CAT_GAIADR3_DIRECT;
	}

	// Carry out the PCC catalog query
	siril_gaiadr3_datalink_query(args->ref_stars, XP_SAMPLED, &datalink_path);

	// Project the catalog using WCS
	siril_catalog_project_with_WCS(args->ref_stars, args->fit, TRUE, FALSE);

	struct phot_config *ps = phot_set_adjusted_for_image(args->fit);
	siril_debug_print("aperture: %2.1f%s\tinner: %2.1f\touter: %2.1f\n", ps->aperture, ps->force_radius?"":" (auto)", ps->inner, ps->outer);
	gint ngood = 0, progress = 0;
	gint errors[PSF_ERR_MAX_VALUE] = { 0 };
	//	Obtain spectrophotometric white balance
	for (int i = 0 ; i < args->ref_stars->nbitems ; i++) {
		// Obtain sampled spectrum between 380-700nm from the xp_sampled FITS
		// Each call to get_xpsampled() opens a new fptr so it is safe to call
		// from multiple threads
		xpsampled xp_sampled = { xpsampled_wl, { 0.0 } };
		get_xpsampled(&xp_sampled, datalink_path, i);

		float ref_flux[3] = { 0.f };
		xpsampled flux_expected = { xpsampled_wl, { 0.0 } };
		// Don't parallelize this loop
		for (int chan = 0 ; chan < 3 ; chan++) {
			multiply_xpsampled(&flux_expected, &response[chan], &xp_sampled);
			ref_flux[i] = integrate_xpsampled(&flux_expected);
		}
		// Obtain the expected ratios of R, G and B
		float max_flux = max(max(ref_flux[0], ref_flux[1]), ref_flux[2]);
		for (int chan = 0 ; chan < 3 ; chan++)
			ref_flux[chan] /= max_flux;

		// Measure actual flux by running peaker around ref star coordinates
		rectangle area = { 0 };
		float flux[3] = { 0.f };
		if (make_selection_around_a_star(args->ref_stars->cat_items[i].x, args->ref_stars->cat_items[i].y, &area, args->fit)) {
			siril_debug_print("star %d is outside image or too close to border\n", i);
			g_atomic_int_inc(errors+PSF_ERR_OUT_OF_WINDOW);
			continue;
		}
		gboolean no_phot = FALSE;
		psf_error error = PSF_NO_ERR;
		for (int chan = 0; chan < 3 && !no_phot; chan ++) {
			psf_star *photometry = psf_get_minimisation(args->fit, chan, &area, TRUE, ps, FALSE, com.pref.starfinder_conf.profile, &error);
			g_atomic_int_inc(errors+error);
			if (!photometry || !photometry->phot_is_valid || error != PSF_NO_ERR)
				no_phot = TRUE;
			else flux[chan] = powf(10.f, -0.4f * (float) photometry->mag);
			if (photometry)
				free_psf(photometry);
		}
		if (no_phot) {
			siril_debug_print("photometry failed for star %d, error %d\n", i, error);
			continue;
		}

		// Compare ref_flux with actual flux
		// (rest is as per PCC)
	}

	free(ps);
	// Apply white balance as per PCC
	int ret = 0;
	if (!ret) {
		ret = apply_photometric_color_correction(args->fit, kw, bg, mins, maxs, norm_channel);
	} else {
		set_progress_bar_data(_("Spectrophotometric Color Calibration failed"), PROGRESS_DONE);
	}

	// Transform image from source profile to working profile using
	// INTENT_ABSOLUTE_COLORIMETRIC
	cmsToneCurve *curve[3], *tonecurve;
	tonecurve = cmsBuildGamma (NULL, 1.00);
	curve[0] = curve[1] = curve[2] = tonecurve;
	cmsHPROFILE source_profile = cmsCreateRGBProfile(&spcc_whitepoint, &primaries, curve);
	cmsFreeToneCurve(tonecurve);

	cmsHPROFILE profile = NULL;
	if (args->fit->icc_profile) {
		profile = copyICCProfile(args->fit->icc_profile);
		if (!fit_icc_is_linear(args->fit)) {
			siril_log_color_message(_("Image color space is nonlinear. It is recommended to "
					"apply photometric color calibration to linear images.\n"), "salmon");
		}
	} else {
		profile = siril_color_profile_linear_from_color_profile(com.icc.working_standard);
		args->fit->icc_profile = copyICCProfile(profile);
		color_manage(args->fit, (args->fit->icc_profile != NULL));
	}

	// Create transform from source profile to image-based profile, clean up the profiles
	// TODO...

	// Apply transform and clean up
	// TODO...

	com.pref.phot_set = backup;
	free(args);
	return GINT_TO_POINTER(ret);
}
