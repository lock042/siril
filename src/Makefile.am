## Process this file with automake to produce Makefile.in

bin_PROGRAMS = siril siril-cli

AM_CPPFLAGS = \
	-DPACKAGE_DATA_DIR='"$(datadir)/${PACKAGE}"' \
	-DPACKAGE_DOC_DIR='"${datarootdir}/doc/${PACKAGE}"' \
	-DLOCALEDIR='"${localedir}"' \
	${OPENCV_CFLAGS} \
	$(FFMPEG_CFLAGS) \
	$(GIO_CFLAGS) \
	$(JSON_GLIB_CFLAGS) \
	$(GTK_CFLAGS) \
	${GTK_MAC_CFLAGS} \
	$(LIBTIFF_CFLAGS) \
	$(LIBPNG_CFLAGS) \
	$(FFTW_CFLAGS) \
	$(CFITSIO_CFLAGS) \
	$(GSL_CFLAGS) \
	$(FFMS2_CFLAGS) \
	$(WCSLIB_CFLAGS) \
	$(LIBCONFIG_CFLAGS) \
	$(LCMS2_CFLAGS)  \
	$(LIBGIT_CFLAGS) \
	$(GTKSOURCEVIEW_CFLAGS) \
	-I$(top_builddir)/subprojects/kplot


siril_resource.c: siril_resource.h
	$(AM_V_GEN)
	glib-compile-resources \
	--sourcedir=$(top_builddir) \
	--sourcedir=$(sourcedir) \
	--target=$@ \
	--generate-source \
	--c-name siril_resource \
	$(srcdir)/siril_resource.xml

siril_common_SOURCES: siril_resource.h

siril_common_SOURCES = \
	core/siril.h \
	core/proto.h \
	algos/anscombe.c \
	algos/anscombe.h \
	algos/astrometry_solver.c \
	algos/astrometry_solver.h \
	algos/background_extraction.c \
	algos/background_extraction.h \
	algos/ccd-inspector.c \
	algos/ccd-inspector.h \
	algos/colors.c \
	algos/colors.h \
	algos/comparison_stars.c \
	algos/comparison_stars.h \
	algos/Def_Math.h \
	algos/Def_Wavelet.h \
	algos/Def_Mem.h \
	algos/extraction.c \
	algos/extraction.h \
	algos/demosaicing.c \
	algos/demosaicing.h \
	algos/fitting.c \
	algos/fitting.h \
	algos/fix_xtrans_af.c \
	algos/fix_xtrans_af.h \
	algos/geometry.c \
	algos/geometry.h \
	algos/io_wave.c \
	algos/median_fast.c \
	algos/noise.c \
	algos/noise.h \
	algos/pave.c \
	algos/photometry.c \
	algos/photometry.h \
	algos/photometric_cc.c \
	algos/photometric_cc.h \
	algos/PSF.c \
	algos/PSF.h \
	algos/quality.c \
	algos/quality_float.c \
	algos/quality.h \
	algos/quantize.c \
	algos/reconstr.c \
	algos/search_objects.c \
	algos/search_objects.h \
	algos/siril_random.cpp \
	algos/siril_wcs.h \
	algos/siril_wcs.c \
	algos/sorting.c \
	algos/sorting.h \
	algos/spcc.c \
	algos/spcc.h \
	algos/star_finder.c \
	algos/star_finder.h \
	algos/statistics.c \
	algos/statistics.h \
	algos/statistics_float.c \
	algos/statistics_float.h \
	algos/transform.c \
	algos/lcms_acceleration/fast_16_tethra.c \
	algos/lcms_acceleration/fast_8_curves.c \
	algos/lcms_acceleration/fast_8_matsh.c \
	algos/lcms_acceleration/fast_8_matsh_sse.c \
	algos/lcms_acceleration/fast_8_tethra.c \
	algos/lcms_acceleration/fast_float_15bits.c \
	algos/lcms_acceleration/fast_float_15mats.c \
	algos/lcms_acceleration/fast_float_cmyk.c \
	algos/lcms_acceleration/fast_float_curves.c \
	algos/lcms_acceleration/fast_float_internal.h \
	algos/lcms_acceleration/fast_float_lab.c \
	algos/lcms_acceleration/fast_float_matsh.c \
	algos/lcms_acceleration/fast_float_separate.c \
	algos/lcms_acceleration/fast_float_sup.c \
	algos/lcms_acceleration/fast_float_tethra.c \
	algos/lcms_acceleration/threaded_core.c \
	algos/lcms_acceleration/threaded_main.c \
	algos/lcms_acceleration/threaded_scheduler.c \
	algos/lcms_acceleration/threaded_split.c \
	compositing/align_rgb.c \
	compositing/align_rgb.h \
	compositing/compositing.c \
	compositing/compositing.h \
	compositing/filters.c \
	compositing/filters.h \
	core/arithm.c \
	core/arithm.h \
	core/command.c \
	core/command.h \
	core/command_def.h \
	core/command_line_processor.c \
	core/command_line_processor.h \
	core/command_list.h \
	core/icc_profile.c \
	core/icc_profile.h \
	core/initfile.c \
	core/initfile.h \
	core/OS_utils.c \
	core/OS_utils.h \
	core/pipe.c \
	core/pipe.h \
	core/preprocess.c \
	core/preprocess.h \
	core/processing.c \
	core/sequence_filtering.c \
	core/sequence_filtering.h \
	core/settings.c \
	core/settings.h \
	core/signals.c \
	core/signals.h \
	core/siril.c \
	core/siril_actions.c \
	core/siril_actions.h \
	core/siril_alloc.c \
	core/siril_alloc.h \
	core/siril_app_dirs.c \
	core/siril_app_dirs.h \
	core/siril_cmd_help.c \
	core/siril_cmd_help.h \
	core/siril_date.c \
	core/siril_date.h \
	core/siril_language.c \
	core/siril_language.h \
	core/siril_log.c \
	core/siril_log.h \
	core/siril_networking.c \
	core/siril_networking.h \
	core/siril_spawn.c \
	core/siril_spawn.h \
	core/siril_update.c \
	core/siril_update.h \
	core/siril_world_cs.c \
	core/siril_world_cs.h \
	core/undo.c \
	core/undo.h \
	core/utils.c \
	drizzle/cdrizzlebox.c \
	drizzle/cdrizzlebox.h \
	drizzle/cdrizzlemap.c \
	drizzzle/cdrizzlemap.h \
	drizzle/cdrizzleutil.c \
	drizzle/cdrizzleutil.h \
	drizzle/driz_portability.h \
	filters/asinh.c \
	filters/asinh.h \
	filters/banding.c \
	filters/banding.h \
	filters/clahe.c \
	filters/clahe.h \
	filters/cosmetic_correction.c \
	filters/cosmetic_correction.h \
	filters/curve_transform.c \
	filters/curve_transform.h \
	filters/deconvolution.c \
	filters/epf.c \
	filters/epf.h \
	filters/fft.c \
	filters/fft.h \
	filters/linear_match.c \
	filters/linear_match.h \
	filters/median.c \
	filters/median.h \
	filters/mtf.c \
	filters/mtf.h \
	filters/ght.c \
	filters/ght.h \
	filters/rgradient.c \
	filters/rgradient.h \
	filters/saturation.c \
	filters/saturation.h \
	filters/scnr.c \
	filters/scnr.h \
	filters/starnet.c \
	filters/starnet.h \
	filters/synthstar.c \
	filters/synthstar.h \
	filters/unpurple.c \
	filters/unpurple.h \
	filters/wavelets.c \
	filters/wavelets.h \
	filters/da3d/DA3D.cpp \
	filters/da3d/DA3D.hpp \
	filters/da3d/DftPatch.hpp \
	filters/da3d/Image.hpp \
	filters/da3d/Utils.cpp \
	filters/da3d/Utils.hpp \
	filters/da3d/WeightMap.cpp \
	filters/da3d/WeightMap.hpp \
	filters/deconvolution/angleSet.cpp \
	filters/deconvolution/angleSet.hpp \
	filters/deconvolution/better_than_std.hpp \
	filters/deconvolution/deconvolution.h \
	filters/deconvolution/deconvolve.cpp \
	filters/deconvolution/deconvolve.hpp \
	filters/deconvolution/edgetaper.hpp \
	filters/deconvolution/estimate_kernel.cpp \
	filters/deconvolution/estimate_kernel.hpp \
	filters/deconvolution/fft.hpp \
	filters/deconvolution/fftw_allocator.hpp \
	filters/deconvolution/gf_estimate.cpp \
	filters/deconvolution/gf_estimate.hpp \
	filters/deconvolution/image_expr.hpp \
	filters/deconvolution/image.hpp \
	filters/deconvolution/labeling.hpp \
	filters/deconvolution/optimization.hpp \
	filters/deconvolution/options.hpp \
	filters/deconvolution/utils.hpp \
	filters/deconvolution/vec2.hpp \
	filters/nlbayes/call_nlbayes.cpp \
	filters/nlbayes/call_nlbayes.h \
	filters/nlbayes/LibImages.cpp \
	filters/nlbayes/LibImages.h \
	filters/nlbayes/LibMatrix.cpp \
	filters/nlbayes/NlBayes.cpp \
	filters/nlbayes/NlBayes.h \
	filters/nlbayes/Utilities.cpp \
	filters/nlbayes/Utilities.h \
	gui/about_dialog.c \
	gui/about_dialog.h \
	gui/annotate.c \
	gui/annotations_pref.c \
	gui/astrometry_solver.c \
	gui/denoisegui.c \
	gui/background_extraction.c \
	gui/callbacks.c \
	gui/callbacks.h \
	gui/compstars.c \
	gui/compstars.h \
	gui/conversion.c \
	gui/conversion.h \
	gui/colors.c \
	gui/colors.h \
	gui/curves.c \
	gui/curves.h \
	gui/cut.c \
	gui/dialog_preview.c \
	gui/dialog_preview.h \
	gui/dialogs.c \
	gui/dialogs.h \
	gui/documentation.c \
	gui/documentation.h \
	gui/epf.c \
	gui/epf.h \
	gui/fix_xtrans_af.c \
	gui/fix_xtrans_af.h \
	gui/git_gui.c \
	gui/histogram.c \
	gui/histogram.h \
	gui/keywords_tree.c \
	gui/keywords_tree.h \
	gui/icc_profile.c \
	gui/icc_profile.h \
	gui/image_display.c \
	gui/image_display.h \
	gui/image_interactions.c \
	gui/image_interactions.h \
	gui/linear_match.c \
	gui/menu_gray_geometry.c \
	gui/menu_gray_geometry.h \
	gui/menu_gray_selection.c \
	gui/merge_cfa.c \
	gui/message_dialog.c \
	gui/message_dialog.h \
	gui/mouse_action_functions.c \
	gui/mouse_action_functions.h \
	gui/mouse_action_preferences.c \
	gui/newdeconv.c \
	gui/nina_light_curve.c \
	gui/nina_light_curve.h \
	gui/open_dialog.c \
	gui/open_dialog.h \
	gui/photometric_cc.c \
	gui/photometric_cc.h \
	gui/plot.c \
	gui/plot.h \
	gui/preferences.c \
	gui/preferences.h \
	gui/progress_and_log.c \
	gui/progress_and_log.h \
	gui/PSF_list.h \
	gui/PSF_list.c \
	gui/python_gui.h \
	gui/python_gui.c \
	gui/registration.h \
	gui/registration.c \
	gui/registration_preview.c \
	gui/remixer.c \
	gui/remixer.h \
	gui/save_dialog.c \
	gui/save_dialog.h \
	gui/script_menu.c \
	gui/script_menu.h \
	gui/sequence_list.c \
	gui/sequence_list.h \
	gui/single_image.c \
	gui/single_image.h \
	gui/siril_css.c \
	gui/siril_css.h \
	gui/siril_intro.c \
	gui/siril_intro.h \
	gui/siril_preview.c \
	gui/siril_preview.h \
	gui/siril_plot.c \
	gui/siril_plot.h \
	gui/siril-window.c \
	gui/siril-window.h \
	gui/splashscreen.c \
	gui/splashscreen.h \
	gui/split_cfa.c \
	gui/stacking.c \
	gui/stacking.h \
	gui/star_finder.c \
	gui/starnetgui.c \
	gui/starnetgui.h \
	gui/statistics_list.c \
	gui/unpurple.c \
	gui/unpurple.h \
	gui/user_polygons.c \
	gui/user_polygons.h \
	gui/utils.c \
	gui/utils.h \
	io/annotation_catalogues.c \
	io/annotation_catalogues.h \
	io/aavso_extended.c \
	io/aavso_extended.h \
	io/Astro-TIFF.c \
	io/Astro-TIFF.h \
	io/kstars/binfile.h \
	io/kstars/byteorder.h \
	io/conversion.c \
	io/conversion.h \
	io/films.c \
	io/films.h \
	io/fits_keywords.c \
	io/fits_keywords.h \
	io/fits_sequence.c \
	io/fits_sequence.h \
	io/FITS_symlink.c \
	io/FITS_symlink.h \
	io/healpix/healpix.cpp \
	io/healpix/healpix_cat.h \
	io/image_format_fits.c \
	io/image_format_fits.h \
	io/image_formats_internal.c \
	io/image_formats_libraries.c \
	io/local_catalogues.c \
	io/local_catalogues.h \
	io/mp4_output.c \
	io/mp4_output.h \
	io/path_parse.c \
	io/path_parse.h \
	io/remote_catalogues.c \
	io/remote_catalogues.h \
	io/seqfile.c \
	io/sequence.c \
	io/sequence.h \
	io/sequence_export.c \
	io/seqwriter.h \
	io/seqwriter.c \
	io/ser.c \
	io/ser.h \
	io/single_image.c \
	io/single_image.h \
	io/siril_catalogues.c \
	io/siril_catalogues.h \
	io/siril_git.c \
	io/siril_git.h \
	io/siril_plot.c \
	io/siril_plot.h \
	io/siril_pythoncommands.h \
	io/siril_pythoncommands.c \
	io/siril_pythonmodule.h \
	io/siril_pythonmodule.c \
	io/spcc_json.c \
	livestacking/livestacking.c \
	livestacking/livestacking.h \
	livestacking/gui.c \
	livestacking/gui.h \
	pixelMath/pixel_math_runner.c \
	pixelMath/pixel_math_runner.h \
	pixelMath/tinyexpr.c \
	pixelMath/tinyexpr.h \
	registration/applyreg.c \
	registration/3stars.c \
	registration/comet.c \
	registration/global.c \
	registration/matching/match.c \
	registration/matching/atpmatch.c \
	registration/matching/misc.c \
	registration/matching/degtorad.h \
	registration/matching/apply_match.c \
	registration/astrometric.c \
	registration/distorsion.c \
	registration/distorsion.h \
	registration/registration.c \
	registration/registration.h \
	registration/shift_methods.c \
	stacking/blending.c \
	stacking/blending.h \
	stacking/median_and_mean.c \
	stacking/rejection_float.c \
	stacking/normalization.c \
	stacking/siril_fit_linear.c \
	stacking/siril_fit_linear.h \
	stacking/stacking.c \
	stacking/stacking.h \
	stacking/stackminmax.c \
	stacking/sum.c \
	stacking/sum.h \
	stacking/upscaling.c \
	rt/rt_algo.cc \
	rt/gauss.cc \
	rt/boxblur.cc

if PLATFORM_WIN32
libdbghelp = -ldbghelp
endif
if !HAVE_YYJSON
YYJSON_BUILDDIR = $(abs_top_builddir)/subprojects/yyjson/build
YYJSON_SRCDIR = $(abs_top_srcdir)/subprojects/yyjson

$(YYJSON_BUILDDIR)/Makefile: $(YYJSON_SRCDIR)/CMakeLists.txt
	$(MKDIR_P) $(YYJSON_BUILDDIR)
	cd $(YYJSON_BUILDDIR) && $(CMAKE) $(CMAKE_GENERATOR) \
		-DCMAKE_INSTALL_PREFIX:PATH=$(prefix) \
		-DCMAKE_C_FLAGS="$(CFLAGS)" \
		$(YYJSON_SRCDIR)

yyjson: $(YYJSON_BUILDDIR)/Makefile
	$(MAKE) -C $(YYJSON_BUILDDIR)

BUILT_SOURCES = yyjson
endif

if !HAVE_YYJSON
AM_CPPFLAGS += -I$(YYJSON_SRCDIR)/src
YYJSON_LIBS = $(YYJSON_BUILDDIR)/libyyjson.a
else
AM_CPPFLAGS += $(LIBYYJSON_CFLAGS)
YYJSON_LIBS = $(LIBYYJSON_LIBS)
endif

if HAVE_LIBRTPROCESS
AM_CPPFLAGS += -DHAVE_LIBRTPROCESS \
						$(LIBRTPROCESS_CFLAGS)
librtprocesslibs = $(LIBRTPROCESS_LIBS)
else
AM_CPPFLAGS += -I$(top_srcdir)/subprojects/librtprocess/src/include
librtprocesslibs = $(top_builddir)/subprojects/librtprocess/build/src/librtprocess.a
endif

AM_CPPFLAGS += -I$(top_srcdir)/subprojects/htmesh
libhtmeshlibs = $(top_builddir)/subprojects/htmesh/build/libhtmesh.a

AM_CPPFLAGS += -I$(top_srcdir)/subprojects/healpix_cxx
libhealpix_cxxlibs = $(top_builddir)/subprojects/healpix_cxx/build/libhealpix_cxx.a

siril_common_LDADD = \
	$(top_builddir)/subprojects/kplot/libkplot.a \
	$(FFMPEG_LIBS) \
	$(GTK_LIBS) \
	${GTK_MAC_LIBS} \
	$(LIBRAW_LIBS) \
	$(LIBHEIF_LIBS) \
	$(LIBTIFF_LIBS) \
	$(LIBPNG_LIBS) \
	$(LIBXISF_LIBS) \
	$(FFTW_LIBS) \
	$(CFITSIO_LIBS) \
	$(LCMS2_LIBS) \
	$(GSL_LIBS) \
	$(JSON_GLIB_LIBS) \
	$(FFMS2_LIBS) \
	$(LIBCONFIG_LIBS) \
	$(LIBCURL) \
	$(WCSLIB_LIBS) \
	$(EXTRA_LIBS) \
	$(LIBATOMIC_LIBS) \
	$(libdbghelp) \
	$(librtprocesslibs) \
	$(LIBGIT_LIBS) \
	$(YYJSON_LIBS) \
	$(GTKSOURCEVIEW_LIBS)

siril_common_SOURCES += \
	algos/demosaicing_rtp.cpp \
	io/avi_pipp/pipp_avi_write.cpp \
	io/avi_pipp/pipp_avi_write.h \
	io/avi_pipp/pipp_avi_write_dib.cpp \
	io/avi_pipp/pipp_avi_write_dib.h \
	io/avi_pipp/pipp_utf8.h \
	io/avi_pipp/pipp_buffer.cpp \
	io/avi_pipp/pipp_buffer.h \
	io/avi_pipp/pipp_video_write.h \
	io/avi_pipp/avi_writer.cpp \
	io/avi_pipp/avi_writer.h \
	io/SirilXISFReader.cpp \
	io/SirilXISFReader.h \
	io/SirilXISFWraper.h \
	opencv/opencv.cpp \
	opencv/opencv.h \
	opencv/guidedfilter.cpp \
	opencv/guidedfilter.h \
	opencv/kombat/kombat.cpp \
	opencv/kombat/kombat.h \
	core/exif.cpp \
	core/exif.h \
	io/kstars/htmesh_wrapper.cpp \
	io/kstars/htmesh_wrapper.h

siril_common_LDADD += $(OPENCV_LIBS) \
		      $(EXIV2_LIBS) \
		      $(libhtmeshlibs) \
		      $(libhealpix_cxxlibs)

siril_SOURCES = $(siril_common_SOURCES) siril_resource.c siril_resource.h main.c
siril_cli_SOURCES = $(siril_common_SOURCES) main-cli.c
siril_LDADD = $(siril_common_LDADD)
siril_cli_LDADD = $(siril_common_LDADD)

if PLATFORM_WIN32
AM_LDFLAGS += -Wl,--export-all-symbols
if !SIRIL_DEBUG
AM_LDFLAGS += -mwindows -lbcrypt
endif
else
if PLATFORM_OSX
AM_LDFLAGS += -rdynamic -Wl,-framework,Security
else
AM_LDFLAGS += -rdynamic
endif
endif

siril_LINK = $(CXXLINK)
