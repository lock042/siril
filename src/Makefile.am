## Process this file with automake to produce Makefile.in

# siril-cli is always built; siril (GUI) is built only when --enable-gtk
# is in effect (the default).  When --disable-gtk is used, GTK_CFLAGS,
# GTKSOURCEVIEW_CFLAGS, GTK_MAC_CFLAGS and their *_LIBS counterparts are
# all empty (PKG_CHECK_MODULES skipped them), so leaving the variable
# references in AM_CPPFLAGS / *_LDADD is harmless.

# siril is only built when --enable-gtk is in effect; EXTRA_PROGRAMS keeps
# automake aware of it so siril_SOURCES are tracked for dist purposes.
bin_PROGRAMS = siril-cli
EXTRA_PROGRAMS = siril
if ENABLE_GTK
bin_PROGRAMS += siril
endif

AM_CPPFLAGS = \
	-DPACKAGE_DATA_DIR='"$(datadir)/${PACKAGE}"' \
	-DPACKAGE_DOC_DIR='"${datarootdir}/doc/${PACKAGE}"' \
	-DLOCALEDIR='"${localedir}"' \
	${OPENCV_CFLAGS} \
	$(FFMPEG_CFLAGS) \
	$(GIO_CFLAGS) \
	$(GLIB_CFLAGS) \
	$(JSON_GLIB_CFLAGS) \
	$(CAIRO_CFLAGS) \
	$(PANGO_CFLAGS) \
	$(PANGOCAIRO_CFLAGS) \
	$(GDK_PIXBUF_CFLAGS) \
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
	$(SQLITE_CFLAGS) \
	$(LIBJXL_CFLAGS) \
	$(LIBJXL_THREADS_CFLAGS) \
	$(GTKSOURCEVIEW_CFLAGS) \
	$(X11_CFLAGS) \
	-I$(top_builddir)/subprojects/kplot

# The resource XML points at gui-gtk4/uifiles/* and gui-gtk4/siril.css.
siril_resource_xml = $(srcdir)/siril_resource_gtk4.xml

siril_resource.c: siril_resource.h
	$(AM_V_GEN)
	glib-compile-resources \
	--sourcedir=$(top_builddir) \
	--sourcedir=$(top_srcdir) \
	--sourcedir=$(top_builddir)/src \
	--sourcedir=$(srcdir) \
	--target=$@ \
	--generate-source \
	--c-name siril_resource \
	$(siril_resource_xml)

# ── Core (GTK-free) sources ─────────────────────────────────────────────────
# Compiled into both siril and siril-cli.  Must not include <gtk/gtk.h> or
# any other GUI-only header (directly or transitively).  This mirrors
# meson's src_files list — keep them in sync.
siril_core_sources = \
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
	algos/demosaicing.c \
	algos/demosaicing.h \
	algos/demosaicing_rtp.cpp \
	algos/extraction.c \
	algos/extraction.h \
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
	core/cut.c \
	core/cut.h \
	core/gui_iface.h \
	core/gui_iface_stubs.c \
	core/icc_profile.c \
	core/icc_profile.h \
	core/initfile.c \
	core/initfile.h \
	core/masks.c \
	core/masks.h \
	core/OS_utils.c \
	core/OS_utils.h \
	core/pipe.c \
	core/pipe.h \
	core/preprocess.c \
	core/preprocess.h \
	core/processing_thread.c \
	core/processing_thread.h \
	core/processing.c \
	core/processing.h \
	core/sequence_filtering.c \
	core/sequence_filtering.h \
	core/settings.c \
	core/settings.h \
	core/signals.c \
	core/signals.h \
	core/siril.c \
	core/siril_alloc.c \
	core/siril_alloc.h \
	core/siril_app_dirs.c \
	core/siril_app_dirs.h \
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
	drizzle/cdrizzlemap.h \
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
	io/aavso_extended.c \
	io/aavso_extended.h \
	io/annotation_catalogues.c \
	io/annotation_catalogues.h \
	io/Astro-TIFF.c \
	io/Astro-TIFF.h \
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
	io/healpix/fluxcache.cpp \
	io/healpix/fluxcache.h \
	io/healpix/gaia_xp_design.c \
	io/healpix/gaia_xp_design.h \
	io/healpix/healpix.cpp \
	io/healpix/healpix_cat.h \
	io/healpix/healpix_image.cpp \
	io/healpix/healpix_image.h \
	io/healpix/xp_continuous.cpp \
	io/healpix/xp_continuous.h \
	io/image_format_fits.c \
	io/image_format_fits.h \
	io/image_formats_internal.c \
	io/image_formats_libraries.c \
	io/kstars/binfile.h \
	io/kstars/byteorder.h \
	io/kstars/htmesh_wrapper.cpp \
	io/kstars/htmesh_wrapper.h \
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
	io/SirilXISFReader.cpp \
	io/SirilXISFReader.h \
	io/SirilXISFWraper.h \
	io/SirilJpegXLWrapper.cpp \
	io/spcc_json.c \
	livestacking/livestacking.c \
	livestacking/livestacking.h \
	opencv/opencv.cpp \
	opencv/opencv.h \
	opencv/guidedfilter.cpp \
	opencv/guidedfilter.h \
	opencv/kombat/kombat.cpp \
	opencv/kombat/kombat.h \
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
	core/exif.cpp \
	core/exif.h \
	rt/rt_algo.cc \
	rt/gauss.cc \
	rt/boxblur.cc

# ── GUI (GTK-dependent) sources ─────────────────────────────────────────────
# Compiled into siril only.  Each file may include <gtk/gtk.h>.
# Mirrors meson's src_files_gui — keep in sync.  The GUI is GTK4-only and
# lives in src/gui-gtk4.
if USE_GTK4
siril_gui_sources = \
	gui-gtk4/about_dialog.c \
	gui-gtk4/about_dialog.h \
	gui-gtk4/annotate.c \
	gui-gtk4/annotations_pref.c \
	gui-gtk4/asinh.c \
	gui-gtk4/astrometry_solver.c \
	gui-gtk4/background_extraction.c \
	gui-gtk4/banding.c \
	gui-gtk4/callbacks.c \
	gui-gtk4/callbacks.h \
	gui-gtk4/catmag.c \
	gui-gtk4/catmag.h \
	gui-gtk4/ccd-inspector.c \
	gui-gtk4/clahe.c \
	gui-gtk4/colors.c \
	gui-gtk4/colors.h \
	gui-gtk4/compositing.c \
	gui-gtk4/compstars.c \
	gui-gtk4/compstars.h \
	gui-gtk4/conversion.c \
	gui-gtk4/conversion.h \
	gui-gtk4/cosmetic_correction.c \
	gui-gtk4/curves.c \
	gui-gtk4/curves.h \
	gui-gtk4/cut.c \
	gui-gtk4/cut.h \
	gui-gtk4/denoisegui.c \
	gui-gtk4/dialog_preview.c \
	gui-gtk4/dialog_preview.h \
	gui-gtk4/dialogs.c \
	gui-gtk4/dialogs.h \
	gui-gtk4/documentation.c \
	gui-gtk4/documentation.h \
	gui-gtk4/epf.c \
	gui-gtk4/epf.h \
	gui-gtk4/fft.c \
	gui-gtk4/file_browser.c \
	gui-gtk4/file_browser.h \
	gui-gtk4/fix_xtrans_af.c \
	gui-gtk4/fix_xtrans_af.h \
	gui-gtk4/geometry.c \
	gui-gtk4/git_gui.c \
	gui-gtk4/gui_iface_impl.c \
	gui-gtk4/gui_state.h \
	gui-gtk4/gtk3_event_compat.h \
	gui-gtk4/histo_display.c \
	gui-gtk4/histo_display.h \
	gui-gtk4/histogram.c \
	gui-gtk4/histogram.h \
	gui-gtk4/histogram_utils.c \
	gui-gtk4/histogram_utils.h \
	gui-gtk4/icc_profile.c \
	gui-gtk4/icc_profile.h \
	gui-gtk4/image_display.c \
	gui-gtk4/image_display.h \
	gui-gtk4/image_interactions.c \
	gui-gtk4/image_interactions.h \
	gui-gtk4/keywords_tree.c \
	gui-gtk4/keywords_tree.h \
	gui-gtk4/linear_match.c \
	gui-gtk4/livestacking_gui.c \
	gui-gtk4/masks_gui.c \
	gui-gtk4/masks_gui.h \
	gui-gtk4/median.c \
	gui-gtk4/median.h \
	gui-gtk4/menu_gray_geometry.c \
	gui-gtk4/menu_gray_geometry.h \
	gui-gtk4/menu_gray_selection.c \
	gui-gtk4/merge_cfa.c \
	gui-gtk4/message_dialog.c \
	gui-gtk4/message_dialog.h \
	gui-gtk4/mouse_action_functions.c \
	gui-gtk4/mouse_action_functions.h \
	gui-gtk4/mouse_action_preferences.c \
	gui-gtk4/newdeconv.c \
	gui-gtk4/nina_light_curve.c \
	gui-gtk4/nina_light_curve.h \
	gui-gtk4/open_dialog.c \
	gui-gtk4/open_dialog.h \
	gui-gtk4/photometric_cc.c \
	gui-gtk4/photometric_cc.h \
	gui-gtk4/pixelmath.c \
	gui-gtk4/plot.c \
	gui-gtk4/plot.h \
	gui-gtk4/preferences.c \
	gui-gtk4/preferences.h \
	gui-gtk4/preprocess.c \
	gui-gtk4/progress_and_log.c \
	gui-gtk4/progress_and_log.h \
	gui-gtk4/PSF_list.h \
	gui-gtk4/PSF_list.c \
	gui-gtk4/python_gui.h \
	gui-gtk4/python_gui.c \
	gui-gtk4/registration.h \
	gui-gtk4/registration.c \
	gui-gtk4/registration_preview.c \
	gui-gtk4/registration_preview.h \
	gui-gtk4/remixer.c \
	gui-gtk4/remixer.h \
	gui-gtk4/rgradient.c \
	gui-gtk4/saturation.c \
	gui-gtk4/saturation.h \
	gui-gtk4/save_dialog.c \
	gui-gtk4/save_dialog.h \
	gui-gtk4/scnr.c \
	gui-gtk4/script_console.c \
	gui-gtk4/script_console.h \
	gui-gtk4/script_menu.c \
	gui-gtk4/script_menu.h \
	gui-gtk4/search_objects.c \
	gui-gtk4/search_objects.h \
	gui-gtk4/sequence_export.c \
	gui-gtk4/sequence_list.c \
	gui-gtk4/sequence_list.h \
	gui-gtk4/single_image.c \
	gui-gtk4/single_image.h \
	gui-gtk4/siril_actions.c \
	gui-gtk4/siril_actions.h \
	gui-gtk4/siril_cmd_help.c \
	gui-gtk4/siril_cmd_help.h \
	gui-gtk4/siril_css.c \
	gui-gtk4/siril_css.h \
	gui-gtk4/siril_intro.c \
	gui-gtk4/siril_intro.h \
	gui-gtk4/siril_preview.c \
	gui-gtk4/siril_preview.h \
	gui-gtk4/siril_plot.c \
	gui-gtk4/siril_plot.h \
	gui-gtk4/siril-window.c \
	gui-gtk4/siril-window.h \
	gui-gtk4/splashscreen.c \
	gui-gtk4/splashscreen.h \
	gui-gtk4/split_cfa.c \
	gui-gtk4/stacking.c \
	gui-gtk4/stacking.h \
	gui-gtk4/star_finder.c \
	gui-gtk4/starnetgui.c \
	gui-gtk4/starnetgui.h \
	gui-gtk4/statistics.c \
	gui-gtk4/statistics_list.c \
	gui-gtk4/ui_files.h \
	gui-gtk4/undo_gui.c \
	gui-gtk4/undo_gui.h \
	gui-gtk4/unpurple.c \
	gui-gtk4/unpurple.h \
	gui-gtk4/user_polygons.c \
	gui-gtk4/user_polygons.h \
	gui-gtk4/utils.c \
	gui-gtk4/utils.h \
	gui-gtk4/wavelets.c \
	gui-gtk4/wavelets.h
endif

# ── Headless stubs ──────────────────────────────────────────────────────────
# Compiled into siril-cli only.  Provides no-op implementations of GUI
# functions still called directly from core code (matches meson's
# siril_headless_lib).
siril_headless_sources = \
	core/headless_stubs.c

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

# ── Common link inputs ──────────────────────────────────────────────────────
# GTK_LIBS, GTKSOURCEVIEW_LIBS and GTK_MAC_LIBS are intentionally NOT
# included here: they are linked only into siril (the GUI executable).
siril_common_LDADD = \
	$(top_builddir)/subprojects/kplot/libkplot.a \
	$(FFMPEG_LIBS) \
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
	$(GLIB_LIBS) \
	$(GIO_LIBS) \
	$(CAIRO_LIBS) \
	$(PANGO_LIBS) \
	$(PANGOCAIRO_LIBS) \
	$(GDK_PIXBUF_LIBS) \
	$(FFMS2_LIBS) \
	$(LIBCONFIG_LIBS) \
	$(LIBCURL) \
	$(WCSLIB_LIBS) \
	$(EXTRA_LIBS) \
	$(LIBATOMIC_LIBS) \
	$(libdbghelp) \
	$(librtprocesslibs) \
	$(LIBGIT_LIBS) \
	$(SQLITE_LIBS) \
	$(LIBJXL_LIBS) \
	$(LIBJXL_THREADS_LIBS) \
	$(YYJSON_LIBS) \
	$(OPENCV_LIBS) \
	$(EXIV2_LIBS) \
	$(libhtmeshlibs) \
	$(libhealpix_cxxlibs)

# ── Targets ─────────────────────────────────────────────────────────────────
# siril-cli is always built; siril is built only when --enable-gtk.

siril_cli_SOURCES = $(siril_core_sources) $(siril_headless_sources) main-cli.c
siril_cli_LDADD = $(siril_common_LDADD)

# Unconditionally defined so they appear in EXTRA_DIST regardless of
# whether siril is built this configuration.
siril_SOURCES = $(siril_core_sources) $(siril_gui_sources) siril_resource.c siril_resource.h main.c
siril_LDADD = $(siril_common_LDADD) $(GTK_LIBS) ${GTK_MAC_LIBS} $(GTKSOURCEVIEW_LIBS) $(X11_LIBS)

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
siril_cli_LINK = $(CXXLINK)
