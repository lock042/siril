################################################################################
# Meson siril file
################################################################################

siril_inc_dir = include_directories('.')

################################################################################
# Generate configuration files

configure_file(
  input : 'git-version.h.in',
  output : 'git-version.h',
  configuration : git_data
)

configure_file(output : 'config.h',
               configuration : conf_data)

################################################################################
# Generate Siril static library with all sources

src_files = [
  'algos/anscombe.c',
  'algos/img_t/better_than_std.hpp',
  'algos/img_t/fftw_allocator.hpp',
  'algos/img_t/image.hpp',
  'algos/img_t/image_boundary.hpp',
  'algos/img_t/image_dft.hpp',
  'algos/img_t/image_expr.hpp',
  'algos/img_t/image_ops.hpp',
  'algos/img_t/img_fits.hpp',
  'algos/img_t/vec2.hpp',
  'algos/astrometry_solver.c',
  'algos/background_extraction.c',
  'algos/ccd-inspector.c',
  'algos/colors.c',
  'algos/comparison_stars.c',
  'algos/demosaicing.c',
  'algos/demosaicing_rtp.cpp',
  'algos/extraction.c',
  'algos/fitting.c',
  'algos/fix_xtrans_af.c',
  'algos/geometry.c',
  'algos/io_wave.c',
  'algos/median_fast.c',
  'algos/noise.c',
  'algos/pave.c',
  'algos/photometry.c',
  'algos/photometric_cc.c',
  'algos/spcc.c',
  'algos/PSF.c',
  'algos/quality.c',
  'algos/quality_float.c',
  'algos/quantize.c',
  'algos/reconstr.c',
  'algos/siril_random.cpp',
  'algos/siril_wcs.c',
  'algos/sorting.c',
  'algos/star_finder.c',
  'algos/statistics.c',
  'algos/statistics_float.c',
  'algos/transform.c',
  'algos/wavelet_denoise.c',
  'algos/search_objects.c',
  'algos/lcms_acceleration/fast_16_tethra.c',
  'algos/lcms_acceleration/fast_8_curves.c',
  'algos/lcms_acceleration/fast_8_matsh.c',
  'algos/lcms_acceleration/fast_8_matsh_sse.c',
  'algos/lcms_acceleration/fast_8_tethra.c',
  'algos/lcms_acceleration/fast_float_15bits.c',
  'algos/lcms_acceleration/fast_float_15mats.c',
  'algos/lcms_acceleration/fast_float_cmyk.c',
  'algos/lcms_acceleration/fast_float_curves.c',
  'algos/lcms_acceleration/fast_float_lab.c',
  'algos/lcms_acceleration/fast_float_matsh.c',
  'algos/lcms_acceleration/fast_float_separate.c',
  'algos/lcms_acceleration/fast_float_sup.c',
  'algos/lcms_acceleration/fast_float_tethra.c',
  'algos/lcms_acceleration/threaded_core.c',
  'algos/lcms_acceleration/threaded_main.c',
  'algos/lcms_acceleration/threaded_scheduler.c',
  'algos/lcms_acceleration/threaded_split.c',

  'compositing/align_rgb.c',
  'compositing/filters.c',

  'core/arithm.c',
  'core/command.c',
  'core/command_line_processor.c',
  'core/exif.cpp',
  'core/icc_profile.c',
  'core/initfile.c',
  'core/masks.c',
  'core/OS_utils.c',
  'core/pipe.c',
  'core/preprocess.c',
  'core/processing.c',
  'core/processing_thread.c',
  'core/sequence_filtering.c',
  'core/settings.c',
  'core/signals.c',
  'core/siril.c',
  'core/siril_alloc.c',
  'core/siril_app_dirs.c',
  'core/siril_date.c',
  'core/siril_language.c',
  'core/gui_iface_stubs.c',
  'core/siril_log.c',
  'core/siril_networking.c',
  'core/siril_spawn.c',
  'core/siril_update.c',
  'core/siril_world_cs.c',
  'core/undo.c',
  'core/utils.c',
  'core/cut.c',

  'drizzle/cdrizzlebox.c',
  'drizzle/cdrizzlemap.c',
  'drizzle/cdrizzleutil.c',

  'filters/asinh.c',
  'filters/banding.c',
  'filters/clahe.c',
  'filters/cosmetic_correction.c',
  'filters/curve_transform.c',
  'filters/deconvolution.c',
  'filters/epf.c',
  'filters/fft.c',
  'filters/ght.c',
  'filters/linear_match.c',
  'filters/median.c',
  'filters/mtf.c',
  'filters/rgradient.c',
  'filters/saturation.c',
  'filters/scnr.c',
  'filters/starnet.c',
  'filters/synthstar.c',
  'filters/unpurple.c',
  'filters/wavelets.c',

  'filters/da3d/DA3D.cpp',
  'filters/da3d/WeightMap.cpp',
  'filters/nlbayes/call_nlbayes.cpp',
  'filters/nlbayes/LibImages.cpp',
  'filters/nlbayes/LibMatrix.cpp',
  'filters/nlbayes/NlBayes.cpp',
  'filters/nlbayes/Utilities.cpp',

  'filters/deconvolution/angleSet.cpp',
  'filters/deconvolution/deconvolution.h',
  'filters/deconvolution/deconvolve.cpp',
  'filters/deconvolution/deconvolve.hpp',
  'filters/deconvolution/edgetaper.hpp',
  'filters/deconvolution/estimate_kernel.cpp',
  'filters/deconvolution/fft.hpp',
  'filters/deconvolution/gf_estimate.cpp',
  'filters/deconvolution/labeling.hpp',
  'filters/deconvolution/optimization.hpp',
  'filters/deconvolution/options.hpp',
  'filters/deconvolution/utils.hpp',

  'io/aavso_extended.c',
  'io/annotation_catalogues.c',
  'io/Astro-TIFF.c',
  'io/avi_pipp/pipp_avi_write.cpp',
  'io/avi_pipp/pipp_avi_write.h',
  'io/avi_preview.c',
  'io/avi_preview.h',
  'io/avi_pipp/pipp_avi_write_dib.cpp',
  'io/avi_pipp/pipp_avi_write_dib.h',
  'io/avi_pipp/pipp_utf8.h',
  'io/avi_pipp/pipp_buffer.cpp',
  'io/avi_pipp/pipp_buffer.h',
  'io/avi_pipp/pipp_video_write.h',
  'io/avi_pipp/avi_writer.cpp',
  'io/avi_pipp/avi_writer.h',
  'io/kstars/htmesh_wrapper.cpp',
  'io/conversion.c',
  'io/films.c',
  'io/fits_keywords.c',
  'io/fits_sequence.c',
  'io/FITS_symlink.c',
  'io/healpix/healpix.cpp',
  'io/healpix/healpix_image.cpp',
  'io/healpix/fluxcache.cpp',
  'io/healpix/xp_continuous.cpp',
  'io/healpix/gaia_xp_design.c',
  'io/image_format_fits.c',
  'io/image_formats_internal.c',
  'io/image_formats_libraries.c',
  'io/local_catalogues.c',
  'io/mp4_output.c',
  'io/path_parse.c',
  'io/remote_catalogues.c',
  'io/seqfile.c',
  'io/sequence.c',
  'io/sequence_export.c',
  'io/seqwriter.c',
  'io/ser.c',
  'io/single_image.c',
  'io/siril_catalogues.c',
  'io/siril_git.c',
  'io/siril_plot.c',
  'io/siril_pythoncommands.c',
  'io/siril_pythonmodule.c',
  'io/SirilXISFReader.cpp',
  'io/SirilJpegXLWrapper.cpp',
  'io/spcc_json.c',

  'livestacking/livestacking.c',

  'opencv/opencv.cpp',
  'opencv/opencv.h',
  'opencv/guidedfilter.cpp',
  'opencv/guidedfilter.h',
  'opencv/kombat/kombat.cpp',
  'opencv/kombat/kombat.h',

  'pixelMath/pixel_math_runner.c',
  'pixelMath/tinyexpr.c',

  'registration/3stars.c',
  'registration/comet.c',
  'registration/global.c',
  'registration/matching/match.c',
  'registration/matching/atpmatch.c',
  'registration/matching/misc.c',
  'registration/matching/apply_match.c',
  'registration/astrometric.c',
  'registration/distorsion.c',
  'registration/registration.c',
  'registration/applyreg.c',
  'registration/shift_methods.c',

  'stacking/blending.c',
  'stacking/median_and_mean.c',
  'stacking/normalization.c',
  'stacking/rejection_float.c',
  'stacking/siril_fit_linear.c',
  'stacking/stacking.c',
  'stacking/stackminmax.c',
  'stacking/sum.c',
  'stacking/upscaling.c',

  'rt/boxblur.cc',
  'rt/gauss.cc',
  'rt/rt_algo.cc'
]

# GUI-only source files (require GTK; excluded when -Dgtk=false).
#
# gui_subdir is fixed to 'gui-gtk4' in the top-level meson.build; the GUI
# is GTK4-only.
gui_filenames = [
  'livestacking_gui.c',
  'about_dialog.c',
  'annotate.c',
  'asinh.c',
  'banding.c',
  'clahe.c',
  'cosmetic_correction.c',
  'fft.c',
  'rgradient.c',
  'saturation.c',
  'scnr.c',
  'annotations_pref.c',
  'astrometry_solver.c',
  'denoisegui.c',
  'background_extraction.c',
  'callbacks.c',
  'undo_gui.c',
  'ccd-inspector.c',
  'colors.c',
  'catmag.c',
  'compstars.c',
  'conversion.c',
  'curves.c',
  'cut.c',
  'dialog_preview.c',
  'dialogs.c',
  'documentation.c',
  'epf.c',
  'gui_iface_impl.c',
  'fix_xtrans_af.c',
  'geometry.c',
  'git_gui.c',
  'histo_display.c',
  'histogram.c',
  'histogram_utils.c',
  'keywords_tree.c',
  'icc_profile.c',
  'image_display.c',
  'image_interactions.c',
  'linear_match.c',
  'masks_gui.c',
  'median.c',
  'menu_gray_geometry.c',
  'menu_gray_selection.c',
  'merge_cfa.c',
  'message_dialog.c',
  'mouse_action_functions.c',
  'mouse_action_preferences.c',
  'newdeconv.c',
  'nina_light_curve.c',
  'open_dialog.c',
  'photometric_cc.c',
  'plot.c',
  'preferences.c',
  'progress_and_log.c',
  'PSF_list.c',
  'python_gui.c',
  'registration.c',
  'registration_preview.c',
  'remixer.c',
  'save_dialog.c',
  'script_console.c',
  'script_menu.c',
  'search_objects.c',
  'sequence_list.c',
  'single_image.c',
  'statistics.c',
  'siril_css.c',
  'siril_intro.c',
  'siril_preview.c',
  'siril_plot.c',
  'siril-window.c',
  'splashscreen.c',
  'split_cfa.c',
  'stacking.c',
  'star_finder.c',
  'starnetgui.c',
  'statistics_list.c',
  'compositing.c',
  'unpurple.c',
  'user_polygons.c',
  'pixelmath.c',
  'preprocess.c',
  'sequence_export.c',
  'siril_actions.c',
  'siril_cmd_help.c',
  'utils.c',
  'wavelets.c',
  'file_browser.c',  # custom GtkDirectoryList-based file picker with preview
]

src_files_gui = []
foreach f : gui_filenames
  src_files_gui += gui_subdir / f
endforeach

# enable_gtk / gui_subdir are defined in the top-level meson.build
#
# Build strategy:
#
#  siril_lib         — core (non-GUI) sources only: src_files.  Contains
#                      gui_iface_stubs.c which provides no-op vtable entries.
#                      Never contains gui_iface_impl.c or headless_stubs.c.
#
#  siril_headless_lib — headless_stubs.c alone; always built.  Provides no-op
#                      implementations of GUI functions still called directly
#                      from core code (pending full gui_iface routing).
#                      Linked on-demand (link_with, not link_whole) so its
#                      symbols are only pulled in when unresolved: in siril-cli
#                      and tests they satisfy the calls from core; in siril the
#                      same symbols are already defined by src_files_gui so the
#                      stubs are never included (no duplicate-symbol risk).
#
#  siril             — link_whole: siril_lib  +  src_files_gui compiled directly.
#                      gui_iface_impl.c (in src_files_gui) provides real GTK
#                      implementations; main.c defines guiinfo gui.  Only built
#                      when -Dgtk=true (the default).
#
#  siril-cli         — link_whole: siril_lib  +  siril_headless_lib (on demand).
#                      gui_iface_impl.c is NEVER linked here; gui_iface_stubs.c
#                      (in siril_lib) provides headless no-ops for all slots.
#                      main-cli.c does NOT define guiinfo gui.  Always built.

# Define the resource file.  It points at gui-gtk4/uifiles/* and
# gui-gtk4/siril.css.
gnome = import('gnome')
resource_xml = 'siril_resource_gtk4.xml'
siril_resources = gnome.compile_resources(
  'siril_resource', resource_xml,
  source_dir : [ '..', '.' ],
  c_name : 'siril_resource')

siril_lib = static_library('siril',
                           siril_resources[0], src_files,
                           dependencies: core_dep,
                           include_directories: siril_inc_dir,
                           c_args : siril_c_flag,
                           cpp_args : siril_cpp_flag)

siril_headless_lib = static_library('siril_headless',
                                    'core/headless_stubs.c',
                                    dependencies: core_dep,
                                    include_directories: siril_inc_dir,
                                    c_args : siril_c_flag,
                                    cpp_args : siril_cpp_flag)

siril_dep = declare_dependency(link_with: [siril_lib, siril_headless_lib],
                               dependencies: externals_dep,
                               include_directories: siril_inc_dir)

# Custom target for glib-compile-resources

################################################################################
# Build gui executable (only when -Dgtk=true, the default)

if enable_gtk
siril_gui = executable('siril', siril_resources[1], 'main.c', src_files_gui, win_res,
                       dependencies: siril_dep,
                       link_whole: siril_lib,
                       link_args : siril_link_arg,
                       c_args : siril_c_flag,
                       cpp_args : siril_cpp_flag,
                       install : true,
                       win_subsystem : 'console',
                      )
endif

################################################################################
# Build cli executable
# siril-cli links siril_lib (core + gui_iface_stubs) and siril_headless_lib
# (on-demand stubs for remaining direct GUI calls from core code).
# gui_iface_impl.c is never linked here.  main-cli.c does not define guiinfo gui.

siril_cli = executable('siril-cli', 'main-cli.c', win_res,
                       dependencies: siril_dep,
                       link_whole: siril_lib,
                       link_args : siril_link_arg,
                       c_args : siril_c_flag,
                       cpp_args : siril_cpp_flag,
                       install : true
                      )

################################################################################
# Tests

subdir('tests')
