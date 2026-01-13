#ifndef SRC_CORE_COMMAND_H_
#define SRC_CORE_COMMAND_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#define MAX_COMMAND_WORDS 50		// max number of words to split in command line input

typedef enum {
	REQ_CMD_NONE = (1 << 0),
	REQ_CMD_SINGLE_IMAGE = (1 << 1),
	REQ_CMD_SEQUENCE = (1 << 2),
	REQ_CMD_FOR_CFA = (1 << 3),
	REQ_CMD_FOR_MONO = (1 << 4),
	REQ_CMD_FOR_RGB = (1 << 5),
	REQ_CMD_SELECTION = (1 << 6),
	REQ_CMD_NO_THREAD = (1 << 10)
} cmd_prerequires;

typedef
struct {
	char *name;
	int nbarg;
	char *usage;
	int (* process)(int);
	char *definition;
	gboolean scriptable;
	cmd_prerequires prerequires;
} command;

extern char *word[MAX_COMMAND_WORDS];  // NULL terminated

gboolean image_cfa_warning_check();
gboolean get_followstar_idle(gpointer user_data);

int	process_addmax(int nb);
int	process_autostretch(int nb);
int	process_autostretch_mask(int nb);
int	process_autoghs(int nb);
int	process_asinh(int nb);

int	process_bg(int nb);
int	process_bgnoise(int nb);
int	process_binarize_mask(int nb);
int	process_binxy(int nb);
int	process_blur_mask(int nb);
int	process_denoise(int nb);
gpointer run_nlbayes_on_fit(gpointer p);
gpointer run_bm3d_on_fit(gpointer p);

int	process_boxselect(int nb);

int	process_calibrate(int nb);
int	process_calibrate_single(int nb);
int	process_capabilities(int nb);
int	process_catsearch(int nb);
int	process_ccm(int nb);
int	process_cd(int nb);
int	process_cdg(int nb);
int	process_crop(int nb);
int	process_clahe(int nb);
int	process_clear(int nb);
int	process_clear_mask(int nb);
int	process_clearstar(int nb);
int	process_close(int nb);
int	process_conesearch(int nb);
int	process_convert(int nb);
int	process_cosme(int nb);

int	process_ddp(int nb);
int	process_disto(int nb);
int	process_dumpheader(int nb);

int	process_entropy(int nb);
int	process_epf(int nb);
int	process_exit(int nb);
int	process_extract(int nb);
int	process_extractGreen(int nb);
int	extract_Ha(extraction_scaling scaling);
int	process_extractHa(int nb);
int	extract_HaOIII(extraction_scaling scaling);
int	process_extractHaOIII(int nb);

int	process_fdiv(int nb);
int	process_feather_mask(int nb);
int	process_fft(int nb);
int	process_fill(int nb);
int	process_ffill(int nb);
int	process_findcompstars(int nb);
int	process_findcosme(int nb);
int	process_findhot(int nb);
int	process_findstar(int nb);
int	process_fix_xtrans(int nb);
int	process_fixbanding(int nb);
int	process_fmedian(int nb);
int	process_fmul(int nb);

int	process_gauss(int nb);
int	process_getref(int nb);
int	process_grey_flat(int nb);

int	process_help(int nb);
int	process_histo(int nb);

int	process_jsonmetadata(int nb);

int	process_icc_assign(int nb);
int	process_icc_convert_to(int nb);
int	process_icc_remove(int nb);
int	process_imoper(int nb);
int	process_invert_mask(int nb);
int	process_inspector(int nb);

int	process_light_curve(int nb);
int	process_limit(int nb);
int	process_link(int nb);
int	process_linear_match(int nb);
int	process_load(int nb);
int	process_load_seq(int nb);
int	process_log(int nb);
int	process_ls(int nb);

int	process_makepsf(int nb);
int process_mask_bitpix(int nb);
int	process_mask_fmul(int nb);
int	process_mask_from_channel(int nb);
int process_mask_from_color(int nb);
int	process_mask_from_lum(int nb);
int	process_mask_from_stars(int nb);
int	process_merge(int nb);
int	process_mirrorx(int nb);
int	process_mirrorx_single(int nb);
int	process_mirrory(int nb);
int	process_mtf(int nb);

int	process_neg(int nb);
int	process_new(int nb);
int	process_nozero(int nb);

int	process_offset(int nb);

int	process_ght(int nb);
int	process_invght(int nb);
int	process_invmodasinh(int nb);
int	process_linstretch(int nb);
int	process_modasinh(int nb);

int	process_offline(int nb);
int	process_online(int nb);

int	process_parse(int nb);
int	process_pcc(int nb);
int	process_platesolve(int nb);
int	process_pm(int nb);
int	process_profile(int nb);
int	process_psf(int nb);
int	process_pwd(int nb);
int	process_pyscript(int nb);

int	process_rebayer(int nb);
int	process_register(int nb);
int	process_resample(int nb);
int	process_reloadscripts(int nb);
int	process_requires(int nb);
int	process_rgbcomp(int nb);
int	process_rgradient(int nb);
int	process_rl(int nb);
int	process_rotate(int nb);
int	process_rotatepi(int nb);

int	process_satu(int nb);
int	process_save(int nb);
#ifdef HAVE_LIBHEIF
int	process_saveavif(int nb);
#endif
int	process_savebmp(int nb);
#ifdef HAVE_LIBHEIF
int	process_saveheif(int nb);
#endif
#ifdef HAVE_LIBJPEG
int	process_savejpg(int nb);
#endif
#ifdef HAVE_LIBJXL
int	process_savejxl(int nb);
#endif
#ifdef HAVE_LIBPNG
int	process_savepng(int nb);
#endif
int	process_savepnm(int nb);
#ifdef HAVE_LIBTIFF
int	process_savetif(int nb);
int	process_starnet(int nb);
#endif
int	process_sb(int nb);
int	process_scnr(int nb);
int	process_search_fct(int nb);
int	process_select(int nb);
int	process_seq_applyreg(int nb);
int	process_seq_clean(int nb);
int	process_seq_cosme(int nb);
int	process_seq_crop(int nb);
int	process_seq_extractHa(int nb);
int	process_seq_extractGreen(int nb);
int	process_seq_extractHaOIII(int nb);
int	process_seq_findstar(int nb);
int	process_seq_fixbanding(int nb);
int	process_seq_ght(int nb);
int	process_seq_header(int nb);
int	process_seq_invght(int nb);
int	process_seq_invmodasinh(int nb);
int	process_seq_linstretch(int nb);
int	process_seq_merge_cfa(int nb);
int	process_seq_modasinh(int nb);
int	process_seq_mtf(int nb);
int	process_seq_profile(int nb);
int	process_seq_psf(int nb);
int	process_seq_resample(int nb);
int	process_seq_rl(int nb);
int	process_seq_sb(int nb);
int	process_seq_split_cfa(int nb);
int	process_seq_starnet(int nb);
int	process_seq_stat(int nb);
int	process_seq_tilt(int nb);
int	process_seq_update_key(int nb);
int	process_seq_wiener(int nb);
int	process_set(int nb);
int	process_set_32bits(int nb);
int	process_set_compress(int nb);
#ifdef _OPENMP
int	process_set_cpu(int nb);
#endif
int	process_set_ext(int nb);
int	process_set_findstar(int nb);
int	process_set_stfbits(int nb);
int	process_set_mag(int nb);
int	process_set_mag_seq(int nb);
int	process_set_mem(int nb);
int	process_set_photometry(int nb);
int	process_set_ref(int nb);
int	process_subsky(int nb);
int	process_spcc(int nb);
int	process_spcc_list(int nb);
int	process_split(int nb);
int	process_split_cfa(int nb);
int	process_stat(int nb);
int	process_stackall(int nb);
int	process_stackone(int nb);
int	process_synthstar(int nb);

int	process_thresh(int nb);
int	process_threshlo(int nb);
int	process_threshhi(int nb);
int	process_tilt(int nb);
int	process_trixel(int nb);

int	process_unclip(int nb);
int	process_unpurple(int nb);
int	process_unset_mag(int nb);
int	process_unset_mag_seq(int nb);
int	process_unselect(int nb);
int	process_unsharp(int nb);
int process_update_key(int nb);

int	process_visu(int nb);

int	process_wavelet(int nb);
int	process_wiener(int nb);
int	process_wrecons(int nb);

/* live stacking specials */
int process_start_ls(int nb);
int process_livestack(int nb);
int process_stop_ls(int nb);

int	process_show(int nb);

#endif
