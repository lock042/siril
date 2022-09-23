#ifndef SRC_CORE_COMMAND_LIST_H_
#define SRC_CORE_COMMAND_LIST_H_


#include "core/siril.h"
#include "core/command.h"
#include "core/command_def.h"

#define MAX_COMMAND_WORDS 50		// max number of words to split in command line input

extern char *word[MAX_COMMAND_WORDS];	// NULL terminated

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

static command commands[] = {
	/* name, nbarg, usage, function pointer, description, scriptable */
	{"addmax", 1, "addmax filename", process_addmax, STR_ADDMAX, FALSE, REQ_CMD_SINGLE_IMAGE},
	{"autostretch", 0, "autostretch [-linked] [shadowsclip [targetbg]]", process_autostretch, STR_AUTOSTRETCH, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"asinh", 1, "asinh [-human] stretch [offset]", process_asinh, STR_ASINH, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"bg", 0, "bg", process_bg, STR_BG, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"bgnoise", 0, "bgnoise", process_bgnoise, STR_BGNOISE, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"boxselect", 0, "boxselect [-clear] [x y width height]", process_boxselect, STR_BOXSELECT, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE | REQ_CMD_NO_THREAD},

	{"capabilities", 0, "capabilities", process_capabilities, STR_CAPABILITIES, TRUE, REQ_CMD_NONE},
	{"cd", 1, "cd directory", process_cd, STR_CD, TRUE, REQ_CMD_NONE},
	{"cdg", 0, "cdg", process_cdg, STR_CDG, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"clahe", 2, "clahe cliplimit tileSize", process_clahe, STR_CLAHE, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"clear", 0, "clear", process_clear, STR_CLEAR, FALSE, REQ_CMD_NONE},
	{"clearstar", 0, "clearstar", process_clearstar, STR_CLEARSTAR, FALSE, REQ_CMD_NONE},
	{"close", 0, "close", process_close, STR_CLOSE, TRUE, REQ_CMD_NONE},
	{"convert", 1, "convert basename [-debayer] [-fitseq] [-ser] [-start=index] [-out=]", process_convert, STR_CONVERT, TRUE, REQ_CMD_NO_THREAD},
	{"convertraw", 1, "convertraw basename [-debayer] [-fitseq] [-ser] [-start=index] [-out=]", process_convertraw, STR_CONVERTRAW, TRUE, REQ_CMD_NO_THREAD},
	{"cosme", 1, "cosme [filename].lst", process_cosme, STR_COSME, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"cosme_cfa", 1, "cosme_cfa [filename].lst", process_cosme, STR_COSME_CFA, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_CFA},
	{"crop", 0, "crop [x y width height]", process_crop, STR_CROP, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"ddp", 3, "ddp level coef sigma", process_ddp, STR_DDP, FALSE, REQ_CMD_SINGLE_IMAGE},
#ifdef _WIN32
	{"dir", 0, "dir", process_ls, STR_LS, FALSE, REQ_CMD_NONE},
#endif
	{"dumpheader", 0, "dumpheader", process_dumpheader, STR_DUMPHEADER, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},

	{"entropy", 0, "entropy", process_entropy, STR_ENTROPY, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"exit", 0, "exit", process_exit, STR_EXIT, TRUE, REQ_CMD_NONE},
	{"extract", 1, "extract NbPlans", process_extract, STR_EXTRACT, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"extract_Ha", 0, "extract_Ha", process_extractHa, STR_EXTRACTHA, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_CFA},
	{"extract_Green", 0, "extract_Green", process_extractGreen, STR_EXTRACTGREEN, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_CFA},
	{"extract_HaOIII", 0, "extract_HaOIII", process_extractHaOIII, STR_EXTRACTHAOIII, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_CFA},

	{"fdiv", 2, "fdiv filename scalar", process_fdiv, STR_FDIV, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"fftd", 2, "fftd modulus phase", process_fft, STR_FFTD, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"ffti", 2, "ffti modulus phase", process_fft, STR_FFTI, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"fill", 1, "fill value [x y width height]", process_fill, STR_FILL, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"fill2", 1, "fill2 value [x y width height]", process_fill2, STR_FILL2, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"find_cosme", 2, "find_cosme cold_sigma hot_sigma", process_findcosme, STR_FIND_COSME, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"find_cosme_cfa", 2, "find_cosme_cfa cold_sigma hot_sigma", process_findcosme, STR_FIND_COSME_CFA, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_CFA | REQ_CMD_NO_THREAD},
	{"find_hot", 3, "find_hot filename cold_sigma hot_sigma", process_findhot, STR_FIND_HOT, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"findstar", 0, "findstar [-out=] [-layer=] [-maxstars=]", process_findstar, STR_FINDSTAR, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"fix_xtrans", 0, "fix_xtrans", process_fix_xtrans, STR_FIXXTRANS, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE | REQ_CMD_FOR_MONO},
	{"fixbanding", 2, "fixbanding amount sigma", process_fixbanding, STR_FIXBANDING, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"fmedian", 2, "fmedian ksize modulation", process_fmedian, STR_FMEDIAN, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"fmul", 1, "fmul scalar", process_fmul, STR_FMUL, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"gauss", 1, "gauss sigma", process_gauss, STR_GAUSS, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"get", 1, "get { -a | -A | variable }", process_set, STR_GET, TRUE, REQ_CMD_NONE},
	{"ght", 5, "ght [-human | -even | -independent] D B LP SP HP [channels]", process_ght, STR_GHT, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"grey_flat", 0, "grey_flat", process_grey_flat, STR_GREY_FLAT, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"help", 0, "help [command]", process_help, STR_HELP, TRUE, REQ_CMD_NONE},
	{"histo", 1, "histo channel (channel=0, 1, 2 with 0: red, 1: green, 2: blue)", process_histo, STR_HISTO, TRUE, REQ_CMD_SINGLE_IMAGE},

	/* commands open filename and current image */
	{"iadd", 1, "iadd filename", process_imoper, STR_IADD, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"idiv", 1, "idiv filename", process_imoper, STR_IDIV, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"inspector", 0, "inspector", process_inspector, STR_INSPECTOR, FALSE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"invght", 5, "invght [-human | -even | -independent] D B LP SP HP [channels]", process_invght, STR_INVGHT, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"invmodasinh", 4, "invmodasinh [-human | -even | -independent] D LP SP HP [channels]", process_invmodasinh, STR_INVMODASINH, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"imul", 1, "imul filename", process_imoper, STR_IMUL, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"isub", 1, "isub filename", process_imoper, STR_ISUB, TRUE, REQ_CMD_SINGLE_IMAGE},


	{"light_curve", 3, "light_curve sequencename channel { -ninastars=file | { -at=x,y | -wcs=ra,dec } { -refat=x,y | -refwcs=ra,dec } ...}", process_light_curve, STR_LIGHTCURVE, TRUE, REQ_CMD_NO_THREAD},
	{"linear_match", 2, "linear_match reference low high", process_linear_match, STR_LMATCH, TRUE, REQ_CMD_SINGLE_IMAGE}, /* logarifies current image */
	{"link", 1, "link basename [-start=index] [-out=]", process_link, STR_LINK, TRUE, REQ_CMD_NO_THREAD},
	{"linstretch", 1, "linstretch BP [channels]", process_linstretch, STR_LINSTRETCH, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"load", 1, "load filename[.ext]", process_load, STR_LOAD, TRUE, REQ_CMD_NONE},
	// specific loads are not required, but could be used to force the
	// extension to a higher priority in case two files with same basename
	// exist (stat_file() manages that priority order for now).
	{"log", 0, "log", process_log, STR_LOG, TRUE, REQ_CMD_SINGLE_IMAGE}, /* logarifies current image */
#ifndef _WIN32
	{"ls", 0, "ls", process_ls, STR_LS, FALSE, REQ_CMD_NONE},
#endif

	{"merge", 3, "merge sequence1 sequence2 [sequence3 ...] output_sequence", process_merge, STR_MERGE, TRUE, REQ_CMD_NONE},
	{"mirrorx", 0, "mirrorx", process_mirrorx, STR_MIRRORX, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"mirrory", 0, "mirrory", process_mirrory, STR_MIRRORY, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"modasinh", 4, "modasinh [-human | -even | -independent] D LP SP HP [channels]", process_modasinh, STR_MODASINH, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"mtf", 3, "mtf low mid high [channels]", process_mtf, STR_MTF, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},

	{"neg", 0, "neg", process_neg, STR_NEG, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"new", 3, "new width height nb_channel", process_new, STR_NEW, FALSE, REQ_CMD_NONE},
	{"nomad", 0, "nomad [limit_magnitude]", process_nomad, STR_NOMAD, FALSE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"nozero", 1, "nozero level", process_nozero, STR_NOZERO, TRUE, REQ_CMD_SINGLE_IMAGE}, /* replaces null values by level */

	{"offset", 1, "offset value", process_offset, STR_OFFSET, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"pcc", 0, "pcc [image_center_coords] [-noflip] [-platesolve] [-focal=] [-pixelsize=]", process_pcc, STR_PCC, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"pm", 1, "pm \"expression\" [-rescale [low] [high]]", process_pm, STR_PM, TRUE, REQ_CMD_NONE},
	{"preprocess", 1, "preprocess sequencename [-bias=filename] [-dark=filename] [-flat=filename] [-cc=dark [siglo sighi] || -cc=bpm bpmfile] [-cfa] [-debayer] [-fix_xtrans] [-equalize_cfa] [-opt] [-prefix=] [-fitseq]", process_preprocess, STR_PREPROCESS, TRUE, REQ_CMD_NONE},
	{"preprocess_single", 1, "preprocess_single imagename [-bias=filename] [-dark=filename] [-flat=filename] [-cfa] [-debayer] [-fix_xtrans] [-equalize_cfa] [-opt] [-prefix=]", process_preprocess_single, STR_PREPROCESS_SINGLE, TRUE, REQ_CMD_NONE},
	{"psf", 0, "psf [channel]", process_psf, STR_PSF, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"register", 1, "register sequence [-2pass] [-noout] [-drizzle] [-prefix=] [-minpairs=] [-transf=] [-layer=] [-maxstars=] [-nostarlist] [-interp=] [-selected]", process_register, STR_REGISTER, TRUE, REQ_CMD_NO_THREAD},
	{"reloadscripts", 0, "reloadscripts", process_reloadscripts, STR_RELOADSCRIPTS, FALSE, REQ_CMD_NONE},
	{"requires", 1, "requires version", process_requires, STR_REQUIRES, TRUE, REQ_CMD_NONE},
	{"resample", 1, "resample factor", process_resample, STR_RESAMPLE, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"rgbcomp", 2, "rgbcomp [-lum=image [rgb_image]] [red green blue] [-out=result_filename]", process_rgbcomp, STR_RGBCOMP, TRUE, REQ_CMD_NONE},
	{"rgradient", 4, "rgradient xc yc dR dalpha", process_rgradient, STR_RGRADIENT, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"rl", 3, "rl sigma corner_radius_boost iterations", process_rl, STR_RL, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"rmgreen", 0, "rmgreen [type]", process_scnr, STR_RMGREEN, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_RGB | REQ_CMD_NO_THREAD},
	{"rotate", 1, "rotate degree [-nocrop]", process_rotate, STR_ROTATE, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"rotatePi", 0, "rotatePi", process_rotatepi, STR_ROTATEPI, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"satu", 1, "satu amount [background_factor [hue_range_index]]", process_satu, STR_SATU, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_RGB | REQ_CMD_NO_THREAD},
	{"save", 1, "save filename", process_save, STR_SAVE, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"savebmp", 1, "savebmp filename", process_savebmp, STR_SAVEBMP, TRUE, REQ_CMD_SINGLE_IMAGE},
#ifdef HAVE_LIBJPEG
	{"savejpg", 1, "savejpg filename [quality]", process_savejpg, STR_SAVEJPG, TRUE, REQ_CMD_SINGLE_IMAGE},
#endif
#ifdef HAVE_LIBPNG
	{"savepng", 1, "savepng filename", process_savepng, STR_SAVEPNG, TRUE, REQ_CMD_SINGLE_IMAGE},
#endif
	{"savepnm", 1, "savepnm filename", process_savepnm, STR_SAVEPNM, TRUE, REQ_CMD_SINGLE_IMAGE},
#ifdef HAVE_LIBTIFF
	{"savetif", 1, "savetif filename", process_savetif, STR_SAVETIF, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"savetif32", 1, "savetif32 filename", process_savetif, STR_SAVETIF32, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"savetif8", 1, "savetif8 filename", process_savetif, STR_SAVETIF8, TRUE, REQ_CMD_SINGLE_IMAGE},
#endif
	{"select", 2, "select from to", process_select, STR_SELECT, FALSE, REQ_CMD_SEQUENCE},
	{"seqapplyreg", 1, "seqapplyreg sequencename [-drizzle] [-interp=] [-layer=] [-framing=] [-prefix=] [-filter-fwhm=value[%]] [-filter-wfwhm=value[%]] [-filter-round=value[%]] [-filter-bkg=value[%]] [-filter-nbstars=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]]", process_seq_applyreg, STR_SEQAPPLYREG, TRUE, REQ_CMD_NO_THREAD},
	{"seqclean", 1, "seqclean sequencename [-reg] [-stat] [-sel]", process_seq_clean, STR_SEQCLEAN, TRUE, REQ_CMD_NONE},
	{"seqextract_Ha", 1, "seqextract_Ha sequencename [-prefix=]", process_seq_extractHa, STR_SEQEXTRACTHA, TRUE, REQ_CMD_NO_THREAD},
	{"seqextract_Green", 1, "seqextract_Green sequencename [-prefix=]", process_seq_extractGreen, STR_SEQEXTRACTGREEN, TRUE, REQ_CMD_NO_THREAD},
	{"seqextract_HaOIII", 1, "seqextract_HaOIII sequencename", process_seq_extractHaOIII, STR_SEQEXTRACTHAOIII, TRUE, REQ_CMD_NO_THREAD},
	{"seqcosme", 2, "seqcosme sequencename [filename].lst [-prefix=]", process_seq_cosme, STR_SEQCOSME, TRUE, REQ_CMD_NONE},
	{"seqcosme_cfa", 2, "seqcosme_cfa sequencename [filename].lst [-prefix=]", process_seq_cosme, STR_SEQCOSME_CFA, TRUE, REQ_CMD_NONE},
	{"seqcrop", 5, "seqcrop sequencename x y width height [-prefix=]", process_seq_crop, STR_SEQCROP, TRUE, REQ_CMD_NO_THREAD},
	{"seqheader", 2, "seqheader sequencename keyword", process_seq_header, STR_SEQHEADER, TRUE, REQ_CMD_NONE},
	{"seqfind_cosme", 3, "seqfind_cosme sequencename cold_sigma hot_sigma [-prefix=]", process_findcosme, STR_SEQFIND_COSME, TRUE, REQ_CMD_NONE},
	{"seqfind_cosme_cfa", 3, "seqfind_cosme_cfa sequencename cold_sigma hot_sigma [-prefix=]", process_findcosme, STR_SEQFIND_COSME_CFA, TRUE, REQ_CMD_NONE},
	{"seqfindstar", 1, "seqfindstar sequencename [-layer=] [-maxstars=]", process_seq_findstar, STR_FINDSTAR, TRUE, REQ_CMD_NONE},
	{"seqmtf", 4, "seqmtf sequencename low mid high [channels] [-prefix=]", process_seq_mtf, STR_SEQMTF, TRUE, REQ_CMD_NONE},
	{"seqpsf", 0, "seqpsf [sequencename channel { -at=x,y | -wcs=ra,dec }]", process_seq_psf, STR_SEQPSF, TRUE, REQ_CMD_NO_THREAD},
	{"seqsplit_cfa", 1, "seqsplit_cfa sequencename [-prefix=]", process_seq_split_cfa, STR_SEQSPLIT_CFA, TRUE, REQ_CMD_NO_THREAD},
	{"seqstat", 2, "seqstat sequencename output [option]", process_seq_stat, STR_SEQSTAT, TRUE, REQ_CMD_NO_THREAD},
	{"seqsubsky", 2, "seqsubsky sequencename { -rbf | degree } [-samples=20] [-tolerance=1.0] [-smooth=0.5] [-prefix=]", process_subsky, STR_SEQSUBSKY, TRUE, REQ_CMD_NONE},
	{"seqtilt", 0, "seqtilt [sequencename]", process_seq_tilt, STR_SEQTILT, TRUE, REQ_CMD_NO_THREAD},
	{"set", 1, "set { -import=inifilepath | variable=value }", process_set, STR_SET, TRUE, REQ_CMD_NONE},
	{"set16bits", 0, "set16bits", process_set_32bits, STR_SET16, TRUE, REQ_CMD_NONE},
	{"set32bits", 0, "set32bits", process_set_32bits, STR_SET32, TRUE, REQ_CMD_NONE},
	{"setcompress", 1, "setcompress 0/1 [-type=] [q]", process_set_compress, STR_SETCOMPRESS, TRUE, REQ_CMD_NONE},
#ifdef _OPENMP
	{"setcpu", 1, "setcpu number", process_set_cpu, STR_SETCPU, TRUE, REQ_CMD_NONE},
#endif
	{"setext", 1, "setext extension", process_set_ext, STR_SETEXT, TRUE, REQ_CMD_NONE},
	{"setfindstar", 0, "setfindstar [reset] [-radius=] [-sigma=] [-roundness=] [-focal=] [-pixelsize=] [-auto=on|off] [-convergence=] [-relax=on|off]" , process_set_findstar, STR_SETFINDSTAR, TRUE, REQ_CMD_NONE},
	{"setmag", 1, "setmag magnitude", process_set_mag, STR_SETMAG, FALSE, REQ_CMD_SINGLE_IMAGE},
	{"setmagseq", 1, "setmagseq magnitude", process_set_mag_seq, STR_SETMAGSEQ, FALSE, REQ_CMD_SEQUENCE},
	{"setmem", 1, "setmem ratio", process_set_mem, STR_SETMEM, TRUE, REQ_CMD_NONE},
	{"setphot", 0, "setphot [-inner=20] [-outer=30] [-aperture=10] [-force_radius=no] [-gain=2.3] [-min_val=0] [-max_val=60000]", process_set_photometry, STR_SETPHOT, TRUE, REQ_CMD_NONE},
	{"setref", 2, "setref sequencename image_number", process_set_ref, STR_SETREF, TRUE, REQ_CMD_NONE},
	{"split", 3, "split fileR fileG fileB", process_split, STR_SPLIT, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_RGB | REQ_CMD_NO_THREAD},
	{"split_cfa", 0, "split_cfa", process_split_cfa, STR_SPLIT_CFA, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_CFA},
	{"stack", 1, "stack sequencename [type] [rejection type] [sigma_low sigma_high] [-nonorm, norm=] [-output_norm] [-rgb_equal] [-out=result_filename] [-filter-fwhm=value[%]] [-filter-wfwhm=value[%]] [-filter-round=value[%]] [-filter-bkg=value[%]] [-filter-nbstars=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]] [-weight_from_noise] [-weight_from_nbstack] [-weight_from_nbstars] [-weight_from_wfwhm] [-fastnorm]", process_stackone, STR_STACK, TRUE, REQ_CMD_NONE},
	{"stackall", 0, "stackall [type] [rejection type] [sigma_low sigma_high] [-nonorm, norm=] [-output_norm] [-filter-fwhm=value[%]] [-filter-wfwhm=value[%]] [-filter-round=value[%]] [-filter-bkg=value[%]] [-filter-nbstars=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]] [-weight_from_noise] [-weight_from_nbstack] [-fastnorm]", process_stackall, STR_STACKALL, TRUE, REQ_CMD_NONE},
#ifdef HAVE_LIBTIFF
	{"starnet", 0, "starnet [-stretch] [-upscale] [-stride=value] [-nostarmask]", process_starnet, STR_STARNET, TRUE, REQ_CMD_SINGLE_IMAGE},
#endif
	{"stat", 0, "stat", process_stat, STR_STAT, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"subsky", 1, "subsky { -rbf | degree } [-samples=20] [-tolerance=1.0] [-smooth=0.5]", process_subsky, STR_SUBSKY, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"synthstar", 0, "synthstar", process_synthstar, STR_SYNTHSTAR, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"threshlo", 1, "threshlo level", process_threshlo, STR_THRESHLO, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"threshhi", 1, "threshi level", process_threshhi, STR_THRESHHI, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"thresh", 2, "thresh lo hi", process_thresh, STR_THRESH, TRUE, REQ_CMD_SINGLE_IMAGE}, /* threshes hi and lo */
	{"tilt", 0, "tilt [clear]", process_tilt, STR_TILT, FALSE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},

	{"unselect", 2, "unselect from to", process_unselect, STR_UNSELECT, FALSE, REQ_CMD_SEQUENCE},
	{"unsetmag", 0, "unsetmag", process_unset_mag, STR_UNSETMAG, FALSE, REQ_CMD_NONE},
	{"unsetmagseq", 0, "unsetmagseq", process_unset_mag_seq, STR_UNSETMAGSEQ, FALSE, REQ_CMD_SINGLE_IMAGE},
	{"unsharp", 2, "unsharp sigma multi", process_unsharp, STR_UNSHARP, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"visu", 2, "visu low high", process_visu, STR_VISU, FALSE, REQ_CMD_SINGLE_IMAGE},

	/* wavelet transform in nbr_plan plans */
	{"wavelet", 1, "wavelet nbr_plan type", process_wavelet, STR_WAVELET, TRUE, REQ_CMD_SINGLE_IMAGE},
	/* reconstruct from wavelet transform and weighs plans with c1, c2, c3... */
	{"wrecons", 2, "wrecons c1 c2 c3 ...", process_wrecons, STR_WRECONS, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"start_ls", 0, "start_ls [-dark=filename] [-flat=filename]", process_start_ls, STR_START_LS, TRUE, REQ_CMD_NO_THREAD},
	{"livestack", 1, "livestack filename", process_livestack, STR_LIVESTACK, TRUE, REQ_CMD_NONE},
	{"stop_ls", 0, "stop_ls", process_stop_ls, STR_STOP_LS, TRUE, REQ_CMD_NONE},

	{"",0,"",0, STR_NONE, FALSE, REQ_CMD_NONE}
};

#endif /* SRC_CORE_COMMAND_LIST_H_ */
