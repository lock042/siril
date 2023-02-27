#ifndef SRC_CORE_COMMAND_LIST_H_
#define SRC_CORE_COMMAND_LIST_H_


#include "core/siril.h"
#include "core/command.h"
#include "core/command_def.h"
#include "core/command_extra.h"

static command commands[] = {
	/* name, nbarg, usage, function pointer, description, scriptable */
	{"addmax", 1, "addmax filename", process_addmax, STR_ADDMAX, FALSE, REQ_CMD_SINGLE_IMAGE},
	{"asinh", 1, "asinh [-human] stretch [offset]", process_asinh, STR_ASINH, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"autoghs", 1, "autoghs [-linked] shadowsclip stretchamount [-b=] [-hp=] [-lp=]", process_autoghs, STR_AUTOGHS, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"autostretch", 0, "autostretch [-linked] [shadowsclip [targetbg]]", process_autostretch, STR_AUTOSTRETCH, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},

	{"bg", 0, "bg", process_bg, STR_BG, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"bgnoise", 0, "bgnoise", process_bgnoise, STR_BGNOISE, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"binxy", 1, "binxy coefficient [-sum]", process_binxy, STR_BINXY, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"boxselect", 0, "boxselect [-clear] [x y width height]", process_boxselect, STR_BOXSELECT, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE | REQ_CMD_NO_THREAD},

	{"calibrate", 1, "calibrate sequencename [-bias=filename] [-dark=filename] [-flat=filename] [-cc=dark [siglo sighi] || -cc=bpm bpmfile] [-cfa] [-debayer] [-fix_xtrans] [-equalize_cfa] [-opt] [-all] [-prefix=] [-fitseq]", process_calibrate, STR_CALIBRATE, TRUE, REQ_CMD_NONE},
	{"calibrate_single", 1, "calibrate_single imagename [-bias=filename] [-dark=filename] [-flat=filename] [-cfa] [-debayer] [-fix_xtrans] [-equalize_cfa] [-opt] [-prefix=]", process_calibrate_single, STR_CALIBRATE_SINGLE, TRUE, REQ_CMD_NONE},
	{"capabilities", 0, "capabilities", process_capabilities, STR_CAPABILITIES, TRUE, REQ_CMD_NONE},
	{"catsearch", 1, "catsearch star_name", process_catsearch, STR_CATSEARCH, TRUE, REQ_CMD_NONE},
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
	{"denoise", 0, "denoise [-nocosmetic] [-mod=m] [ -vst | -da3d | -sos=n [-rho=r] ] [-indep]", process_denoise, STR_DENOISE, TRUE, REQ_CMD_SINGLE_IMAGE},
#ifdef _WIN32
	{"dir", 0, "dir", process_ls, STR_LS, FALSE, REQ_CMD_NONE},
#endif
	{"dumpheader", 0, "dumpheader", process_dumpheader, STR_DUMPHEADER, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},

	{"entropy", 0, "entropy", process_entropy, STR_ENTROPY, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"exit", 0, "exit", process_exit, STR_EXIT, TRUE, REQ_CMD_NONE},
	{"extract", 1, "extract NbPlans", process_extract, STR_EXTRACT, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"extract_Green", 0, "extract_Green", process_extractGreen, STR_EXTRACTGREEN, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_CFA},
	{"extract_Ha", 0, "extract_Ha", process_extractHa, STR_EXTRACTHA, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_CFA},
	{"extract_HaOIII", 0, "extract_HaOIII [-resample=]", process_extractHaOIII, STR_EXTRACTHAOIII, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_CFA},

	{"fdiv", 2, "fdiv filename scalar", process_fdiv, STR_FDIV, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"ffill", 1, "ffill value [x y width height]", process_ffill, STR_FFILL, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"fftd", 2, "fftd modulus phase", process_fft, STR_FFTD, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"ffti", 2, "ffti modulus phase", process_fft, STR_FFTI, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"fill", 1, "fill value [x y width height]", process_fill, STR_FILL, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"find_cosme", 2, "find_cosme cold_sigma hot_sigma", process_findcosme, STR_FIND_COSME, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"find_cosme_cfa", 2, "find_cosme_cfa cold_sigma hot_sigma", process_findcosme, STR_FIND_COSME_CFA, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_CFA | REQ_CMD_NO_THREAD},
	{"find_hot", 3, "find_hot filename cold_sigma hot_sigma", process_findhot, STR_FIND_HOT, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"findstar", 0, "findstar [-out=] [-layer=] [-maxstars=]", process_findstar, STR_FINDSTAR, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"fix_xtrans", 0, "fix_xtrans", process_fix_xtrans, STR_FIXXTRANS, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE | REQ_CMD_FOR_MONO},
	{"fixbanding", 2, "fixbanding amount sigma [-vertical]", process_fixbanding, STR_FIXBANDING, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"fmedian", 2, "fmedian ksize modulation", process_fmedian, STR_FMEDIAN, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"fmul", 1, "fmul scalar", process_fmul, STR_FMUL, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"gauss", 1, "gauss sigma", process_gauss, STR_GAUSS, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"get", 1, "get { -a | -A | variable }", process_set, STR_GET, TRUE, REQ_CMD_NONE},
	{"getref", 1, "getref sequencename", process_getref, STR_GETREF, TRUE, REQ_CMD_NONE},
	{"ght", 1, "ght -D= [-B=] [-LP=] [-SP=] [-HP=] [-human | -even | -independent] [channels]", process_ght, STR_GHT, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"grey_flat", 0, "grey_flat", process_grey_flat, STR_GREY_FLAT, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"help", 0, "help [command]", process_help, STR_HELP, TRUE, REQ_CMD_NONE},
	{"histo", 1, "histo channel (channel=0, 1, 2 with 0: red, 1: green, 2: blue)", process_histo, STR_HISTO, TRUE, REQ_CMD_SINGLE_IMAGE},

	/* commands open filename and current image */
	{"iadd", 1, "iadd filename", process_imoper, STR_IADD, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"idiv", 1, "idiv filename", process_imoper, STR_IDIV, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"imul", 1, "imul filename", process_imoper, STR_IMUL, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"inspector", 0, "inspector", process_inspector, STR_INSPECTOR, FALSE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"invght", 1, "invght -D= [-B=] [-LP=] [-SP=] [-HP=] [-human | -even | -independent] [channels]", process_invght, STR_INVGHT, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"invmodasinh", 1, "-D= [-LP=] [-SP=] [-HP=] [-human | -even | -independent] [channels]", process_invmodasinh, STR_INVMODASINH, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"invmtf", 3, "invmtf low mid high [channels]", process_mtf, STR_INVMTF, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"isub", 1, "isub filename", process_imoper, STR_ISUB, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"jsonmetadata", 1, "jsonmetadata FITS_file [-stats_from_loaded] [-nostats] [-out=]", process_jsonmetadata, STR_JSONMETADATA, TRUE, REQ_CMD_NONE},

	{"light_curve", 3, "light_curve sequencename channel [-autoring] { -at=x,y | -wcs=ra,dec } { -refat=x,y | -refwcs=ra,dec } ...\n"
				"light_curve sequencename channel [-autoring] -ninastars=file", process_light_curve, STR_LIGHTCURVE, TRUE, REQ_CMD_NO_THREAD},
	{"linear_match", 2, "linear_match reference low high", process_linear_match, STR_LMATCH, TRUE, REQ_CMD_SINGLE_IMAGE}, /* logarifies current image */
	{"link", 1, "link basename [-start=index] [-out=]", process_link, STR_LINK, TRUE, REQ_CMD_NO_THREAD},
	{"linstretch", 1, "linstretch -BP= [channels]", process_linstretch, STR_LINSTRETCH, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"livestack", 1, "livestack filename", process_livestack, STR_LIVESTACK, TRUE, REQ_CMD_NONE},
	{"load", 1, "load filename[.ext]", process_load, STR_LOAD, TRUE, REQ_CMD_NONE},
	// specific loads are not required, but could be used to force the
	// extension to a higher priority in case two files with same basename
	// exist (stat_file() manages that priority order for now).
	{"log", 0, "log", process_log, STR_LOG, TRUE, REQ_CMD_SINGLE_IMAGE}, /* logarifies current image */
#ifndef _WIN32
	{"ls", 0, "ls", process_ls, STR_LS, FALSE, REQ_CMD_NONE},
#endif
	{"makepsf", 1, "makepsf clear [-ks=]\n"
				"makepsf load filename\n"
				"makepsf save\n"
				"makepsf blind [-l0] [-si] [-multiscale] [-lambda=] [-comp=] [-ks=]\n"
				"makepsf stars [-sym] [-ks=]\n"
				"makepsf manual { -gaussian | -moffat | -disc | -airy } [-fwhm=] [-angle=] [-ratio=] [-beta=] [-dia=] [-fl=] [-wl=] [-pixelsize=] [-obstruct=] [-ks=]", process_makepsf, STR_MAKEPSF, TRUE, REQ_CMD_NONE},
	{"merge", 3, "merge sequence1 sequence2 [sequence3 ...] output_sequence", process_merge, STR_MERGE, TRUE, REQ_CMD_NONE},
	{"merge_cfa", 5, "merge_cfa file_CFA0 file_CFA1 file_CFA2 file_CFA3 bayerpattern", process_rebayer, STR_REBAYER, TRUE, REQ_CMD_NONE},
	{"mirrorx", 0, "mirrorx [-bottomup]", process_mirrorx, STR_MIRRORX, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"mirrorx_single", 1, "mirrorx_single image", process_mirrorx_single, STR_MIRRORX_SINGLE, TRUE, REQ_CMD_NONE},
	{"mirrory", 0, "mirrory", process_mirrory, STR_MIRRORY, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"modasinh", 1, "modasinh -D= [-LP=] [-SP=] [-HP=] [-human | -even | -independent] [channels]", process_modasinh, STR_MODASINH, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"mtf", 3, "mtf low mid high [channels]", process_mtf, STR_MTF, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},

	{"neg", 0, "neg", process_neg, STR_NEG, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"new", 3, "new width height nb_channel", process_new, STR_NEW, FALSE, REQ_CMD_NONE},
	{"nomad", 0, "nomad [limit_magnitude] [-catalog=] [-photo]", process_nomad, STR_NOMAD, FALSE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"nozero", 1, "nozero level", process_nozero, STR_NOZERO, TRUE, REQ_CMD_SINGLE_IMAGE}, /* replaces null values by level */

	{"offset", 1, "offset value", process_offset, STR_OFFSET, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"parse", 1, "parse str [-r]", process_parse, STR_PARSE, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"pcc", 0, "pcc [image_center_coords] [-noflip] [-platesolve] [-focal=] [-pixelsize=] [-limitmag=[+-]] [-catalog=] [-downscale]", process_pcc, STR_PCC, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_RGB },
	{"platesolve", 0, "platesolve [image_center_coords] [-noflip] [-platesolve] [-focal=] [-pixelsize=] [-limitmag=[+-]] [-catalog=] [-localasnet] [-downscale]", process_pcc, STR_PLATESOLVE, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"pm", 1, "pm \"expression\" [-rescale [low] [high]]", process_pm, STR_PM, TRUE, REQ_CMD_NONE},
	{"preprocess", 1, "preprocess sequencename [-bias=filename] [-dark=filename] [-flat=filename] [-cc=dark [siglo sighi] || -cc=bpm bpmfile] [-cfa] [-debayer] [-fix_xtrans] [-equalize_cfa] [-opt] [-all] [-prefix=] [-fitseq]", process_preprocess, STR_CALIBRATE, TRUE, REQ_CMD_NONE},
	{"preprocess_single", 1, "preprocess_single imagename [-bias=filename] [-dark=filename] [-flat=filename] [-cfa] [-debayer] [-fix_xtrans] [-equalize_cfa] [-opt] [-prefix=]", process_preprocess_single, STR_CALIBRATE_SINGLE, TRUE, REQ_CMD_NONE},
	{"psf", 0, "psf [channel]", process_psf, STR_PSF, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"register", 1, "register sequencename [-2pass] [-noout] [-drizzle] [-prefix=] [-minpairs=] [-transf=] [-layer=] [-maxstars=] [-nostarlist] [-interp=] [-noclamp] [-selected]", process_register, STR_REGISTER, TRUE, REQ_CMD_NO_THREAD},
	{"reloadscripts", 0, "reloadscripts", process_reloadscripts, STR_RELOADSCRIPTS, FALSE, REQ_CMD_NONE},
	{"requires", 1, "requires version", process_requires, STR_REQUIRES, TRUE, REQ_CMD_NONE},
	{"resample", 1, "resample { factor | -width= | -height= } [-interp=] [-noclamp]", process_resample, STR_RESAMPLE, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"rgbcomp", 2, "rgbcomp red green blue [-out=result_filename]\n"
				"rgbcomp -lum=image { rgb_image | red green blue } [-out=result_filename]", process_rgbcomp, STR_RGBCOMP, TRUE, REQ_CMD_NONE},
	{"rgradient", 4, "rgradient xc yc dR dalpha", process_rgradient, STR_RGRADIENT, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"rl", 0, "rl [-alpha=] [-iters=] [-stop=] [-gdstep=] [-tv] [-fh] [-mul]", process_rl, STR_RL, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"rmgreen", 0, "rmgreen [-nopreserve] [type] [amount]", process_scnr, STR_RMGREEN, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_RGB | REQ_CMD_NO_THREAD},
	{"rotate", 1, "rotate degree [-nocrop] [-interp=] [-noclamp]", process_rotate, STR_ROTATE, TRUE, REQ_CMD_SINGLE_IMAGE},
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
	{"savetif", 1, "savetif filename [-astro]", process_savetif, STR_SAVETIF, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"savetif32", 1, "savetif32 filename [-astro]", process_savetif, STR_SAVETIF32, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"savetif8", 1, "savetif8 filename [-astro]", process_savetif, STR_SAVETIF8, TRUE, REQ_CMD_SINGLE_IMAGE},
#endif
	{"sb", 0, "sb [-alpha=] [-iters=]", process_sb, STR_SB, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"select", 3, "select sequencename from to", process_select, STR_SELECT, TRUE, REQ_CMD_NONE},
	{"seqapplyreg", 1, "seqapplyreg sequencename [-drizzle] [-interp=] [-noclamp] [-layer=] [-framing=] [-prefix=] [-filter-fwhm=value[%|k]] [-filter-wfwhm=value[%|k]] [-filter-round=value[%|k]] [-filter-bkg=value[%|k]] [-filter-nbstars=value[%|k]] [-filter-quality=value[%|k]] [-filter-incl[uded]]", process_seq_applyreg, STR_SEQAPPLYREG, TRUE, REQ_CMD_NO_THREAD},
	{"seqclean", 1, "seqclean sequencename [-reg] [-stat] [-sel]", process_seq_clean, STR_SEQCLEAN, TRUE, REQ_CMD_NONE},
	{"seqcosme", 2, "seqcosme sequencename [filename].lst [-prefix=]", process_seq_cosme, STR_SEQCOSME, TRUE, REQ_CMD_NONE},
	{"seqcosme_cfa", 2, "seqcosme_cfa sequencename [filename].lst [-prefix=]", process_seq_cosme, STR_SEQCOSME_CFA, TRUE, REQ_CMD_NONE},
	{"seqcrop", 5, "seqcrop sequencename x y width height [-prefix=]", process_seq_crop, STR_SEQCROP, TRUE, REQ_CMD_NO_THREAD},
	{"seqextract_Green", 1, "seqextract_Green sequencename [-prefix=]", process_seq_extractGreen, STR_SEQEXTRACTGREEN, TRUE, REQ_CMD_NO_THREAD},
	{"seqextract_Ha", 1, "seqextract_Ha sequencename [-prefix=]", process_seq_extractHa, STR_SEQEXTRACTHA, TRUE, REQ_CMD_NO_THREAD},
	{"seqextract_HaOIII", 1, "seqextract_HaOIII sequencename [-resample=]", process_seq_extractHaOIII, STR_SEQEXTRACTHAOIII, TRUE, REQ_CMD_NO_THREAD},
	{"seqght", 2, "seqght sequence -D= [-B=] [-LP=] [-SP=] [-HP=] [-human | -even | -independent] [channels]", process_seq_ght, STR_SEQGHT, TRUE, REQ_CMD_NONE},
	{"seqheader", 2, "seqheader sequencename keyword", process_seq_header, STR_SEQHEADER, TRUE, REQ_CMD_NONE},
	{"seqfind_cosme", 3, "seqfind_cosme sequencename cold_sigma hot_sigma [-prefix=]", process_findcosme, STR_SEQFIND_COSME, TRUE, REQ_CMD_NONE},
	{"seqfind_cosme_cfa", 3, "seqfind_cosme_cfa sequencename cold_sigma hot_sigma [-prefix=]", process_findcosme, STR_SEQFIND_COSME_CFA, TRUE, REQ_CMD_NONE},
	{"seqfindstar", 1, "seqfindstar sequencename [-layer=] [-maxstars=]", process_seq_findstar, STR_SEQFINDSTAR, TRUE, REQ_CMD_NONE},
	{"seqfixbanding", 3, "seqfixbanding sequencename amount sigma [-prefix=] [-vertical]", process_seq_fixbanding, STR_SEQFIXBANDING, TRUE, REQ_CMD_NONE},
	{"seqinvght", 2, "seqinvght sequence -D= [-B=] [-LP=] [-SP=] [-HP=] [-human | -even | -independent] [channels]", process_seq_invght, STR_SEQINVGHT, TRUE, REQ_CMD_NONE},
	{"seqinvmodasinh", 2, "seqinvmodasinh sequence -D= [-LP=] [-SP=] [-HP=] [-human | -even | -independent] [channels]", process_seq_invmodasinh, STR_SEQINVMODASINH, TRUE, REQ_CMD_NONE},
	{"seqlinstretch", 2, "seqlinstretch sequence -BP= [channels]", process_seq_linstretch, STR_SEQLINSTRETCH, TRUE, REQ_CMD_NONE},
	{"seqmerge_cfa", 2, "seqmerge_cfa sequencename bayerpattern [-prefixin=] [-prefixout=]", process_seq_merge_cfa, STR_SEQMERGE_CFA, TRUE, REQ_CMD_NO_THREAD},
	{"seqmodasinh", 2, "seqmodasinh sequence -D= [-LP=] [-SP=] [-HP=] [-human | -even | -independent] [channels]", process_seq_modasinh, STR_SEQINVMODASINH, TRUE, REQ_CMD_NONE},
	{"seqmtf", 4, "seqmtf sequencename low mid high [channels] [-prefix=]", process_seq_mtf, STR_SEQMTF, TRUE, REQ_CMD_NONE},
	{"seqpsf", 0, "seqpsf [sequencename channel { -at=x,y | -wcs=ra,dec }]", process_seq_psf, STR_SEQPSF, TRUE, REQ_CMD_NO_THREAD},
	{"seqplatesolve", 0, "seqplatesolve sequencename [image_center_coords] [-noflip] [-platesolve] [-focal=] [-pixelsize=] [-limitmag=[+-]] [-catalog=] [-downscale]", process_pcc, STR_SEQPLATESOLVE, TRUE, REQ_CMD_NO_THREAD},
	{"seqrl", 1, "seqrl sequencename [-alpha=] [-iters=] [-stop=] [-gdstep=] [-tv] [-fh] [-mul] [-seq]", process_seq_rl, STR_SEQRL, TRUE, REQ_CMD_NONE},
	{"seqsb", 1, "sb sequencename [-alpha=] [-iters=]", process_seq_sb, STR_SEQSB, TRUE, REQ_CMD_NONE},
	{"seqsplit_cfa", 1, "seqsplit_cfa sequencename [-prefix=]", process_seq_split_cfa, STR_SEQSPLIT_CFA, TRUE, REQ_CMD_NO_THREAD},
#ifdef HAVE_LIBTIFF
	{"seqstarnet", 0, "seqstarnet sequencename [-stretch] [-upscale] [-stride=value] [-nostarmask]", process_seq_starnet, STR_SEQSTARNET, TRUE, REQ_CMD_NONE},
#endif
	{"seqstat", 2, "seqstat sequencename output [option] [-cfa]", process_seq_stat, STR_SEQSTAT, TRUE, REQ_CMD_NO_THREAD},
	{"seqsubsky", 2, "seqsubsky sequencename { -rbf | degree } [-samples=20] [-tolerance=1.0] [-smooth=0.5] [-prefix=]", process_subsky, STR_SEQSUBSKY, TRUE, REQ_CMD_NONE},
	{"seqtilt", 1, "seqtilt sequencename", process_seq_tilt, STR_SEQTILT, TRUE, REQ_CMD_NO_THREAD},
	{"seqwiener", 1, "wiener sequencename [-alpha=]", process_seq_wiener, STR_SEQWIENER, TRUE, REQ_CMD_NONE},
	{"set", 1, "set { -import=inifilepath | variable=value }", process_set, STR_SET, TRUE, REQ_CMD_NONE},
	{"set16bits", 0, "set16bits", process_set_32bits, STR_SET16, TRUE, REQ_CMD_NONE},
	{"set32bits", 0, "set32bits", process_set_32bits, STR_SET32, TRUE, REQ_CMD_NONE},
	{"setcompress", 1, "setcompress 0/1 [-type=] [q]", process_set_compress, STR_SETCOMPRESS, TRUE, REQ_CMD_NONE},
#ifdef _OPENMP
	{"setcpu", 1, "setcpu number", process_set_cpu, STR_SETCPU, TRUE, REQ_CMD_NONE},
#endif
	{"setext", 1, "setext extension", process_set_ext, STR_SETEXT, TRUE, REQ_CMD_NONE},
	{"setfindstar", 0, "setfindstar [reset] [-radius=] [-sigma=] [-roundness=] [-focal=] [-pixelsize=] [-convergence=] [ [-gaussian] | [-moffat] ] [-minbeta=] [-relax=on|off] [-minA=] [-maxA=] [-maxR=]", process_set_findstar, STR_SETFINDSTAR, TRUE, REQ_CMD_NONE},
	{"setmag", 1, "setmag magnitude", process_set_mag, STR_SETMAG, FALSE, REQ_CMD_SINGLE_IMAGE},
	{"seqsetmag", 1, "seqsetmag magnitude", process_set_mag_seq, STR_SEQSETMAG, FALSE, REQ_CMD_SEQUENCE},
	{"setmem", 1, "setmem ratio", process_set_mem, STR_SETMEM, TRUE, REQ_CMD_NONE},
	{"setphot", 0, "setphot [-inner=20] [-outer=30] [-aperture=10] [-force_radius=no] [-gain=2.3] [-min_val=0] [-max_val=60000]", process_set_photometry, STR_SETPHOT, TRUE, REQ_CMD_NONE},
	{"setref", 2, "setref sequencename image_number", process_set_ref, STR_SETREF, TRUE, REQ_CMD_NONE},
	{"show", 1, "show [-clear] [{ -list=file.csv | [name] RA Dec }]", process_show, STR_SHOW, FALSE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"solsys", 0, "solsys [-mag=20.0]", process_sso, STR_SSO, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"split", 3, "split file1 file2 file3 [-hsl | -hsv | -lab]", process_split, STR_SPLIT, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_RGB | REQ_CMD_NO_THREAD},
	{"split_cfa", 0, "split_cfa", process_split_cfa, STR_SPLIT_CFA, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_FOR_CFA},
	{"stack", 1, "stack seqfilename\n"
			"stack seqfilename { sum | min | max } [filtering] [-output_norm] [-out=filename]\n"
			"stack seqfilename { med | median } [-nonorm, -norm=] [filtering] [-fastnorm] [-rgb_equal] [-output_norm] [-out=filename]\n"
			"stack seqfilename { rej | mean } [rejection type] [sigma_low sigma_high]  [-rejmap[s]] [-nonorm, -norm=] [filtering] [-fastnorm] [ -weight_from_noise | -weight_from_nbstack | -weight_from_wfwhm | -weight_from_nbstars ] [-rgb_equal] [-output_norm] [-out=filename]", process_stackone, STR_STACK, TRUE, REQ_CMD_NONE},
	{"stackall", 0, "stackall\n"
			"stackall { sum | min | max } [filtering]\n"
			"stackall { med | median } [-nonorm, norm=] [-filter-incl[uded]]\n"
			"stackall { rej | mean } [rejection type] [sigma_low sigma_high] [-nonorm, norm=] [filtering] [ -weight_from_noise | -weight_from_wfwhm | -weight_from_nbstars | -weight_from_nbstack ] [-rgb_equal] [-out=filename]", process_stackall, STR_STACKALL, TRUE, REQ_CMD_NONE},
#ifdef HAVE_LIBTIFF
	{"starnet", 0, "starnet [-stretch] [-upscale] [-stride=value] [-nostarmask]", process_starnet, STR_STARNET, TRUE, REQ_CMD_SINGLE_IMAGE},
#endif
	{"start_ls", 0, "start_ls [-dark=filename] [-flat=filename] [-rotate] [-32bits]", process_start_ls, STR_START_LS, TRUE, REQ_CMD_NO_THREAD},
	{"stat", 0, "stat [-cfa] [main]", process_stat, STR_STAT, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},
	{"stop_ls", 0, "stop_ls", process_stop_ls, STR_STOP_LS, TRUE, REQ_CMD_NONE},
	{"subsky", 1, "subsky { -rbf | degree } [-samples=20] [-tolerance=1.0] [-smooth=0.5]", process_subsky, STR_SUBSKY, TRUE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_NO_THREAD},
	{"synthstar", 0, "synthstar", process_synthstar, STR_SYNTHSTAR, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"threshlo", 1, "threshlo level", process_threshlo, STR_THRESHLO, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"threshhi", 1, "threshi level", process_threshhi, STR_THRESHHI, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"thresh", 2, "thresh lo hi", process_thresh, STR_THRESH, TRUE, REQ_CMD_SINGLE_IMAGE}, /* threshes hi and lo */
	{"tilt", 0, "tilt [clear]", process_tilt, STR_TILT, FALSE, REQ_CMD_SINGLE_IMAGE | REQ_CMD_SEQUENCE},

	{"unclipstars", 0, "unclipstars", process_unclip, STR_SYNTHSTARUNCLIP, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"unselect", 3, "unselect sequencename from to", process_unselect, STR_UNSELECT, TRUE, REQ_CMD_NONE},
	{"unsetmag", 0, "unsetmag", process_unset_mag, STR_UNSETMAG, FALSE, REQ_CMD_NONE},
	{"sequnsetmag", 0, "sequnsetmag", process_unset_mag_seq, STR_SEQUNSETMAG, FALSE, REQ_CMD_SEQUENCE},
	{"unsharp", 2, "unsharp sigma multi", process_unsharp, STR_UNSHARP, TRUE, REQ_CMD_SINGLE_IMAGE},

	{"visu", 2, "visu low high", process_visu, STR_VISU, FALSE, REQ_CMD_SINGLE_IMAGE},

	/* wavelet transform in nbr_plan plans */
	{"wavelet", 1, "wavelet nbr_layers type", process_wavelet, STR_WAVELET, TRUE, REQ_CMD_SINGLE_IMAGE},
	{"wiener", 0, "wiener [-alpha=]", process_wiener, STR_WIENER, TRUE, REQ_CMD_SINGLE_IMAGE},
	/* reconstruct from wavelet transform and weighs plans with c1, c2, c3... */
	{"wrecons", 2, "wrecons c1 c2 c3 ...", process_wrecons, STR_WRECONS, TRUE, REQ_CMD_SINGLE_IMAGE},
	EXTRA_COMMANDS

	{"",0,"",0, STR_NONE, FALSE, REQ_CMD_NONE}
};

#endif /* SRC_CORE_COMMAND_LIST_H_ */
