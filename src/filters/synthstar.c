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

#include <assert.h>
#include <math.h>
#include "core/siril.h"
#include "core/gui_iface.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "algos/colors.h"
#include "algos/median_fast.h"
#include "algos/star_finder.h"
#include "algos/PSF.h"
#include "algos/extraction.h"
#include "algos/siril_random.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "filters/synthstar.h"
#include "opencv/opencv.h"
#include "core/op_descriptors.h"
#include "core/nde/nde_history.h"

/* ---- captured effective star list (NDE phase 4.5 Convention 2) ------------
 *
 * synthstar and unpurple both derive their output from com.stars (or, when the
 * list is empty, from a fresh star detection that reads com.pref — so param-only
 * replay would be nondeterministic).  To make replay deterministic we stash the
 * EFFECTIVE star list that was actually consumed as a `stars=` blob (BGE
 * effective_samples pattern): a `:`-separated list of per-star tuples, fields in
 * the FIXED order `x,y,A,fwhmx,fwhmy,beta,sat,profile` with %.9g floats and ints
 * for sat (has_saturated) and profile (starprofile).  Only the fields the
 * consuming hooks read are stored.  Shared by synthstar (which owns the codec)
 * and unpurple (which includes synthstar.h). */

/* Serialize the consumed star list into "x,y,A,fwhmx,fwhmy,beta,sat,profile"
 * tuples joined by ':'.  Returns NULL when the list is empty. */
gchar *synthstar_stars_to_blob(psf_star **stars, int nb_stars) {
	if (!stars || nb_stars < 1)
		return NULL;
	GString *s = g_string_new(NULL);
	char b[8][G_ASCII_DTOSTR_BUF_SIZE];
	gboolean first = TRUE;
	for (int i = 0; i < nb_stars; i++) {
		const psf_star *st = stars[i];
		if (!st)
			continue;
		g_ascii_formatd(b[0], sizeof b[0], "%.9g", st->xpos);
		g_ascii_formatd(b[1], sizeof b[1], "%.9g", st->ypos);
		g_ascii_formatd(b[2], sizeof b[2], "%.9g", st->A);
		g_ascii_formatd(b[3], sizeof b[3], "%.9g", st->fwhmx);
		g_ascii_formatd(b[4], sizeof b[4], "%.9g", st->fwhmy);
		g_ascii_formatd(b[5], sizeof b[5], "%.9g", st->beta);
		g_string_append_printf(s, "%s%s,%s,%s,%s,%s,%s,%d,%d",
				first ? "" : ":", b[0], b[1], b[2], b[3], b[4], b[5],
				st->has_saturated ? 1 : 0, (int)st->profile);
		first = FALSE;
	}
	if (!s->len) {
		g_string_free(s, TRUE);
		return NULL;
	}
	return g_string_free(s, FALSE);
}

/* Reconstruct a NULL-terminated psf_star array from a `stars=` blob; *nb_out
 * (if non-NULL) receives the count.  NULL for an empty/absent blob.  Free with
 * free_fitted_stars().  Only the stashed fields are set; everything else is
 * zero (calloc) — the hooks read exactly the stashed fields. */
psf_star **synthstar_stars_from_blob(const char *blob, int *nb_out) {
	if (nb_out)
		*nb_out = 0;
	if (!blob || !*blob)
		return NULL;
	gchar **tuples = g_strsplit(blob, ":", -1);
	guint n = g_strv_length(tuples);
	if (n < 1) {
		g_strfreev(tuples);
		return NULL;
	}
	psf_star **stars = malloc((size_t)(n + 1) * sizeof(psf_star *));
	if (!stars) {
		g_strfreev(tuples);
		return NULL;
	}
	int count = 0;
	for (guint i = 0; i < n; i++) {
		gchar **f = g_strsplit(tuples[i], ",", 8);
		if (g_strv_length(f) == 8) {
			psf_star *st = new_psf_star();
			if (st) {
				st->xpos  = g_ascii_strtod(f[0], NULL);
				st->ypos  = g_ascii_strtod(f[1], NULL);
				st->A     = g_ascii_strtod(f[2], NULL);
				st->fwhmx = g_ascii_strtod(f[3], NULL);
				st->fwhmy = g_ascii_strtod(f[4], NULL);
				st->beta  = g_ascii_strtod(f[5], NULL);
				st->has_saturated = (g_ascii_strtoll(f[6], NULL, 10) != 0);
				st->profile = (starprofile)g_ascii_strtoll(f[7], NULL, 10);
				stars[count++] = st;
			}
		}
		g_strfreev(f);
	}
	stars[count] = NULL;
	if (count < 1) {
		free(stars);
		return NULL;
	}
	if (nb_out)
		*nb_out = count;
	return stars;
}

/* ---- detection-parameter codec (DELEGATED provenance) ---------------------
 * Serialize the star_finder_params that a delegated (auto-detect) star op
 * consumed, so replay re-detects with the SAME parameters rather than the
 * user's current preferences.  All fields are detection-affecting. */
void synthstar_conf_to_kv(GString *kv, const star_finder_params *sf) {
	nde_kv_add_int(kv, "sf_radius", sf->radius);
	nde_kv_add_double(kv, "sf_sigma", sf->sigma);
	nde_kv_add_double(kv, "sf_roundness", sf->roundness);
	nde_kv_add_double(kv, "sf_focal", sf->focal_length);
	nde_kv_add_double(kv, "sf_pixsize", sf->pixel_size_x);
	nde_kv_add_int(kv, "sf_convergence", sf->convergence);
	nde_kv_add_bool(kv, "sf_relax", sf->relax_checks);
	/* on-disk value: enum order is frozen by the NDE format — do not reorder */
	nde_kv_add_int(kv, "sf_profile", sf->profile);
	nde_kv_add_double(kv, "sf_min_beta", sf->min_beta);
	nde_kv_add_double(kv, "sf_min_A", sf->min_A);
	nde_kv_add_double(kv, "sf_max_A", sf->max_A);
	nde_kv_add_double(kv, "sf_max_r", sf->max_r);
}

gboolean synthstar_conf_from_kv(GHashTable *kv, star_finder_params *out) {
	gint64 radius, convergence, profile;
	gboolean relax;
	if (!nde_kv_get_int(kv, "sf_radius", &radius) ||
	    !nde_kv_get_double(kv, "sf_sigma", &out->sigma) ||
	    !nde_kv_get_double(kv, "sf_roundness", &out->roundness) ||
	    !nde_kv_get_double(kv, "sf_focal", &out->focal_length) ||
	    !nde_kv_get_double(kv, "sf_pixsize", &out->pixel_size_x) ||
	    !nde_kv_get_int(kv, "sf_convergence", &convergence) ||
	    !nde_kv_get_bool(kv, "sf_relax", &relax) ||
	    !nde_kv_get_int(kv, "sf_profile", &profile) ||
	    !nde_kv_get_double(kv, "sf_min_beta", &out->min_beta) ||
	    !nde_kv_get_double(kv, "sf_min_A", &out->min_A) ||
	    !nde_kv_get_double(kv, "sf_max_A", &out->max_A) ||
	    !nde_kv_get_double(kv, "sf_max_r", &out->max_r))
		return FALSE;
	out->radius = (int)radius;
	out->convergence = (int)convergence;
	out->relax_checks = relax;
	out->profile = (starprofile)profile;
	return TRUE;
}

/* Minimal destructor-first params struct for synthstar (was paramless).  Only
 * carries the stashed effective star list so replay can reinstall it. */
static void free_synthstar_data(void *p) {
	struct synthstar_data *d = p;
	if (!d)
		return;
	g_free(d->stars_blob);
	free(d);
}

struct synthstar_data *new_synthstar_data(void) {
	struct synthstar_data *d = calloc(1, sizeof(struct synthstar_data));
	if (d)
		d->destroy_fn = free_synthstar_data;
	return d;
}

/* Emit either the EXPLICIT pinned list or the DELEGATED detection parameters
 * (never both — provenance is decided at capture). */
static void synthstar_provenance_to_kv(GString *kv, const struct synthstar_data *d) {
	if (!d)
		return;
	if (d->star_auto) {
		nde_kv_add_bool(kv, "stars_auto", TRUE);
		synthstar_conf_to_kv(kv, &d->star_conf);
	} else if (d->stars_blob && *d->stars_blob) {
		nde_kv_add_str(kv, "stars", d->stars_blob);
	}
}

/* Read provenance into @d.  FALSE (→ non-replayable) when neither an explicit
 * stars= list nor a complete delegated conf is present. */
static gboolean synthstar_provenance_from_kv(GHashTable *kv, struct synthstar_data *d) {
	gboolean auto_d = FALSE;
	nde_kv_get_bool(kv, "stars_auto", &auto_d);
	if (auto_d) {
		if (!synthstar_conf_from_kv(kv, &d->star_conf))
			return FALSE;
		d->star_auto = TRUE;
		return TRUE;
	}
	const char *stars = nde_kv_get_str(kv, "stars");
	if (!stars || !*stars)
		return FALSE;
	d->stars_blob = g_strdup(stars);
	return TRUE;
}

static gchar *synthstar_serialize(gconstpointer user) {
	GString *kv = nde_kv_start();
	synthstar_provenance_to_kv(kv, user);
	return nde_kv_end(kv);
}

static gpointer synthstar_deserialize(const gchar *blob, int version) {
	if (version > op_desc_synthstar.version)
		return NULL;
	GHashTable *kv = nde_kv_parse(blob);
	struct synthstar_data *d = new_synthstar_data();
	/* Neither provenance present ⇒ a pre-feature record (or corrupt): stays
	 * honestly non-replayable. */
	if (d && !synthstar_provenance_from_kv(kv, d)) {
		free_synthstar_data(d);
		d = NULL;
	}
	g_hash_table_unref(kv);
	return d;
}

/* replay_pre: EXPLICIT ⇒ reinstall com.stars from the pinned blob so the hook
 * consumes it unchanged.  DELEGATED ⇒ nothing to do here: the hook re-detects
 * with the recorded conf (passed via the data struct), ignoring com.stars.
 * Replay runs in the exclusive job slot, so no concurrent star op is mid-
 * flight. */
static int synthstar_replay_pre(gpointer user, GHashTable *kv, fits *target) {
	(void)target; (void)kv;
	struct synthstar_data *d = user;
	if (d->star_auto)
		return 0;
	int nb = 0;
	psf_star **arr = synthstar_stars_from_blob(d->stars_blob, &nb);
	if (!arr || nb < 1)
		return 1;
	replace_com_stars(arr);   /* installs + frees any previous list */
	return 0;
}

/* Op descriptors — single source of truth for these ops (op_descriptor.h) */
const op_descriptor op_desc_synthstar = {
	.id = "star.synthstar", .version = 1,
	.image_hook = synthstar_image_hook,
	.log_hook = synthstar_log_hook,
	.description = N_("Synthetic stars"),
	.mem_ratio = 0.0f,
	.flags = 0,
	.serialize = synthstar_serialize, .deserialize = synthstar_deserialize,
	.replay_pre = synthstar_replay_pre,
};

/* star.unclip is Tier A via DELEGATED provenance: reprofile_saturated_stars()
 * always auto-detects per channel (it never consumes com.stars), reading
 * com.pref.starfinder_conf.  We record those detection parameters and re-run
 * detection at replay — faithful to "find stars automatically" intent and
 * re-flowing upstream amends — rather than pinning a star list.  No
 * replay_pre: the conf is passed to the hook via the data struct. */
static gchar *unclip_serialize(gconstpointer user) {
	const struct synthstar_data *d = user;
	GString *kv = nde_kv_start();
	nde_kv_add_bool(kv, "stars_auto", TRUE);
	if (d)
		synthstar_conf_to_kv(kv, &d->star_conf);
	return nde_kv_end(kv);
}

static gpointer unclip_deserialize(const gchar *blob, int version) {
	if (version > op_desc_unclip.version)
		return NULL;
	GHashTable *kv = nde_kv_parse(blob);
	struct synthstar_data *d = new_synthstar_data();
	if (d) {
		d->star_auto = TRUE;
		if (!synthstar_conf_from_kv(kv, &d->star_conf)) {
			free_synthstar_data(d);
			d = NULL;
		}
	}
	g_hash_table_unref(kv);
	return d;
}

const op_descriptor op_desc_unclip = {
	.id = "star.unclip", .version = 1,
	.image_hook = unclip_image_hook,
	.log_hook = unclip_log_hook,
	.description = N_("Unclip stars"),
	.mem_ratio = 0.0f,
	.flags = 0,
	.serialize = unclip_serialize, .deserialize = unclip_deserialize,
};

int generate_synthstars(fits *fit, struct synthstar_data *prov,
                        const star_finder_params *conf_override);
int reprofile_saturated_stars(fits *fit, star_finder_params *conf_out,
                              const star_finder_params *conf_override);

void makeairy(float *psf, const int size, const float lum, const float xoff, const float yoff, float wavelength, float aperture, float focal_length, float pixel_size, float obstruction) {
	// wavelength is given in nm; pixel size in microns; aperture and focal length are given in mm. Convert all to metres for consistency
	// obstruction is central obstruction ratio given between [0, 1[
	wavelength *= 1.e-9f;
	aperture *= 1.e-3f;
	focal_length *= 1.e-3f;
	pixel_size *= 1.e-6f;
	const int halfpsfdim = (size - 1) / 2;
	// obstruction = 0.5; // for testing purposes - can be removed
	float obscorr = (obstruction > 0.f) ? 1.f / pow(1 - obstruction * obstruction, 2.f) : 1.f;

	// Following the formulae at the Wikipedia "Airy disk" article
	const float constant = (2.f * G_PI * (aperture / 2.f) / wavelength) * (1.f / focal_length);
	for (int x = -halfpsfdim; x <= halfpsfdim; x++) {
		for (int y = -halfpsfdim; y <= halfpsfdim; y++) {
			float xf = (x - xoff + 0.5f) * pixel_size;
			float yf = (y - yoff - 0.5f) * pixel_size;
			float q = constant * sqrtf(xf * xf + yf * yf);
			float bessel = j1(q);
			if (obstruction == 0.f) {
				psf[(x + halfpsfdim) + (y + halfpsfdim) * size] = (q != 0.f) ? lum * pow(2.f * bessel / q, 2.f) : lum;
			} else {
				psf[(x + halfpsfdim) + (y + halfpsfdim) * size] = (q != 0.f) ? lum * obscorr * pow(2.f / q * (bessel - obstruction * j1(obstruction * q)), 2.f) : lum;
			}
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

void makemoffat(float *psf, const int size, const float fwhm, const float lum, const float xoff,
				const float yoff, const float beta, const float ratio, const float angle) {
	float anglerad = angle * G_PI / 180.f;
	const float alpha = 0.6667f * fwhm;
	const float alphax = alpha;
	const float alphay = alpha / ratio;
	const int halfpsfdim = (size - 1) / 2;
	float a = powf(cosf(anglerad)/alphax, 2.f) + powf(sinf(anglerad)/alphay, 2.f);
	float b = powf(sinf(anglerad)/alphax, 2.f) + powf(cosf(anglerad)/alphay, 2.f);
	float c = 2.f * sinf(anglerad) * cosf(anglerad) * (1.f/(alphax * alphax) - 1.f/(alphay * alphay));
	for (int x = -halfpsfdim; x <= halfpsfdim; x++) {
		for (int y = -halfpsfdim; y <= halfpsfdim; y++) {
			float xf = (x - xoff + 0.5f);
			float yf = (y - yoff - 0.5f);
			psf[(x + halfpsfdim) + ((y + halfpsfdim) * size)] = lum
					* powf(1.0f + (a * xf * xf) + (b * yf * yf) + (c * xf * yf),
							-beta);
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

void makegaussian(float *psf, int size, float fwhm, float lum, float xoffset, float yoffset, float ratio, float angle) {
	int halfpsfdim = (size - 1) / 2;
	float anglerad = angle * G_PI / 180.f;
	float sigmax = fwhm / _2_SQRT_2_LOG2;
	float sigmay = fwhm / (ratio * _2_SQRT_2_LOG2);
	float tssx = 2 * sigmax * sigmax;
	float tssy = 2 * sigmay * sigmay;
	float a = powf(cosf(anglerad), 2.f) / tssx + powf(sinf(anglerad), 2.f) / tssy;
	float b = sinf(2 * anglerad) / (2 * tssx) - sinf(2 * anglerad) / (2 * tssy);
	float c = powf(sinf(anglerad), 2.f) / tssx + powf(cosf(anglerad), 2.f) / tssy;
	for (int x = -halfpsfdim; x <= halfpsfdim; x++) {
		for (int y = -halfpsfdim; y <= halfpsfdim; y++) {
			float xf = (x - xoffset + 0.5f);
			float yf = (y - yoffset - 0.5f);
			psf[(x + halfpsfdim) + ((y + halfpsfdim) * size)] = lum
					* expf(-(a * xf * xf + 2 * b * xf * yf + c * yf * yf));
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

void makedisc(float *psf, int size, float width, float lum, float xoffset, float yoffset) {
	int halfpsfdim = (size - 1) / 2;
	float radius = width / 2.f;
	// maxranditer big enough to get a good random survey of each pixel,
	// not enough to cause slowness with large kernels.
	const int maxranditer = 10000;
	float radiussq = radius * radius;
	float solidradsq = powf(radius - 0.5f, 2.f); // radius less than which pixel value is 1.f
	float zeroradsq =powf(radius + 0.5f, 2.f); // radius greater than which pixel value is 0.f
	for (int x = -halfpsfdim; x <= halfpsfdim; x++) {
		for (int y = -halfpsfdim; y <= halfpsfdim; y++) {
			float pixradsq = powf((x - xoffset + 0.5f), 2.f) + powf((y - yoffset - 0.5f), 2.f);
			if (pixradsq < solidradsq) {
				psf[(x + halfpsfdim) + ((y + halfpsfdim) * size)] = 1.f;
			} else if (pixradsq > zeroradsq) {
				psf[(x + halfpsfdim) + ((y + halfpsfdim) * size)] = 0.f;
			} else {
				int count = 0;
				for (int randiter = 0 ; randiter < maxranditer; randiter++) {
					float xrandoff = siril_random_float();
					float yrandoff = siril_random_float();
					float xf = x - xoffset + xrandoff + 0.5f;
					float yf = y - yoffset + yrandoff - 0.5f;
					if ((xf * xf + yf * yf) < radiussq)
						count++;
				}
				psf[(x + halfpsfdim) + ((y + halfpsfdim) * size)] = (lum * count) / maxranditer;
			}
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

static void add_star_to_rgb_buffer(const float *H, const float *S, const float *psfL, int size, float *Hsynth, float *Ssynth, float *Lsynth, int x, int y, int dimx, int dimy) {
	int halfpsfdim = (size - 1) / 2;
	int xx, yy;
#define EPSILON 1e-30
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static) collapse(2) num_threads(com.max_thread) private(xx, yy) if(com.max_thread > 1)
#endif
	for (int psfx = 0; psfx < size; psfx++) {
		for (int psfy = 0; psfy < size; psfy++) {
			xx = x + psfx - halfpsfdim;
			yy = y + psfy - halfpsfdim;
			if (xx > 0 && xx < dimx && yy > 0 && yy < dimy) {
				Hsynth[xx + ((dimy - yy) * dimx)] =
						H[xx + ((dimy - yy) * dimx)];
				Ssynth[xx + ((dimy - yy) * dimx)] =
						S[xx + ((dimy - yy) * dimx)];
				Lsynth[xx + ((dimy - yy) * dimx)] += psfL[psfx + (psfy * size)];
			}
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
	return;
}

static void add_star_to_mono_buffer(const float *psfL, int size, float *Lsynth, int x, int y, int dimx, int dimy) {
	const int halfpsfdim = (size - 1) / 2;
	int xx, yy;
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static) collapse(2) num_threads(com.max_thread) private(xx, yy) if(com.max_thread > 1)
#endif
	for (int psfx = 0; psfx < size; psfx++) {
		for (int psfy = 0; psfy < size; psfy++) {
			xx = x + psfx - halfpsfdim;
			yy = y + psfy - halfpsfdim;
			if (xx > 0 && xx < dimx && yy > 0 && yy < dimy) {
				Lsynth[xx + ((dimy - yy) * dimx)] += psfL[psfx + (psfy * size)];
			}
		}
	}
#ifdef _OPENMP
#pragma omp barrier
#endif

	return;
}

static void replace_sat_star_in_buffer(const float *psfL, int size, float *Lsynth, int x, int y, int dimx, int dimy, float sat, float bg, float noise) {
	float* buf = calloc(1, size * size * sizeof(float));
	float* resbuf = malloc(size * size * sizeof(float));
	const int halfpsfdim = (size - 1) / 2;
	int xx, yy;
// Blend synthetic data into Lsynth and make a copy in buf for filtering the join
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
{
#pragma omp for simd schedule(static) collapse(2) private(xx, yy)
#endif
	for (int psfx = 0; psfx < size; psfx++) {
		for (int psfy = 0; psfy < size; psfy++) {
			xx = x + psfx - halfpsfdim;
			yy = y + psfy - halfpsfdim;
// Note the bounds, xx >= 0 and < dimx but yy > 0 and <= dimy
// This is correct, it is because Lsynth gets indexed by xx but by dimy - yy
// Same comment applies below in the copy back
			if (xx >= 0 && xx < dimx && yy > 0 && yy <= dimy) {
				float orig = Lsynth[xx + ((dimy - yy) * dimx)];
				float synth = psfL[psfx + (psfy * size)];
				float synthfactor = (orig < sat) ? 1.f -
					((sat - orig) / sat) : 1.f;
				Lsynth[xx + ((dimy - yy) * dimx)] += max(synth * synthfactor
					+ orig * (1 - synthfactor) - sat, 0.f);
				buf[psfx + psfy * size] = Lsynth[xx + ((dimy - yy) * dimx)];
			}
		}
	}
	memcpy(resbuf, buf, size * size * sizeof(float));

// Carry out median blur of middle part, storing the result in resbuf
// in order not to overwrite data in buf that is still needed as input
// for remaining pixel calculations
#ifdef _OPENMP
#pragma omp for simd schedule(static) collapse(2)
#endif
	for (int i = halfpsfdim * 2/3 ; i < halfpsfdim * 4/3 ; i++) {
		for (int j = halfpsfdim * 2/3 ; j < halfpsfdim * 4/3 ; j++) {
			int il = i - 1;
			int iu = i + 1;
			int jl = (j - 1) * size;
			int jm = j * size;
			int ju = (j + 1) * size;
			int offx = i - halfpsfdim;
			int offy = j - halfpsfdim;
			int rad = halfpsfdim * 2/3;
			// Only blur within a circle of radius rad
			if (offx * offx + offy * offy <= rad * rad)
				resbuf[i + jm] = median9f(
					buf[il + jl],
					buf[i + jl],
					buf[iu + jl],
					buf[il + jm],
					buf[i + jm],
					buf[iu + jm],
					buf[il + ju],
					buf[i + ju],
					buf[iu + ju]);
		}
	}

// Copy resbuf back into Lsynth
#ifdef _OPENMP
#pragma omp for simd schedule(static) collapse(2)
#endif
	for (int psfx = size * 2/3; psfx < size * 2/3; psfx++) {
		for (int psfy = size * 2/3; psfy < size * 4/3; psfy++) {
			xx = x + psfx - halfpsfdim;
			yy = y + psfy - halfpsfdim;
			if (xx >= 0 && xx < dimx && yy > 0 && yy <= dimy) {
				Lsynth[xx + ((dimy - yy) * dimx)] = resbuf[psfx + psfy * size];
			}
		}
	}
#ifdef _OPENMP
#pragma omp barrier
}
#endif
	free(buf);
	free(resbuf);
	return;
}

int starcount(psf_star **stars) {
	int i = 0;
	if (!(stars)) {
		return 0;
	} else {
		while (stars[i]) {
			i++;
		}
	}
	return i;
}

/* generic_image_worker hooks */
int synthstar_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct synthstar_data *d = (struct synthstar_data *)args->user;
	/* Capture (nde_replay FALSE): generate_synthstars records the provenance
	 * into @d (EXPLICIT pinned list, or DELEGATED conf) for the worker's NDE
	 * capture that runs after the hook returns.  DELEGATED replay: pass the
	 * recorded conf as an override so detection reproduces it. */
	const star_finder_params *ov =
			(args->nde_replay && d && d->star_auto) ? &d->star_conf : NULL;
	return generate_synthstars(fit, d, ov) < 0 ? 1 : 0;
}

gchar *synthstar_log_hook(gpointer p, log_hook_detail detail) {
	return g_strdup(_("Synthetic stars"));
}

int unclip_image_hook(struct generic_img_args *args, fits *fit, int threads) {
	struct synthstar_data *d = (struct synthstar_data *)args->user;
	/* unclip always auto-detects (DELEGATED).  Capture: record the conf used
	 * so the serializer can emit it.  Replay: reproduce it via the override. */
	const star_finder_params *ov =
			(args->nde_replay && d && d->star_auto) ? &d->star_conf : NULL;
	star_finder_params used = { 0 };
	int rc = reprofile_saturated_stars(fit, d ? &used : NULL, ov);
	if (d && !args->nde_replay) {
		d->star_auto = TRUE;
		d->star_conf = used;
	}
	return rc < 0 ? 1 : 0;
}

gchar *unclip_log_hook(gpointer p, log_hook_detail detail) {
	return g_strdup(_("Unclip stars"));
}

int generate_synthstars(fits *fit, struct synthstar_data *prov,
                        const star_finder_params *conf_override) {
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	gui_iface.set_progress(PROGRESS_RESET, _("Star synthesis (full star mask creation): processing..."));
	gboolean is_RGB = TRUE;
	gboolean is_32bit = TRUE;
	gboolean stars_needs_freeing = FALSE;
	float norm = 1.0f, invnorm = 1.0f;
	int nb_stars = 0;
	psf_star **stars = NULL;
	/* DELEGATED replay: ignore com.stars and re-detect with the recorded conf
	 * (installed transiently below).  Capture / EXPLICIT replay: consume
	 * com.stars, auto-detecting only when it is empty. */
	gboolean force_detect = (conf_override != NULL);
	gboolean auto_detected = FALSE;

	if (!force_detect) {
		// Private, reader-locked copy of com.stars: the star-rendering loop below
		// runs on a worker thread and must not deref a list another thread may free.
		stars = snapshot_com_stars(&nb_stars);
		if (stars)
			stars_needs_freeing = TRUE;
	}
	int comstar_count = nb_stars;   // 0 when force_detect

	if (comstar_count < 1) {
		auto_detected = TRUE;
		// snapshot_com_stars() can return a non-NULL but empty array (first
		// duplicate_psf OOM); free it before findstar_worker overwrites stars.
		if (stars_needs_freeing) {
			free_fitted_stars(stars);
			stars = NULL;
			stars_needs_freeing = FALSE;
		}
		// Set up starfinder_data structure
		struct starfinder_data *sf_data = calloc(1, sizeof(struct starfinder_data));
		if (!sf_data) {
			siril_log_error(_("Memory allocation failed\n"));
			gui_iface.set_progress(PROGRESS_RESET, PROGRESS_TEXT_RESET);
			// snapshot already freed above; stars_needs_freeing is FALSE here.
			return -1;
		}

		sf_data->im.fit = fit;
		sf_data->im.from_seq = NULL;
		sf_data->im.index_in_seq = -1;
		sf_data->layer = (fit->naxes[2] == 1) ? 0 : 1;
		sf_data->max_stars_fitted = MAX_STARS;
		sf_data->selection = (rectangle){0, 0, 0, 0}; // no selection
		sf_data->save_eqcoords = FALSE;
		sf_data->ref_wcs = NULL;
		sf_data->stars = &stars;
		sf_data->nb_stars = &nb_stars;
		sf_data->threading = MULTI_THREADED;
		sf_data->update_GUI = FALSE;
		sf_data->process_all_images = FALSE;
		sf_data->already_in_thread = TRUE;
		sf_data->keep_stars = FALSE;

		// Call the worker function (with the recorded conf installed when replaying)
		int retval = 0;
		if (conf_override)
			WITH_STARFINDER_CONF(conf_override,
					retval = GPOINTER_TO_INT(findstar_worker(sf_data)));
		else
			retval = GPOINTER_TO_INT(findstar_worker(sf_data));
		free(sf_data);

		if (retval != 0 || !stars) {
			siril_log_error(_("Star detection failed\n"));
			gui_iface.set_progress(PROGRESS_RESET, PROGRESS_TEXT_RESET);
			if (stars)
				free_fitted_stars(stars);
			return -1;
		}
		stars_needs_freeing = TRUE;
	}

	if (nb_stars < 1 || !stars) {
		siril_log_error(_("No stars detected in the image.\n"));
		if (stars_needs_freeing)
			free_fitted_stars(stars);
		gui_iface.set_progress(PROGRESS_RESET, PROGRESS_TEXT_RESET);
		return -1;
	} else {
		siril_log_message(_("Synthesizing %d stars...\n"), nb_stars);
	}

	/* NDE Convention 2: record provenance for replay (capture only — an
	 * override implies we are already replaying).  EXPLICIT (com.stars used):
	 * pin the effective list.  DELEGATED (auto-detected): record the detection
	 * parameters and DROP the list, so replay re-detects against amended
	 * upstream pixels. */
	if (prov && !force_detect) {
		if (auto_detected) {
			prov->star_auto = TRUE;
			prov->star_conf = com.pref.starfinder_conf;
			g_clear_pointer(&prov->stars_blob, g_free);
		} else {
			prov->star_auto = FALSE;
			g_free(prov->stars_blob);
			prov->stars_blob = synthstar_stars_to_blob(stars, nb_stars);
		}
	}

	if (fit->type == DATA_USHORT) {
		is_32bit = FALSE;
		norm = get_normalized_value(fit);
		invnorm = 1.0f / norm;
	}
	if (fit->naxes[2] != 3)
		is_RGB = FALSE;

	int dimx = fit->naxes[0];
	int dimy = fit->naxes[1];
	int npixels = dimx * dimy;
	gboolean buf_needs_freeing = FALSE;
	// Regardless of 16/32bit store the data in a buffer, converting if needed
	float *buf[3];
	if (is_RGB) {
		if (is_32bit) {
			buf[RLAYER] = fit->fpdata[RLAYER];
			buf[GLAYER] = fit->fpdata[GLAYER];
			buf[BLAYER] = fit->fpdata[BLAYER];
		} else {
			buf[RLAYER] = (float*) calloc(npixels, sizeof(float));
			buf[GLAYER] = (float*) calloc(npixels, sizeof(float));
			buf[BLAYER] = (float*) calloc(npixels, sizeof(float));
			buf_needs_freeing = TRUE;
			for (size_t i = 0; i < npixels; i++) {
				buf[RLAYER][i] = (float) fit->pdata[RLAYER][i] * invnorm;
				buf[GLAYER][i] = (float) fit->pdata[GLAYER][i] * invnorm;
				buf[BLAYER][i] = (float) fit->pdata[BLAYER][i] * invnorm;
			}
		}
	} else { // mono
		if (is_32bit)
			buf[RLAYER] = fit->fdata;
		else {
			buf[RLAYER] = (float*) calloc(npixels, sizeof(float));
			buf_needs_freeing = TRUE;
			for (size_t i = 0; i < npixels; i++)
				buf[RLAYER][i] = (float) fit->data[i] * invnorm;
		}
	}

	// Normalize the buffer to avoid issues with colorspace conversion
	float bufmax = 1.f;
	for (size_t chan = 0; chan < fit->naxes[2]; chan++)
		for (size_t i = 0; i < npixels; i++)
			if (buf[chan][i] > bufmax)
				bufmax = buf[chan][i];
	for (size_t chan = 0; chan < fit->naxes[2]; chan++)
		for (size_t i = 0; i < npixels; i++)
			buf[chan][i] /= bufmax;

	float *H = NULL, *S = NULL, *Hsynth = NULL, *Ssynth = NULL, *Lsynth, junk;
	Lsynth = (float*) calloc(npixels, sizeof(float));

	// For RGB images, convert pixel colour data from fit into H and S arrays. L is irrelevant as we will synthesize L.
	if (is_RGB) {
		H = (float*) calloc(npixels, sizeof(float));
		S = (float*) calloc(npixels, sizeof(float));
		Hsynth = (float*) calloc(npixels, sizeof(float));
		Ssynth = (float*) calloc(npixels, sizeof(float));
		for (size_t i = 0; i < npixels; i++) {
			rgb_to_hsl_float_sat(buf[RLAYER][i], buf[GLAYER][i],
					buf[BLAYER][i], 0.f, &H[i], &S[i], &junk);
		}
	}
	// Calculate average Moffat beta
	size_t moffat_count = 0;
	double avg_moffat_beta = 0.;

	gboolean stopcalled = FALSE;
	// Synthesize a PSF for each star in the star array s, based on its measured parameters
	gboolean gaussian = TRUE;
	if (stars[0]->profile == PSF_MOFFAT_BFREE)
		gaussian = FALSE;
	if (!gaussian) {
		for (size_t n = 0 ; n < nb_stars ; n++) {
			moffat_count++;
			avg_moffat_beta += stars[n]->beta;
		}
		avg_moffat_beta /= moffat_count;
		siril_log_debug("# Moffat profile stars: %zd, average beta = %.3f\n", moffat_count, avg_moffat_beta);
	}
	for (int n = 0; n < nb_stars; n++) {
		// Check if stop has been pressed
		if (!processing_should_continue())
			stopcalled = TRUE;
		gui_iface.set_progress((double) n / (double) nb_stars, NULL);
		if (!stopcalled) {
			float lum = (float) stars[n]->A;
			if (lum < 0.0f)
				lum = 0.0f;
			if (!is_32bit)
				lum *= invnorm;
			assert(lum >= 0.0f);
			float xoff = (float) stars[n]->xpos - (int) stars[n]->xpos;
			float yoff = (float) stars[n]->ypos - (int) stars[n]->ypos;
			int size = (int) 5 * max(stars[n]->fwhmx, stars[n]->fwhmy); // This is big enough that even under extreme stretching the synthesized psf tails off smoothly
			if (!gaussian)
				size *= 10 / stars[n]->beta; // Increase the PSF size markedly for low beta stars
			if (!(size % 2))
				size++;
			if (size > 1024) // protect against excessive memory allocations due to bad parameters;
							 // 100px should be more than enough for the fwhm of even a very saturated star
				continue;
			float minfwhm = min(stars[n]->fwhmx, stars[n]->fwhmy);

			// Synthesize the luminance profile and add to the star mask in HSL colourspace
			float *psfL = (float*) calloc(size * size, sizeof(float));
			if (!psfL) // May happen if size is excessively large because of a bad fwhm value
				continue;
			float beta = 8.f;
			if (!gaussian) {
				if (stars[n]->beta > 0.0) {
					beta=stars[n]->beta;
					minfwhm = min(stars[n]->fwhmx, stars[n]->fwhmy);
				} else if (moffat_count > 0)
					beta = avg_moffat_beta;
			}
			if (stars[n]->has_saturated || gaussian)
				makegaussian(psfL, size, minfwhm, lum, xoff, yoff, 1.f, 0.f);
			else
				makemoffat(psfL, size, minfwhm, lum, xoff, yoff, beta, 1.f, 0.f);
			if (is_RGB)
				add_star_to_rgb_buffer(H, S, psfL, size, Hsynth, Ssynth, Lsynth,
						(float) stars[n]->xpos, (float) stars[n]->ypos, dimx, dimy);
			else
				add_star_to_mono_buffer(psfL, size, Lsynth, (float) stars[n]->xpos,
						(float) stars[n]->ypos, dimx, dimy);
			free(psfL);
		}
	}
	// Stars are only freed if they were *not* taken from com.stars: if the
	// user has made a specific selection of stars, we want to leave that
	// selection intact.
	if (stars_needs_freeing)
		free_fitted_stars(stars);

	// Construct the RGB from synthetic L (and for RGB images, also the H and S values from the orginal image thus giving our synthesized stars the correct colour)
	if (!stopcalled) {
		if (is_RGB) {
			float *R, *G, *B;
			R = (float*) calloc(npixels, sizeof(float));
			G = (float*) calloc(npixels, sizeof(float));
			B = (float*) calloc(npixels, sizeof(float));
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread) if (com.max_thread > 1)
{
#endif
			float bufmaxx = 1.f;
			for (size_t i = 0; i < npixels; i++)
				if (Lsynth[i] > bufmaxx)
					bufmaxx = Lsynth[i];
#ifdef _OPENMP
#pragma omp for simd schedule(static)
#endif
			for (size_t i = 0; i < npixels; i++)
				Lsynth[i] /= bufmaxx;

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (size_t n = 0; n < npixels; n++) {
				hsl_to_rgb_float_sat(Hsynth[n], Ssynth[n], Lsynth[n], &R[n],
						&G[n], &B[n]);
				// Trap NaNs and infinities
				R[n] = (isnan(R[n]) || isinf(R[n]) || R[n] < 0.f) ? 0.f : R[n];
				G[n] = (isnan(G[n]) || isinf(G[n]) || G[n] < 0.f) ? 0.f : G[n];
				B[n] = (isnan(B[n]) || isinf(B[n]) || B[n] < 0.f) ? 0.f : B[n];
			}
			if (is_32bit) {
#ifdef _OPENMP
#pragma omp for simd schedule(static)
#endif
				for (size_t n = 0; n < npixels; n++) {
					fit->fpdata[RLAYER][n] = R[n];
					fit->fpdata[GLAYER][n] = G[n];
					fit->fpdata[BLAYER][n] = B[n];
				}
			} else {
#ifdef _OPENMP
#pragma omp for simd schedule(static)
#endif
				for (size_t n = 0; n < npixels; n++) {
					fit->pdata[RLAYER][n] = roundf_to_WORD(R[n] * norm);
					fit->pdata[GLAYER][n] = roundf_to_WORD(G[n] * norm);
					fit->pdata[BLAYER][n] = roundf_to_WORD(B[n] * norm);
				}
			}
#ifdef _OPENMP
}
#endif
			// Free memory
			free(R);
			free(G);
			free(B);
		} else {
		// Mono image. Populate the L values into the fits WORD or float data array
#ifdef _OPENMP
#pragma omp parallel
			{
#endif
				if (is_32bit) {
#ifdef _OPENMP
#pragma omp for simd schedule(static,8)
#endif
					for (size_t n = 0; n < npixels; n++) {
						fit->fdata[n] = (float) Lsynth[n];
					}
					if (com.pref.force_16bit) {
						const size_t ndata = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];
						fit_replace_buffer(fit, float_buffer_to_ushort(fit->fdata, ndata), DATA_USHORT);
					}
				} else {
#ifdef _OPENMP
#pragma omp for simd schedule(static,8)
#endif
					for (size_t n = 0; n < npixels; n++) {
						fit->data[n] = roundf_to_WORD(Lsynth[n] * norm);
					}
				}
#ifdef _OPENMP
			}
#endif
		}
	}
	if (H != NULL)
		free(H);
	if (S != NULL)
		free(S);
	if (Hsynth != NULL)
		free(Hsynth);
	if (Ssynth != NULL)
		free(Ssynth);
	free(Lsynth);
	if (buf_needs_freeing) {
		if (is_RGB) {
			for (size_t i = 0; i < 3; i++)
				free(buf[i]);
		} else
			free(buf[RLAYER]);
	}
	update_filter_information(fit, "StarMask", TRUE);
	/* No notify_gfit_data_modified() / gfit_modified_update_gui() here:
	 * generic_image_worker performs both universally when args->fit == gfit. */
	gettimeofday(&t_end, NULL);
	show_time_msg(t_start, t_end, "Execution time");
	gui_iface.set_progress(PROGRESS_RESET, PROGRESS_TEXT_RESET);
	return 0;
}

// Fix up saturated stars only

int reprofile_saturated_stars(fits *fit, star_finder_params *conf_out,
                              const star_finder_params *conf_override) {
	struct timeval t_start, t_end;
	gettimeofday(&t_start, NULL);
	/* unclip always auto-detects (DELEGATED provenance).  Capture: report the
	 * conf used so it can be recorded.  Replay: install the recorded conf
	 * transiently around each channel's detection. */
	if (conf_out)
		*conf_out = conf_override ? *conf_override : com.pref.starfinder_conf;
	char *msg = siril_log_info(_("Star synthesis (desaturating clipped star profiles): processing...\n"));
	msg[strlen(msg) - 1] = '\0';
	gui_iface.set_progress(PROGRESS_RESET, msg);
	gboolean is_RGB = (fit->naxes[2] == 3) ? TRUE : FALSE;
	gboolean is_32bit = TRUE;
	float norm = 1.0f, invnorm = 1.0f;
	if (fit->type == DATA_USHORT) {
		is_32bit = FALSE;
		norm = (float) get_normalized_value(fit);
		invnorm = 1.0f / norm;
	}
	siril_log_debug("norm %f, invnorm %f\n", (float) norm, (float) invnorm);
	int dimx = fit->naxes[0];
	int dimy = fit->naxes[1];
	int count = dimx * dimy;
	float *buf[3];

	buf[RLAYER] = malloc(count * sizeof(float));
	if (is_RGB) {
		buf[GLAYER] = malloc(count * sizeof(float));
		buf[BLAYER] = malloc(count * sizeof(float));
		if (is_32bit) {
			memcpy(buf[RLAYER], fit->fpdata[RLAYER], fit->rx * fit->ry * sizeof(float));
			memcpy(buf[GLAYER], fit->fpdata[GLAYER], fit->rx * fit->ry * sizeof(float));
			memcpy(buf[BLAYER], fit->fpdata[BLAYER], fit->rx * fit->ry * sizeof(float));
		} else {
			for (size_t i = 0; i < count; i++) {
				buf[RLAYER][i] = (float) fit->pdata[RLAYER][i] * invnorm;
				buf[GLAYER][i] = (float) fit->pdata[GLAYER][i] * invnorm;
				buf[BLAYER][i] = (float) fit->pdata[BLAYER][i] * invnorm;
			}
		}
	} else { // mono
		if (is_32bit)
			memcpy(buf[RLAYER], fit->fdata, fit->rx * fit->ry * sizeof(float));
		else {
			buf[RLAYER] = malloc(count * sizeof(float));
			for (size_t i = 0; i < count; i++)
				buf[RLAYER][i] = (float) fit->data[i] * invnorm;
		}
	}

	// Set up starfinder_data structure once, we will reuse it for each channel
	struct starfinder_data sf_data = { 0 };
	sf_data.im.fit = fit;
	sf_data.im.from_seq = NULL;
	sf_data.im.index_in_seq = -1;
	sf_data.max_stars_fitted = MAX_STARS;
	sf_data.selection = (rectangle){0, 0, 0, 0}; // no selection
	sf_data.save_eqcoords = FALSE;
	sf_data.ref_wcs = NULL;
	sf_data.threading = MULTI_THREADED;
	sf_data.update_GUI = FALSE;
	sf_data.process_all_images = FALSE;
	sf_data.already_in_thread = TRUE;
	sf_data.keep_stars = FALSE;

	// Synthesize a PSF for each saturated star in the star array, based on its measured parameters.
	// To fix saturated star profiles we have to do this for each color channel as we can't rely on
	// the hue and saturation within the saturated area, whereas the profiles will be accurate.
	gboolean stopcalled = FALSE;
	for (size_t chan = 0; chan < fit->naxes[2]; chan++) {
		if (stopcalled)
			break;

		psf_star **stars = NULL;
		int nb_stars = 0;

		// Update only the channel-specific fields
		sf_data.layer = chan;
		sf_data.stars = &stars;
		sf_data.nb_stars = &nb_stars;

		// Call the worker function (with the recorded conf installed when replaying)
		int retval = 0;
		if (conf_override)
			WITH_STARFINDER_CONF(conf_override,
					retval = GPOINTER_TO_INT(findstar_worker(&sf_data)));
		else
			retval = GPOINTER_TO_INT(findstar_worker(&sf_data));

		if (retval != 0 || !stars) {
			siril_log_error(_("Star detection failed for channel %u\n"), chan);
			if (stars)
				free_fitted_stars(stars);
			continue; // Skip this channel but continue with others
		}

		int sat_stars = 0;
		siril_log_message(_("Star synthesis: desaturating stars in channel %u...\n"), chan);
		double total = fit->naxes[2] * nb_stars;
		for (size_t n = 0; n < nb_stars; n++) {
			// Check if stop has been pressed
			if (!processing_should_continue())
				stopcalled = TRUE;
			gui_iface.set_progress((double) (n * fit->naxes[2] + chan) / total, NULL);
			if (stars[n]->has_saturated && !stopcalled) {
				float lum = (float) stars[n]->A;
				float bg = (float) stars[n]->B;
				float sat = (float) stars[n]->sat;
				if (lum < 0.0f)
					lum = 0.0f;
				if (!is_32bit) {
					lum *= invnorm;
					bg *= invnorm;
					sat *= invnorm;
				}
				assert(lum >= 0.0f);
				float xoff = (float) stars[n]->xpos - (int) stars[n]->xpos;
				float yoff = (float) stars[n]->ypos - (int) stars[n]->ypos;
				int size = 5.f * max(stars[n]->fwhmx, stars[n]->fwhmy); // This is big enough that it should cover the saturated parts of the star
				if (!(size % 2))
					size++;
				if (size > 1024)
					size = 1024; // Protection against bad star params
				float ratio = stars[n]->fwhmx / stars[n]->fwhmy;
				float angle = (float) stars[n]->angle;

				float *psfL = (float*) calloc(size * size, sizeof(float));
				if (!psfL)
					continue;
				makegaussian(psfL, size, stars[n]->fwhmx, (lum - bg), xoff, yoff, ratio, angle);

				// Replace the part of the profile above the sat threshold
				replace_sat_star_in_buffer(psfL, size, buf[chan],
						(float) stars[n]->xpos, (float) stars[n]->ypos, dimx,
						dimy, sat, bg, 0.f);
				free(psfL);
				sat_stars++;
			}
		}
		free_fitted_stars(stars);
		siril_log_message(_("Star synthesis: %d stars desaturated\n"), sat_stars);
	}

	// Desaturating stars will take their peak brightness over 1.f so we need to rescale the values
	// of all pixels by a factor of (1 / maxbuf) where maxbuf is the maximum subpixel value across all channels
	if (!stopcalled) {
		float bufmax = 1.f;
		for (size_t chan = 0; chan < fit->naxes[2]; chan++)
			for (size_t i = 0; i < count; i++)
				if (buf[chan][i] > bufmax)
					bufmax = buf[chan][i];
		if (bufmax > 1.f) {
			float invbufmax = 1.f / bufmax;
			siril_log_message(_("Remapping output to floating point range 0.0 to 1.0\n"));
			for (size_t chan = 0; chan < fit->naxes[2]; chan++)
				for (size_t i = 0; i < count; i++)
					buf[chan][i] *= invbufmax;
		}
		if (is_32bit) {
			if (is_RGB) {
				memcpy(fit->fpdata[RLAYER], buf[RLAYER], fit->rx * fit->ry * sizeof(float));
				memcpy(fit->fpdata[GLAYER], buf[GLAYER], fit->rx * fit->ry * sizeof(float));
				memcpy(fit->fpdata[BLAYER], buf[BLAYER], fit->rx * fit->ry * sizeof(float));
			}
			else {
				memcpy(fit->fdata, buf[RLAYER], fit->rx * fit->ry * sizeof(float));
			}
		} else {
			for (size_t n = 0; n < count; n++) {
				if (is_RGB) {
					fit->pdata[RLAYER][n] = roundf_to_WORD(buf[RLAYER][n] * norm);
					fit->pdata[GLAYER][n] = roundf_to_WORD(buf[GLAYER][n] * norm);
					fit->pdata[BLAYER][n] = roundf_to_WORD(buf[BLAYER][n] * norm);
				}
				else {
					fit->data[n] = roundf_to_WORD(buf[RLAYER][n] * norm);
				}
			}
		}
	}
	if (is_RGB) {
		for (size_t i = 0; i < 3; i++)
			free(buf[i]);
	} else
		free(buf[RLAYER]);

	/* No notify_gfit_data_modified() / gfit_modified_update_gui() here:
	 * generic_image_worker performs both universally when args->fit == gfit. */
	gettimeofday(&t_end, NULL);
	show_time_msg(t_start, t_end, "Execution time");
	gui_iface.set_progress(PROGRESS_RESET, PROGRESS_TEXT_RESET);
	return 0;
}
