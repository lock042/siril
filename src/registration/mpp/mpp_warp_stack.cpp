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
 *
 * This implementation of multipoint registration & stacking is based on
 * PlanetarySystemStacker by Rolf Hempel:
 *     https://github.com/Rolf-Hempel/PlanetarySystemStacker
 */

/*
 * Warp-field stacking engine (MPP_STACK_WARP).
 *
 * The patch-blend engine (mpp_stack.cpp) places each AP's patch rigidly and
 * mosaics the per-AP buffers with a blend window; any disagreement between
 * neighbouring APs is resolved by blending, and what blending can't hide shows
 * up organised on the AP lattice. This engine removes the lattice from the
 * architecture: per frame it interpolates the per-AP shifts into a smooth dense
 * displacement field D and the per-(frame,AP) selection weights into a dense
 * weight field W, warps the frame once with cv::remap through
 *     map = derot_base + global_offset*S - D*S
 * and accumulates acc += (W*mu) (.) warped, wsum += W*mu.
 *
 * It is also the derotation engine: a smooth planetary rotation is a spatially
 * varying displacement that the rigid per-AP patch engine cannot represent in
 * one resample. When a per-frame derotation map is supplied (`derot`), it
 * replaces the identity base of the remap so derotation and AP-residual
 * alignment fold into a single interpolation, and its emission-angle map `mu`
 * down-weights foreshortened near-limb samples while off-disk points (mu = 0)
 * drop out. With no `derot` this is the pure de-warp engine (identity base,
 * mu = 1).
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

#include <opencv2/imgproc.hpp>

#include "registration/mpp.h"
#include "registration/mpp/mpp_stack_priv.hpp"

#include "core/gui_iface.h"
#include "core/processing.h"
#include "core/siril_log.h"

namespace mpp {

/* Down-weight floor for foreshortened (near-limb) samples: the visible disk
 * keeps at least this fraction of weight so the limb is not starved to noise,
 * while off-disk points (mu = 0) are hard-zeroed. */
static constexpr float DEROT_MU_FLOOR = 0.2f;

cv::Mat stack_float_to_uint16(const cv::Mat &buf, int num_layers, int bitdepth) {
	const int C = num_layers;
	const int buf_type = (C == 3) ? CV_32FC3 : CV_32F;
	const double in_max = (bitdepth == 8) ? 255.0 : 65535.0;
	cv::Mat normalised;
	buf.convertTo(normalised, buf_type, 1.0 / in_max);
	cv::min(normalised, 1.0f, normalised);
	cv::max(normalised, 0.0f, normalised);

	const int out_type = (C == 3) ? CV_16UC3 : CV_16U;
	cv::Mat out(normalised.rows, normalised.cols, out_type);
	for (int y = 0; y < normalised.rows; ++y) {
		const float *src = normalised.ptr<float>(y);
		uint16_t *dst = out.ptr<uint16_t>(y);
		const int n = normalised.cols * C;
		for (int x = 0; x < n; ++x) {
			const double v = (double) src[x] * 65535.0 + 0.5;
			const int iv = (int) v;
			dst[x] = (uint16_t) std::clamp(iv, 0, 65535);
		}
	}
	return out;
}

WarpStackResult stack_warp_apply_streamed(const FrameProvider &provider,
                                          int num_frames, int num_layers,
                                          const mpp_aps_t &aps,
                                          const APQualities &apq,
                                          const mpp_shifts_t *shifts,
                                          const std::vector<FrameOffset> &offsets,
                                          const std::vector<double> &frame_brightness,
                                          const std::vector<int> &quality_sorted_idx,
                                          const cv::Vec4i &intersection,
                                          const mpp_config_t &cfg,
                                          const int *included,
                                          int max_threads,
                                          bool provider_thread_safe,
                                          size_t mem_budget_bytes,
                                          const DerotMapProvider *derot) {
	WarpStackResult out;
	const int M = aps.count;
	const int N = num_frames;
	const int C = num_layers;
	const double S = std::max(1.0, cfg.drizzle_scale);
	const int dim_y = intersection[1] - intersection[0];
	const int dim_x = intersection[3] - intersection[2];
	const int DY = (int) std::lround(dim_y * S);
	const int DX = (int) std::lround(dim_x * S);
	if (M <= 0 || N <= 0 || dim_y <= 0 || dim_x <= 0
	 || (int) offsets.size() != N)
		return out;

	/* ---- static geometry: coarse interpolation grid ---- */
	const int step = mpp_cfg_step_size(&cfg);
	const int cs = std::max(8, step / 3);
	const int Gy = std::max(2, (int) std::lround((double) dim_y / cs));
	const int Gx = std::max(2, (int) std::lround((double) dim_x / cs));
	auto node_pos = [](int g, int dim, int G) {
		return ((double) g + 0.5) * (double) dim / (double) G - 0.5;
	};

	/* Per-AP scatter weight = the AP's own Hann blend window over its patch
	 * bounds (frame-edge extensions included): at an AP centre every other AP
	 * is outside its own window, so the interpolated displacement there is
	 * exactly that AP's measured shift, and overlap zones blend displacements
	 * with the patch engine's spatial profile. */
	std::vector<std::vector<float>> ap_wy(M), ap_wx(M);
	for (int a = 0; a < M; ++a) {
		const auto &r = aps.records[a];
		ap_wy[a] = stack_one_dim_weight(r.patch_y_low, r.patch_y_high, r.y,
		                                r.patch_y_low == 0,
		                                r.patch_y_high == dim_y);
		ap_wx[a] = stack_one_dim_weight(r.patch_x_low, r.patch_x_high, r.x,
		                                r.patch_x_low == 0,
		                                r.patch_x_high == dim_x);
	}

	struct NodeAP { int ap; float lambda; };
	std::vector<std::vector<NodeAP>> node_aps((size_t) Gy * Gx);
	for (int gy = 0; gy < Gy; ++gy) {
		const int ny = (int) std::lround(node_pos(gy, dim_y, Gy));
		for (int gx = 0; gx < Gx; ++gx) {
			const int nx = (int) std::lround(node_pos(gx, dim_x, Gx));
			auto &list = node_aps[(size_t) gy * Gx + gx];
			double lsum = 0.0;
			for (int a = 0; a < M; ++a) {
				const auto &r = aps.records[a];
				if (ny < r.patch_y_low || ny >= r.patch_y_high
				 || nx < r.patch_x_low || nx >= r.patch_x_high)
					continue;
				const double l = (double) ap_wy[a][ny - r.patch_y_low]
				               * (double) ap_wx[a][nx - r.patch_x_low];
				if (l <= 0.0) continue;
				list.push_back({a, (float) l});
				lsum += l;
			}
			if (lsum > 1.0) {
				const float inv = (float) (1.0 / lsum);
				for (auto &e : list) e.lambda *= inv;
			}
		}
	}

	/* AP neighbour lists for the robust shift clamp. */
	const int nbr_radius = (step * 8) / 5;   /* 1.6 × grid step */
	std::vector<std::vector<int>> ap_nbrs(M);
	for (int a = 0; a < M; ++a)
		for (int b = a + 1; b < M; ++b)
			if (std::abs(aps.records[a].y - aps.records[b].y) <= nbr_radius
			 && std::abs(aps.records[a].x - aps.records[b].x) <= nbr_radius) {
				ap_nbrs[a].push_back(b);
				ap_nbrs[b].push_back(a);
			}
	const double clamp_thr = std::max(2.0,
	    0.5 * (double) cfg.alignment_points_search_width);

	/* Static coverage → does the canvas need the background blend? */
	bool holes = false;
	{
		cv::Mat covc(Gy, Gx, CV_32F);
		for (int gy = 0; gy < Gy; ++gy) {
			float *row = covc.ptr<float>(gy);
			for (int gx = 0; gx < Gx; ++gx) {
				double s = 0.0;
				for (const auto &e : node_aps[(size_t) gy * Gx + gx])
					s += e.lambda;
				row[gx] = (float) s;
			}
		}
		cv::Mat cov, hm;
		cv::resize(covc, cov, cv::Size(DX, DY), 0, 0, cv::INTER_CUBIC);
		cv::compare(cov, 0.05, hm, cv::CMP_LT);
		holes = cv::countNonZero(hm) > 0;
	}

	std::set<int> top_for_bg;
	for (int i = 0; i < apq.stack_size && i < (int) quality_sorted_idx.size(); ++i)
		top_for_bg.insert(quality_sorted_idx[i]);

	double median_brightness = 0.0;
	if (cfg.frames_normalization) {
		std::vector<double> sorted_b(frame_brightness);
		std::sort(sorted_b.begin(), sorted_b.end());
		const size_t n = sorted_b.size();
		median_brightness = (n % 2 == 1)
		    ? sorted_b[n / 2]
		    : 0.5 * (sorted_b[n / 2 - 1] + sorted_b[n / 2]);
	}

	/* Shared read-only identity maps (used as the remap base when there is no
	 * derotation, and to recover the residual map otherwise). */
	cv::Mat rowIdx(DY, DX, CV_32F), colIdx(DY, DX, CV_32F);
	for (int y = 0; y < DY; ++y) {
		float *ry = rowIdx.ptr<float>(y);
		float *cx = colIdx.ptr<float>(y);
		for (int x = 0; x < DX; ++x) { ry[x] = (float) y; cx[x] = (float) x; }
	}

	std::vector<int> work;
	work.reserve(N);
	for (int f = 0; f < N; ++f)
		if (!included || included[f]) work.push_back(f);
	const int Wn = (int) work.size();
	const int n_included = Wn > 0 ? Wn : 1;

	const int nt = max_threads > 0 ? max_threads : 1;
	const bool par = provider_thread_safe && nt > 1;
	int n_threads = par ? nt : 1;
	if (par && mem_budget_bytes > 0) {
		const size_t canvas = (size_t) DY * DX * sizeof(float);
		const size_t per_thread = canvas * ((size_t) (holes ? 3 : 2) * C + 10);
		int fit = (int) (mem_budget_bytes / std::max<size_t>(1, per_thread));
		if (fit < 1) fit = 1;
		if (fit < n_threads) n_threads = fit;
	}

	std::vector<cv::Mat> acc(C);
	std::vector<cv::Mat> bg(C);
	cv::Mat wsum;

	gint cur_nb = 0;
	gint cancelled = 0;
	gint oom = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads) if (n_threads > 1)
#endif
	{
		bool local_ok = true;
		std::vector<cv::Mat> lacc(C), lbg(C);
		cv::Mat lwsum;
		try {
			for (int c = 0; c < C; ++c) {
				lacc[c] = cv::Mat::zeros(DY, DX, CV_32F);
				if (holes) lbg[c] = cv::Mat::zeros(DY, DX, CV_32F);
			}
			lwsum = cv::Mat::zeros(DY, DX, CV_32F);
		} catch (const std::exception &) {
			local_ok = false;
			g_atomic_int_set(&oom, 1);
		}
		std::vector<float> wA(M, 0.0f);
		std::vector<float> syA(M, 0.0f), sxA(M, 0.0f);
		std::vector<uint8_t> okA(M, 0);
		cv::Mat Wc(Gy, Gx, CV_32F), Dyc(Gy, Gx, CV_32F), Dxc(Gy, Gx, CV_32F);

#ifdef _OPENMP
#pragma omp for schedule(dynamic) nowait
#endif
		for (int wi = 0; wi < Wn; ++wi) {
			if (g_atomic_int_get(&cancelled) || g_atomic_int_get(&oom)) continue;
			if (!processing_should_continue()) { g_atomic_int_set(&cancelled, 1); continue; }
			const int f = work[wi];
			try {

			const cv::Mat frame_raw = provider(f);
			if (frame_raw.empty()) {
				gui_iface.set_progress(0.5 + 0.5 * (double) g_atomic_int_add(&cur_nb, 1)
				                       / (double) n_included, NULL);
				continue;
			}
			cv::Mat frame_f32;
			if (cfg.frames_normalization) {
				const double sc = median_brightness / (frame_brightness[f] + 1e-7);
				frame_raw.convertTo(frame_f32, CV_32F, sc);
			} else {
				frame_raw.convertTo(frame_f32, CV_32F);
			}
			cv::Mat frame_drizzled;
			if (S != 1.0)
				cv::resize(frame_f32, frame_drizzled,
				           cv::Size((int) std::lround(frame_f32.cols * S),
				                    (int) std::lround(frame_f32.rows * S)),
				           0, 0, cv::INTER_LINEAR);
			else
				frame_drizzled = frame_f32;
			const int Hd = frame_drizzled.rows;
			const int Wd_f = frame_drizzled.cols;

			/* Per-frame AP weights and shifts. */
			const auto &used = apq.used_alignment_points[f];
			for (const auto &u : used) {
				wA[u.ap] = u.weight;
				if (shifts) {
					const size_t soff = (size_t) (f * M + u.ap) * 2;
					syA[u.ap] = (float) shifts->shifts[soff + 0];
					sxA[u.ap] = (float) shifts->shifts[soff + 1];
					okA[u.ap] = shifts->success[(size_t) f * M + u.ap];
				} else {
					syA[u.ap] = 0.0f;
					sxA[u.ap] = 0.0f;
					okA[u.ap] = 1;
				}
			}

			/* Robust clamp: snap an AP deviating from its neighbours' median by
			 * more than the threshold onto that median, computed against the
			 * original values and committed afterwards (order-independent). */
			std::vector<std::pair<int, cv::Vec2f>> repl;
			std::vector<float> vy, vx;
			for (const auto &u : used) {
				const int a = u.ap;
				if (!okA[a]) continue;
				vy.clear(); vx.clear();
				for (int b : ap_nbrs[a]) {
					if (wA[b] <= 0.0f || !okA[b]) continue;
					vy.push_back(syA[b]);
					vx.push_back(sxA[b]);
				}
				if (vy.size() < 3) continue;
				std::nth_element(vy.begin(), vy.begin() + vy.size() / 2, vy.end());
				std::nth_element(vx.begin(), vx.begin() + vx.size() / 2, vx.end());
				const float my = vy[vy.size() / 2];
				const float mx = vx[vx.size() / 2];
				if (std::abs(syA[a] - my) > clamp_thr
				 || std::abs(sxA[a] - mx) > clamp_thr)
					repl.emplace_back(a, cv::Vec2f(my, mx));
			}
			for (const auto &r : repl) {
				syA[r.first] = r.second[0];
				sxA[r.first] = r.second[1];
			}

			/* Coarse node fields. */
			for (int gy = 0; gy < Gy; ++gy) {
				float *wrow = Wc.ptr<float>(gy);
				float *yrow = Dyc.ptr<float>(gy);
				float *xrow = Dxc.ptr<float>(gy);
				for (int gx = 0; gx < Gx; ++gx) {
					double wsum_n = 0.0, den = 0.0, ny = 0.0, nx = 0.0;
					for (const auto &e : node_aps[(size_t) gy * Gx + gx]) {
						const double lw = (double) e.lambda * wA[e.ap];
						if (lw <= 0.0) continue;
						wsum_n += lw;
						if (!okA[e.ap]) continue;
						den += lw;
						ny += lw * syA[e.ap];
						nx += lw * sxA[e.ap];
					}
					wrow[gx] = (float) wsum_n;
					yrow[gx] = den > 1e-12 ? (float) (ny / den) : 0.0f;
					xrow[gx] = den > 1e-12 ? (float) (nx / den) : 0.0f;
				}
			}

			/* Reset scratch for the next frame. */
			for (const auto &u : used) {
				wA[u.ap] = 0.0f;
				syA[u.ap] = sxA[u.ap] = 0.0f;
				okA[u.ap] = 0;
			}

			/* Dense weight + displacement fields. */
			cv::Mat Wd, Dyd, Dxd;
			cv::resize(Wc, Wd, cv::Size(DX, DY), 0, 0, cv::INTER_CUBIC);
			cv::resize(Dyc, Dyd, cv::Size(DX, DY), 0, 0, cv::INTER_CUBIC);
			cv::resize(Dxc, Dxd, cv::Size(DX, DY), 0, 0, cv::INTER_CUBIC);
			cv::max(Wd, 0.0f, Wd);   /* cubic overshoot */

			/* Optional per-frame derotation base. dmx/dmy give, for each
			 * epoch-canvas pixel, the source position in the MPP-aligned frame
			 * canvas of the same surface point; dmu is the emission cosine
			 * (0 off-disk). With no provider the base is the identity. */
			cv::Mat dmx, dmy, dmu;
			const bool have_derot = derot && (*derot)(f, DY, DX, dmx, dmy, dmu);

			cv::Mat mapy, mapx;
			cv::add(have_derot ? dmy : rowIdx, cv::Scalar(offsets[f].dy * S), mapy);
			cv::add(have_derot ? dmx : colIdx, cv::Scalar(offsets[f].dx * S), mapx);
			cv::scaleAdd(Dyd, -S, mapy, mapy);
			cv::scaleAdd(Dxd, -S, mapx, mapx);

			/* Zero the weight where the sample leaves the source frame. */
			cv::Mat vy8, vx8, v8, vf;
			cv::inRange(mapy, 0.0, (double) (Hd - 1), vy8);
			cv::inRange(mapx, 0.0, (double) (Wd_f - 1), vx8);
			cv::bitwise_and(vy8, vx8, v8);
			v8.convertTo(vf, CV_32F, 1.0 / 255.0);
			cv::multiply(Wd, vf, Wd);

			/* Fold the derotation foreshortening weight: max(mu, floor) on the
			 * visible disk, hard zero off-disk. This also prevents the -1
			 * off-disk map sentinel (now offset into range) from contributing. */
			if (have_derot) {
				cv::Mat wmu, onmask, onf;
				cv::max(dmu, DEROT_MU_FLOOR, wmu);
				cv::compare(dmu, 0.0f, onmask, cv::CMP_GT);
				onmask.convertTo(onf, CV_32F, 1.0 / 255.0);
				cv::multiply(wmu, onf, wmu);
				cv::multiply(Wd, wmu, Wd);
			}

			cv::Mat warped;
			cv::remap(frame_drizzled, warped, mapx, mapy,
			          cv::INTER_LANCZOS4, cv::BORDER_REPLICATE);

			cv::Mat tmp;
			if (C == 1) {
				cv::multiply(warped, Wd, tmp);
				lacc[0] += tmp;
				if (holes && top_for_bg.find(f) != top_for_bg.end())
					lbg[0] += warped;
			} else {
				std::vector<cv::Mat> chans;
				cv::split(warped, chans);
				const bool in_bg = holes && top_for_bg.find(f) != top_for_bg.end();
				for (int c = 0; c < C; ++c) {
					cv::multiply(chans[c], Wd, tmp);
					lacc[c] += tmp;
					if (in_bg) lbg[c] += chans[c];
				}
			}
			lwsum += Wd;

			gui_iface.set_progress(0.5 + 0.5 * (double) g_atomic_int_add(&cur_nb, 1)
			                       / (double) n_included, NULL);

			} catch (const std::exception &) {
				g_atomic_int_set(&oom, 1);
				continue;
			}
		}

#ifdef _OPENMP
#pragma omp critical(mpp_warp_stack_reduce)
#endif
		{
			if (local_ok && !g_atomic_int_get(&oom)) {
				try {
					if (wsum.empty()) {
						wsum = lwsum;
						for (int c = 0; c < C; ++c) {
							acc[c] = lacc[c];
							if (holes) bg[c] = lbg[c];
						}
					} else {
						wsum += lwsum;
						for (int c = 0; c < C; ++c) {
							acc[c] += lacc[c];
							if (holes) bg[c] += lbg[c];
						}
					}
				} catch (const std::exception &) {
					g_atomic_int_set(&oom, 1);
				}
			}
		}
	}

	if (g_atomic_int_get(&oom)) { out.oom = true; return out; }
	if (cancelled || wsum.empty()) { out.cancelled = true; return out; }

	/* Compose acc / wsum, then the PSS soft background blend where AP coverage
	 * thins (fg ramp = wsum / (blend_t × stack_size)). */
	cv::Mat denom;
	cv::max(wsum, 1e-30f, denom);
	std::vector<cv::Mat> mean_ch(C);
	for (int c = 0; c < C; ++c)
		cv::divide(acc[c], denom, mean_ch[c]);
	if (holes && apq.stack_size > 0) {
		cv::Mat fg;
		const double blend_denom = cfg.stack_frames_background_blend_threshold
		                         * (double) apq.stack_size;
		wsum.convertTo(fg, CV_32F, 1.0 / blend_denom);
		cv::min(fg, 1.0f, fg);
		cv::max(fg, 0.0f, fg);
		for (int c = 0; c < C; ++c) {
			cv::Mat bgm, diff;
			bg[c].convertTo(bgm, CV_32F, 1.0 / (double) apq.stack_size);
			cv::subtract(mean_ch[c], bgm, diff);
			cv::multiply(diff, fg, diff);
			cv::add(bgm, diff, mean_ch[c]);
		}
	}

	cv::Mat buf;
	if (C == 1)
		buf = mean_ch[0];
	else
		cv::merge(mean_ch, buf);
	out.image = stack_float_to_uint16(buf, C, cfg.bitdepth);
	return out;
}

}  // namespace mpp
