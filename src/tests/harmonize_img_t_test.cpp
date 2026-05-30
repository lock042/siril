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

/*
 * Unit tests for the shared img_t image type and the harmonization ops
 * (branch harmonize_img_t). Built as a Criterion C++ test.
 */

#include <criterion/criterion.h>
#include "filters/deconvolution/image.hpp"
#include "filters/deconvolution/image_expr.hpp"
#include "filters/deconvolution/image_ops.hpp"
#include "filters/deconvolution/image_dft.hpp"

// The Siril core globals are not linked into the unit-test binary (see the
// other tests in this directory), so provide local definitions. img_t reads
// com.max_thread for OpenMP work-splitting and references sequence_is_running.
cominfo com;
int sequence_is_running = 0;

// img_t uses com.max_thread for its OpenMP work splitting. In a unit-test
// process the global is zero-initialised and Siril is never started, so pin it
// to a sane value before any img_t op runs.
static void setup_threads(void) {
	com.max_thread = 1;
}

TestSuite(harmonize_img_t, .init = setup_threads);

Test(harmonize_img_t, construct_and_planar_index) {
	img_t<float> a(3, 2, 1); // w=3, h=2, 1 channel
	a.set_value(0.f);
	a(0, 0) = 1.f;
	a(2, 1) = 5.f;
	// planar single-channel index: data[x + w*y]
	cr_assert_float_eq(a.data[0], 1.f, 1e-6, "(0,0) maps to data[0]");
	cr_assert_float_eq(a.data[5], 5.f, 1e-6, "(2,1) maps to data[2 + 3*1] = data[5]");
	cr_assert_eq(a.w, 3);
	cr_assert_eq(a.h, 2);
	cr_assert_eq(a.d, 1);
}

Test(harmonize_img_t, planar_multichannel_layout) {
	img_t<float> a(2, 2, 3);
	a.set_value(0.f);
	// channel 1, pixel (1,0) -> data[x + w*(y + h*dd)] = data[1 + 2*(0 + 2*1)] = data[5]
	a(1, 0, 1) = 7.f;
	cr_assert_float_eq(a.data[5], 7.f, 1e-6, "planar channel offset");
}

Test(harmonize_img_t, sum_reduction) {
	img_t<float> a(4, 4, 1);
	a.set_value(2.f);
	cr_assert_float_eq(a.sum(), 32.f, 1e-4, "sum of 16 cells of 2.0");
}

/* ---- imgops: symmetric reflection ---- */

Test(imgops, symmetric_coordinate_half_sample) {
	// in range: identity
	cr_assert_eq(imgops::symmetric_coordinate(0, 5), 0);
	cr_assert_eq(imgops::symmetric_coordinate(4, 5), 4);
	// left of edge: -1 -> 0, -2 -> 1 (edge not duplicated)
	cr_assert_eq(imgops::symmetric_coordinate(-1, 5), 0);
	cr_assert_eq(imgops::symmetric_coordinate(-2, 5), 1);
	// right of edge: size -> size-1, size+1 -> size-2
	cr_assert_eq(imgops::symmetric_coordinate(5, 5), 4);
	cr_assert_eq(imgops::symmetric_coordinate(6, 5), 3);
}

Test(imgops, reflect_whole_sample_matches_deconv_inline) {
	// whole-sample reflection used by process_in_slices: -1 -> 1, size -> size-2
	cr_assert_eq(imgops::reflect_whole_sample(0, 5), 0);
	cr_assert_eq(imgops::reflect_whole_sample(4, 5), 4);
	cr_assert_eq(imgops::reflect_whole_sample(-1, 5), 1);
	cr_assert_eq(imgops::reflect_whole_sample(-2, 5), 2);
	cr_assert_eq(imgops::reflect_whole_sample(5, 5), 3); // 2*5 - 5 - 2
	cr_assert_eq(imgops::reflect_whole_sample(6, 5), 2); // 2*5 - 6 - 2
}

/* ---- imgops: monochrome ---- */

Test(imgops, is_monochrome) {
	img_t<float> mono(3, 3, 1);
	mono.set_value(0.5f);
	cr_assert(imgops::is_monochrome(mono), "single channel is monochrome");

	img_t<float> grey3(2, 2, 3);
	grey3.set_value(0.3f); // all channels equal everywhere
	cr_assert(imgops::is_monochrome(grey3), "equal channels => monochrome");

	img_t<float> colour(2, 2, 3);
	colour.set_value(0.3f);
	colour(1, 1, 2) = 0.9f; // one channel differs at one pixel
	cr_assert_not(imgops::is_monochrome(colour), "differing channel => not monochrome");
}

Test(imgops, make_monochrome_averages_channels) {
	img_t<float> c(1, 1, 3);
	c(0, 0, 0) = 0.0f;
	c(0, 0, 1) = 0.6f;
	c(0, 0, 2) = 0.9f;
	img_t<float> m = imgops::make_monochrome(c);
	cr_assert_eq(m.d, 1);
	cr_assert_float_eq(m(0, 0), 0.5f, 1e-6, "mean of 0,0.6,0.9 = 0.5");
}

/* ---- imgops: symmetric padding ---- */

Test(imgops, pad_symmetric_reflects_and_roundtrips) {
	img_t<float> in(3, 2, 1);
	// fill with distinct values: value = x + 10*y
	for (int y = 0; y < 2; ++y)
		for (int x = 0; x < 3; ++x)
			in(x, y) = (float)(x + 10 * y);

	const int b = 2;
	img_t<float> p = imgops::pad_symmetric(in, b);
	cr_assert_eq(p.w, 3 + 2 * b);
	cr_assert_eq(p.h, 2 + 2 * b);
	// interior is the original
	for (int y = 0; y < 2; ++y)
		for (int x = 0; x < 3; ++x)
			cr_assert_float_eq(p(x + b, y + b), in(x, y), 1e-6, "interior preserved");
	// left border reflects column 0,1 -> at padded x=b-1 we expect in col 0
	cr_assert_float_eq(p(b - 1, b), in(0, 0), 1e-6, "x=-1 reflects to col 0");
	cr_assert_float_eq(p(b - 2, b), in(1, 0), 1e-6, "x=-2 reflects to col 1");

	// unpad recovers the original exactly
	img_t<float> u = imgops::unpad(p, b);
	cr_assert_eq(u.w, 3);
	cr_assert_eq(u.h, 2);
	for (int y = 0; y < 2; ++y)
		for (int x = 0; x < 3; ++x)
			cr_assert_float_eq(u(x, y), in(x, y), 1e-6, "pad->unpad identity");
}

/* ---- imgops: tiling ---- */

Test(imgops, compute_tiling_cases) {
	imgops::tiling_t t1 = imgops::compute_tiling(8, 8, 4); // square, 4 tiles
	cr_assert_eq(t1.ny, 2);
	cr_assert_eq(t1.nx, 2);

	imgops::tiling_t t2 = imgops::compute_tiling(8, 8, 1);
	cr_assert_eq(t2.ny, 1);
	cr_assert_eq(t2.nx, 1);
}

Test(imgops, split_merge_roundtrip_uniform_weight) {
	img_t<float> in(6, 4, 1);
	for (int y = 0; y < 4; ++y)
		for (int x = 0; x < 6; ++x)
			in(x, y) = (float)(x + 10 * y);

	const int pad = 1;
	imgops::tiling_t tiling = imgops::compute_tiling(in.h, in.w, 4);
	std::vector<img_t<float>> tiles = imgops::split_tiles(in, pad, pad, tiling);

	// build {value, weight=1} pairs from the split tiles
	std::vector<std::pair<img_t<float>, img_t<float>>> merge_in;
	for (img_t<float>& t : tiles) {
		img_t<float> w(t.w, t.h, 1);
		w.set_value(1.f);
		merge_in.emplace_back(t, std::move(w)); // copy-construct value from tile
	}

	img_t<float> out = imgops::merge_tiles(merge_in, in.w, in.h, pad, pad, tiling);
	cr_assert_eq(out.w, 6);
	cr_assert_eq(out.h, 4);
	for (int y = 0; y < 4; ++y)
		for (int x = 0; x < 6; ++x)
			cr_assert_float_eq(out(x, y), in(x, y), 1e-5,
			                   "split->merge reconstructs original");
}

Test(imgops, merge_tiles_weighted_average) {
	// single tile covering the whole image: result = value / weight
	std::vector<std::pair<img_t<float>, img_t<float>>> merge_in;
	img_t<float> val(2, 2, 1);
	val.set_value(10.f);
	img_t<float> wgt(2, 2, 1);
	wgt.set_value(2.f);
	merge_in.emplace_back(std::move(val), std::move(wgt));

	imgops::tiling_t tiling{1, 1};
	img_t<float> out = imgops::merge_tiles(merge_in, 2, 2, 0, 0, tiling);
	for (int y = 0; y < 2; ++y)
		for (int x = 0; x < 2; ++x)
			cr_assert_float_eq(out(x, y), 5.f, 1e-6, "10/2 = 5");
}

/* ---- img_expr_t: axis reduce / broadcast ---- */

Test(img_expr, reduce_axis_d_matches_reduce_d) {
	img_t<float> a(2, 2, 3);
	for (int c = 0; c < 3; ++c)
		for (int y = 0; y < 2; ++y)
			for (int x = 0; x < 2; ++x)
				a(x, y, c) = (float)(c + 1); // channel values 1,2,3

	img_t<float> r;
	r.resize(2, 2, 1);
	r.map(reduce_axis(a, AXIS_D,
	                  std::function<float(float, float)>([](float u, float v){ return u + v; }), 0.f));
	for (int y = 0; y < 2; ++y)
		for (int x = 0; x < 2; ++x)
			cr_assert_float_eq(r(x, y), 6.f, 1e-6, "1+2+3 = 6 along channels");
}

Test(img_expr, reduce_axis_x_row_sums) {
	img_t<float> a(3, 2, 1);
	// row 0: 1,2,3 ; row 1: 10,20,30
	a(0,0)=1; a(1,0)=2; a(2,0)=3;
	a(0,1)=10; a(1,1)=20; a(2,1)=30;
	img_t<float> r;
	r.resize(1, 2, 1);
	r.map(reduce_axis(a, AXIS_X,
	                  std::function<float(float, float)>([](float u, float v){ return u + v; }), 0.f));
	cr_assert_eq(r.w, 1);
	cr_assert_eq(r.h, 2);
	cr_assert_float_eq(r(0, 0), 6.f, 1e-6, "row 0 sum");
	cr_assert_float_eq(r(0, 1), 60.f, 1e-6, "row 1 sum");
}

Test(img_expr, mean_axis_x) {
	img_t<float> a(4, 1, 1);
	a(0,0)=2; a(1,0)=4; a(2,0)=6; a(3,0)=8;
	img_t<float> m = imgops::mean_axis(a, AXIS_X);
	cr_assert_eq(m.w, 1);
	cr_assert_float_eq(m(0, 0), 5.f, 1e-6, "(2+4+6+8)/4 = 5");
}

Test(img_expr, broadcast_axis_repeats) {
	img_t<float> col(1, 2, 1);
	col(0, 0) = 7.f;
	col(0, 1) = 9.f;
	img_t<float> b;
	b.resize(3, 2, 1);
	b.map(broadcast_axis(col, AXIS_X, 3));
	for (int x = 0; x < 3; ++x) {
		cr_assert_float_eq(b(x, 0), 7.f, 1e-6, "row 0 broadcast");
		cr_assert_float_eq(b(x, 1), 9.f, 1e-6, "row 1 broadcast");
	}
}

Test(img_expr, center_data_via_reduce_broadcast) {
	// NL-Bayes style: group(nSimP=3 along x, pixels=2 along y); subtract the
	// per-pixel mean across patches (mean along x) broadcast back over x.
	img_t<float> group(3, 2, 1);
	group(0,0)=1; group(1,0)=2; group(2,0)=3;   // row 0 mean = 2
	group(0,1)=10; group(1,1)=30; group(2,1)=50; // row 1 mean = 30
	img_t<float> baricenter = imgops::mean_axis(group, AXIS_X);
	img_t<float> centered;
	centered.resize(3, 2, 1);
	centered.map(group - broadcast_axis(baricenter, AXIS_X, group.w));
	// each row of centered should now sum to ~0
	cr_assert_float_eq(centered(0,0) + centered(1,0) + centered(2,0), 0.f, 1e-5, "row0 centered");
	cr_assert_float_eq(centered(0,1) + centered(1,1) + centered(2,1), 0.f, 1e-5, "row1 centered");
	cr_assert_float_eq(centered(1,0), 0.f, 1e-6, "2 - 2 = 0");
	cr_assert_float_eq(centered(2,1), 20.f, 1e-5, "50 - 30 = 20");
}

/* ---- imgops: patch DFT ---- */

Test(imgops, dft_patch_dc_component) {
	imgops::dft_patch p(8, 8, 1);
	p.space().set_value(3.f); // constant -> DC = sum, all other bins ~0
	p.to_freq();
	// DC term = sum of all samples = 3 * 64
	cr_assert_float_eq(p.freq(0, 0).real(), 3.f * 64.f, 1e-2, "DC = sum");
	cr_assert_float_eq(p.freq(0, 0).imag(), 0.f, 1e-3, "DC imag ~ 0");
	cr_assert_eq(p.fw(), 5); // 8/2 + 1
}

Test(imgops, dft_patch_roundtrip) {
	imgops::dft_patch p(16, 12, 2);
	// fill with a deterministic pattern
	for (int c = 0; c < 2; ++c)
		for (int y = 0; y < 12; ++y)
			for (int x = 0; x < 16; ++x)
				p.space(x, y, c) = std::sin(0.3f * x + 0.2f * y + c) + 0.5f * c;

	img_t<float> original = p.space(); // copy

	p.to_freq();
	p.to_space();

	for (int c = 0; c < 2; ++c)
		for (int y = 0; y < 12; ++y)
			for (int x = 0; x < 16; ++x)
				cr_assert_float_eq(p.space(x, y, c), original(x, y, c), 1e-4,
				                   "forward+inverse is identity");
}
