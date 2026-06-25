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
 * Multi-source frame layout (Option B, C4). The MPP "_streamed" engines consume
 * a FrameProvider = std::function<cv::Mat(int idx)> over a single contiguous
 * frame index space. To combine several sequences into one stack we lay their
 * frames end to end in a global index space and dispatch each global index back
 * to its (sequence, local frame). This header is pure index arithmetic plus a
 * thin provider adaptor — no I/O, no Siril types — so it unit-tests on its own.
 */

#ifndef SRC_REGISTRATION_MPP_MULTISOURCE_HPP_
#define SRC_REGISTRATION_MPP_MULTISOURCE_HPP_

#include <algorithm>
#include <functional>
#include <vector>

#include <opencv2/core.hpp>

namespace mpp {

/* The concatenation of several sequences' frames into one global index space.
 * `offsets` holds cumulative frame counts: offsets[k] is the first global index
 * of sequence k, offsets[count()] is the total. The order in which sequences
 * are add()ed is the order they occupy the global space. */
struct MultiSourceLayout {
	std::vector<int> offsets{0};

	/* Append a sequence with `n` frames; returns its sequence index. */
	int add(int n) {
		const int k = static_cast<int>(offsets.size()) - 1;
		offsets.push_back(offsets.back() + n);
		return k;
	}

	int count() const { return static_cast<int>(offsets.size()) - 1; }
	int total() const { return offsets.back(); }
	int frames_of(int seq) const { return offsets[seq + 1] - offsets[seq]; }
	int base_of(int seq) const { return offsets[seq]; }

	/* Map a global index to its (sequence, local frame). Returns false (and
	 * leaves the outputs untouched) for an out-of-range index. */
	bool locate(int global, int *seq, int *local) const {
		if (global < 0 || global >= total())
			return false;
		/* last offset <= global */
		const int k = static_cast<int>(
		    std::upper_bound(offsets.begin(), offsets.end(), global)
		    - offsets.begin()) - 1;
		*seq = k;
		*local = global - offsets[k];
		return true;
	}
};

/* Build a FrameProvider over a layout. `read(seq, local)` returns the frame Mat
 * for sequence index `seq`, local frame `local`; an out-of-range global index
 * yields an empty Mat (which the streamed engines treat as an excluded frame,
 * same as a read failure). The return type is exactly mpp's FrameProvider
 * (std::function<cv::Mat(int)>). */
inline std::function<cv::Mat(int)> make_multisource_provider(
        MultiSourceLayout layout,
        std::function<cv::Mat(int seq, int local)> read) {
	return [layout = std::move(layout), read = std::move(read)](int global) -> cv::Mat {
		int s, l;
		if (!layout.locate(global, &s, &l))
			return cv::Mat();
		return read(s, l);
	};
}

}  // namespace mpp

#endif /* SRC_REGISTRATION_MPP_MULTISOURCE_HPP_ */
