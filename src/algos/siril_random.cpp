/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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

/* SIRIL RANDOM NUMBER GENERATION
 * The random number generators given below use libgsl's implementation
 * of the Mersenne Twister mt19937 algorithm. This is a high quality
 * RNG with good statistical random qualities (it passes the DIEHARD
 * statistical tests) and high performance. The functions provided here
 * should be used in place of rand() which does not produce good quality
 * random numbers.
 */

#include "core/siril.h"
#undef TBYTE
#include "core/proto.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <type_traits>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#include <bcrypt.h>
#else
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#endif

#if defined(__APPLE__)
#include <Security/SecRandom.h>
#endif

static gsl_rng *r;

// Attempts to use the system-specific high quality entropy source as a seed
// If this fails, it falls back to using time(). This is sufficient for Siril.
static uint32_t siril_rng_seed() {
    uint32_t seed = 0;

#if defined(_WIN32) || defined(_WIN64)
    if (BCryptGenRandom(NULL, (PUCHAR)&seed, sizeof(seed), BCRYPT_USE_SYSTEM_PREFERRED_RNG) != 0) {
        siril_debug_print("Failed to generate random seed on Windows.\n");
        goto use_timeofday;
    }
#elif defined(__APPLE__)
    if (SecRandomCopyBytes(kSecRandomDefault, sizeof(seed), (uint8_t*)&seed) != errSecSuccess) {
        siril_debug_print("Failed to generate random seed on macOS.\n");
        goto use_timeofday;
    }
#else
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd == -1) {
        siril_debug_print("Failed to open /dev/urandom");
        goto use_timeofday;
    }
    if (read(fd, &seed, sizeof(seed)) != sizeof(seed)) {
        siril_debug_print("Failed to read from /dev/urandom");
        close(fd);
        goto use_timeofday;
    }
    close(fd);
#endif
    return seed;

use_timeofday:

	return (uint32_t) time(NULL);
}

extern "C" void siril_initialize_rng() {
	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, siril_rng_seed());
}

extern "C" double siril_random_double() {
	return gsl_rng_uniform(r);
}

extern "C" float siril_random_float() {
	return (float) gsl_rng_uniform(r);
}

extern "C" WORD siril_random_WORD() {
	return round_to_WORD(USHRT_MAX_DOUBLE * gsl_rng_uniform(r));
}

extern "C" BYTE siril_random_BYTE() {
	return round_to_BYTE(UCHAR_MAX_DOUBLE * gsl_rng_uniform(r));
}

extern "C" uint32_t siril_random_uint() {
	return (uint32_t) gsl_rng_get(r);
}

extern "C" int32_t siril_random_int() {
	return (int32_t) ((uint32_t) gsl_rng_get(r));
}

template <typename T>
T siril_templated_random() {
    return static_cast<T>(gsl_rng_uniform(r));
}

// Explicit template instantiations
template float siril_templated_random<float>();
template double siril_templated_random<double>();
template uint32_t siril_templated_random<uint32_t>();
template int32_t siril_templated_random<int32_t>();

template <typename T>
T siril_templated_random_max() {
    if (std::is_floating_point<T>::value) {
        return static_cast<T>(1.0);
    } else {
        return static_cast<T>(gsl_rng_max(r));
    }
}

// Explicit template instantiations
template float siril_templated_random_max<float>();
template double siril_templated_random_max<double>();
template uint32_t siril_templated_random_max<uint32_t>();
template int32_t siril_templated_random_max<int32_t>();
