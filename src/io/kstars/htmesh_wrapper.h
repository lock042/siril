#pragma once

#ifdef __cplusplus
#include <cstdint>

extern "C" {
	uint32_t get_htm_index_for_coords(double ra, double dec, int levels);
}

#else

#include <stdint.h>
uint32_t get_htm_index_for_coords(double ra, double dec, int levels);
#endif
