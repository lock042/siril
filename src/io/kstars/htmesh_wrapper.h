#pragma once

#ifdef __cplusplus

extern "C" {
	int get_htm_index_for_coords(double ra, double dec, int levels);
	int get_htm_indices_around_target(double ra, double dec, double radius, int levels, int **trixels, int *nb_trixels);
}

#else

int get_htm_index_for_coords(double ra, double dec, int levels);
int get_htm_indices_around_target(double ra, double dec, double radius, int levels, int **trixels, int *nb_trixels);

#endif
