#pragma once

#ifdef __cplusplus

extern "C" {
	int get_htm_index_for_coords(double ra, double dec, int levels);
	int get_htm_indices_around_target(double ra, double dec, double radius, int levels, int **trixels, int *nb_trixels);
	void get_vertices_for_index(int index, int level, double *ra1, double *dec1, double *ra2, double *dec2, double *ra3, double *dec3);
	int get_htm_indices_around_rectangle(double ra[4], double dec[4], int levels, int **trixels, int *nb_trixels);
}

#else

int get_htm_index_for_coords(double ra, double dec, int levels);
int get_htm_indices_around_target(double ra, double dec, double radius, int levels, int **trixels, int *nb_trixels);
void get_vertices_for_index(int index, int level, double *ra1, double *dec1, double *ra2, double *dec2, double *ra3, double *dec3);
int get_htm_indices_around_rectangle(double ra[4], double dec[4], int levels, int **trixels, int *nb_trixels);

#endif
