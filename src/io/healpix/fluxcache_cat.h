#ifndef FLUX_CACHE_C_H
#define FLUX_CACHE_C_H

#ifdef __cplusplus
extern "C" {
#endif

void flux_cache_close_all(void);
void flux_cache_purge(int days);
void flux_cache_get_stats(uint32_t *file_count, uint64_t *total_bytes, uint64_t *total_rows);
int* get_cached_l8_healpixels(int *count);

#ifdef __cplusplus
}
#endif

#endif
