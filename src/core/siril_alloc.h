#ifndef OPTIMAL_ALLOC_H
#define OPTIMAL_ALLOC_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>

// SIMD alignment requirement
// 64 bytes covers AVX-512, which is the most demanding current SIMD instruction set
#define SIMD_ALIGNMENT 64

#ifdef _WIN32

// Threshold for switching to OS-level allocation on Windows
#define LARGE_ALLOC_THRESHOLD (128 * 1024)

// Header stored before each allocation to track allocation method
typedef struct {
	void* original_ptr;    // Original pointer from OS (for freeing)
	size_t size;           // Size of allocation
	unsigned char method;  // 0 = siril_malloc, 1 = OS allocation
	unsigned char padding[7]; // Padding for alignment
} AllocHeader;

#else

#include <stdlib.h>

#endif /* _WIN32 */

/**
	* Check if a pointer is properly aligned for SIMD operations.
	* Useful for debugging and verification on Windows.
	*
	* @param ptr Pointer to check
	* @return 1 if aligned to SIMD_ALIGNMENT, 0 otherwise
	*/
static inline int optimal_is_simd_aligned(void* ptr) {
	return ((uintptr_t)ptr & (SIMD_ALIGNMENT - 1)) == 0;
}

void* siril_malloc(size_t size);
void* siril_calloc(size_t nmemb, size_t size);
void siril_free(void* ptr);
void* siril_realloc(void* ptr, size_t new_size);

#endif /* OPTIMAL_ALLOC_H */
