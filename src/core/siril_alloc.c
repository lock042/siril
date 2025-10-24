#include "core/siril_alloc.h"

/**
* Allocate memory optimally for the given size.
* On POSIX systems, uses malloc (which automatically switches to mmap for large allocations).
* On Windows, small allocations use malloc/free, large allocations use VirtualAlloc
* to avoid heap fragmentation, with SIMD-aligned pointers (64-byte).
*
* @param size Number of bytes to allocate
* @return Pointer to allocated memory, or NULL on failure
*/
void* siril_malloc(size_t size) {
#ifdef _WIN32
	if (size >= LARGE_ALLOC_THRESHOLD) {
		// Large allocation - use VirtualAlloc with SIMD alignment
		// VirtualAlloc returns 64KB-aligned memory, which is already SIMD-aligned
		// We need extra space for header + potential alignment padding
		size_t page_size = 4096;
		size_t total_size = size + sizeof(AllocHeader) + SIMD_ALIGNMENT;
		size_t aligned_size = (total_size + page_size - 1) & ~(page_size - 1);

		LPVOID os_ptr = VirtualAlloc(NULL,                    // Let OS choose address
									aligned_size,            // Size rounded to page
									MEM_COMMIT | MEM_RESERVE, // Reserve and commit
									PAGE_READWRITE);         // Read/write access

		if (os_ptr == NULL) {
			return NULL;  // Allocation failed
		}

		// Calculate where the user pointer should be (after header, aligned)
		// Start with position after header
		uintptr_t after_header = (uintptr_t)os_ptr + sizeof(AllocHeader);

		// Round up to next SIMD alignment boundary
		uintptr_t user_ptr = (after_header + SIMD_ALIGNMENT - 1) & ~(SIMD_ALIGNMENT - 1);

		// Store header just before the aligned user pointer
		AllocHeader* header = (AllocHeader*)(user_ptr - sizeof(AllocHeader));
		header->original_ptr = os_ptr;  // Save original pointer for VirtualFree
		header->size = aligned_size;    // Store actual allocated size
		header->method = 1;              // Mark as OS allocation

		// Return SIMD-aligned pointer
		return (void*)user_ptr;
	} else {
		// Small allocation - use malloc (fast, low overhead)
		// No SIMD alignment guarantee for small allocations
		size_t total_size = size + sizeof(AllocHeader);
		void* os_ptr = malloc(total_size);

		if (os_ptr == NULL) {
			return NULL;  // Allocation failed
		}

		// For small allocations, header is at the start, user data follows
		AllocHeader* header = (AllocHeader*)os_ptr;
		header->original_ptr = os_ptr;  // Same as os_ptr for malloc
		header->size = total_size;      // Store actual allocated size
		header->method = 0;              // Mark as malloc allocation

		// Return pointer after header (not aligned)
		return (void*)((char*)os_ptr + sizeof(AllocHeader));
	}
#else
	// POSIX: glibc automatically uses mmap for large allocations (≥128KB)
	// and provides page-aligned memory, which is already SIMD-aligned
	return malloc(size);
#endif
}

/**
* Allocate and zero-initialize memory optimally for the given size.
* Same behavior as optimal_alloc, but memory is zeroed.
*
* @param nmemb Number of elements
* @param size Size of each element in bytes
* @return Pointer to allocated and zeroed memory, or NULL on failure
*/
void* siril_calloc(size_t nmemb, size_t size) {
#ifdef _WIN32
	// Check for overflow
	if (nmemb != 0 && size > SIZE_MAX / nmemb) {
		return NULL;  // Would overflow
	}

	size_t total = nmemb * size;
	void* ptr = siril_malloc(total);

	if (ptr != NULL) {
		// Zero the memory
		memset(ptr, 0, total);
	}

	return ptr;
#else
	// POSIX: use standard calloc
	return calloc(nmemb, size);
#endif
}

/**
* Free memory allocated with optimal_alloc or optimal_calloc.
* On POSIX systems, uses free.
* On Windows, automatically uses the correct deallocation method.
*
* @param ptr Pointer returned by optimal_alloc/optimal_calloc (can be NULL)
*/
void siril_free(void* ptr) {
	if (ptr == NULL) {
		return;  // NULL pointer is safe to free
	}

#ifdef _WIN32
	// Get the header that's stored just before the user's pointer
	AllocHeader* header = (AllocHeader*)((char*)ptr - sizeof(AllocHeader));

	if (header->method == 1) {
		// OS allocation - use VirtualFree
		VirtualFree(header->original_ptr,  // Original pointer from VirtualAlloc
				0,                      // Size must be 0 for MEM_RELEASE
				MEM_RELEASE);           // Release both reserved and committed
	} else {
		// malloc allocation - use free
		free(header->original_ptr);
	}
#else
	// POSIX: use standard free
	free(ptr);
#endif
}

/**
* Reallocate memory allocated with optimal_alloc.
* On POSIX systems, uses realloc.
* On Windows, handles transitions between allocation methods if size crosses threshold,
* and for large allocations, the returned pointer is SIMD-aligned.
*
* @param ptr Pointer returned by optimal_alloc (can be NULL)
* @param new_size New size in bytes
* @return Pointer to reallocated memory, or NULL on failure
*/
void* siril_realloc(void* ptr, size_t new_size) {
#ifdef _WIN32
	if (ptr == NULL) {
		// realloc(NULL, size) is equivalent to malloc(size)
		return siril_malloc(new_size);
	}

	if (new_size == 0) {
		// realloc(ptr, 0) is equivalent to free(ptr)
		siril_free(ptr);
		return NULL;
	}

	// Get the header to find current allocation method
	AllocHeader* header = (AllocHeader*)((char*)ptr - sizeof(AllocHeader));

	// Calculate old usable size (size minus header and any alignment padding)
	// For simplicity, we'll just copy conservatively based on new_size
	size_t old_size;
	if (header->method == 0) {
		// malloc allocation - size includes header
		old_size = header->size - sizeof(AllocHeader);
	} else {
		// OS allocation - we allocated extra for alignment, but we don't know
		// exactly how much the user originally requested. Be conservative.
		old_size = header->size - sizeof(AllocHeader) - SIMD_ALIGNMENT;
	}

	// Check if we're crossing the threshold boundary
	int old_is_large = (old_size >= LARGE_ALLOC_THRESHOLD);
	int new_is_large = (new_size >= LARGE_ALLOC_THRESHOLD);

	if (old_is_large == new_is_large && header->method == 0) {
		// Both small, both using malloc - use realloc directly
		void* old_malloc_ptr = header->original_ptr;
		void* new_malloc_ptr = realloc(old_malloc_ptr, new_size + sizeof(AllocHeader));

		if (new_malloc_ptr == NULL) {
			return NULL;  // Realloc failed, original pointer still valid
		}

		// Update header
		header = (AllocHeader*)new_malloc_ptr;
		header->original_ptr = new_malloc_ptr;
		header->size = new_size + sizeof(AllocHeader);
		header->method = 0;

		return (void*)((char*)new_malloc_ptr + sizeof(AllocHeader));
	} else {
		// Crossing threshold or using OS allocation - allocate new, copy, free old
		void* new_ptr = siril_malloc(new_size);

		if (new_ptr == NULL) {
			return NULL;  // Allocation failed, original pointer still valid
		}

		// Copy old data (copy the minimum of old and new size)
		size_t copy_size = (old_size < new_size) ? old_size : new_size;
		memcpy(new_ptr, ptr, copy_size);

		// Free old allocation
		siril_free(ptr);

		return new_ptr;
	}
#else
	// POSIX: use standard realloc
	return realloc(ptr, new_size);
#endif
}

#endif /* OPTIMAL_ALLOC_H */
