// SPDX-License-Identifier: GPL-3.0-or-later
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


#ifdef _WIN32
#include <windows.h>
#else
#include <stdlib.h>
#include <string.h>
#endif

#include <stddef.h>
#include "core/siril.h"
	
/**
 * Allocates a block of memory of the specified size.
 * On Windows, uses VirtualAlloc for large allocations.
 * On other platforms, uses standard malloc.
 * 
 * @param size Number of bytes to allocate
 * @return Pointer to allocated memory, or NULL on failure
 */
void* siril_malloc(size_t size) {
	if (size == 0) {
		return NULL;
	}

#ifdef _WIN32
	// Use VirtualAlloc on Windows for better handling of large allocations
	// MEM_COMMIT | MEM_RESERVE allocates and commits memory in one step
	// PAGE_READWRITE allows read/write access
	void* ptr = VirtualAlloc(NULL, size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
	siril_debug_print("Calling VirtualAlloc\n");
	return ptr;  // Returns NULL on failure, just like malloc
#else
	return malloc(size);
#endif
}

/**
 * Allocates a block of memory for an array and initializes all bytes to zero.
 * On Windows, uses VirtualAlloc which already zeros memory.
 * On other platforms, uses standard calloc.
 * 
 * @param num Number of elements
 * @param size Size of each element in bytes
 * @return Pointer to allocated and zeroed memory, or NULL on failure
 */
void* siril_calloc(size_t num, size_t size) {
	if (num == 0 || size == 0) {
		return NULL;
	}

	// Check for overflow
	size_t total_size = num * size;
	if (total_size / num != size) {
		return NULL;  // Overflow occurred
	}

#ifdef _WIN32
	// VirtualAlloc already zeros the allocated memory
	void* ptr = VirtualAlloc(NULL, total_size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
	siril_debug_print("Calling VirtualAlloc\n");
	return ptr;  // Memory is already zeroed, returns NULL on failure
#else
	return calloc(num, size);
#endif
}

/**
 * Frees memory allocated by siril_malloc or siril_calloc.
 * Must be used to free memory allocated by siril_malloc/siril_calloc.
 * 
 * @param ptr Pointer to memory to free (can be NULL)
 */
void siril_free(void* ptr) {
	if (ptr == NULL) {
		return;
	}

#ifdef _WIN32
	// VirtualFree with MEM_RELEASE frees the entire region
	VirtualFree(ptr, 0, MEM_RELEASE);
#else
	free(ptr);
#endif
}

// /**
//  * Reallocates memory to a new size.
//  * On Windows, allocates new memory and copies data.
//  * On other platforms, uses standard realloc.
//  * 
//  * Note: On Windows, this requires knowing the old size for proper copying.
//  * For production use, consider maintaining a size tracking mechanism.
//  * 
//  * @param ptr Pointer to previously allocated memory
//  * @param new_size New size in bytes
//  * @return Pointer to reallocated memory, or NULL on failure
//  */
// void* siril_realloc(void* ptr, size_t new_size) {
// 	if (ptr == NULL) {
// 		return siril_malloc(new_size);
// 	}
	
// 	if (new_size == 0) {
// 		siril_free(ptr);
// 		return NULL;
// 	}

// #ifdef _WIN32
// 	// VirtualAlloc doesn't support reallocation directly
// 	// Allocate new block and copy data
// 	void* new_ptr = siril_malloc(new_size);
// 	if (new_ptr == NULL) {
// 		return NULL;
// 	}
	
// 	// Query the old allocation size
// 	MEMORY_BASIC_INFORMATION mbi;
// 	if (VirtualQuery(ptr, &mbi, sizeof(mbi))) {
// 		size_t old_size = mbi.RegionSize;
// 		size_t copy_size = (old_size < new_size) ? old_size : new_size;
// 		memcpy(new_ptr, ptr, copy_size);
// 	}
	
// 	siril_free(ptr);
// 	return new_ptr;
// #else
// 	return realloc(ptr, new_size);
// #endif
// }
