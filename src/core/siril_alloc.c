#include "core/siril_alloc.h"

#ifdef _WIN32
#include <windows.h>
#endif
/**
* Allocate memory optimally for the given size.
* On POSIX systems, uses siril_malloc (which automatically switches to mmap for large allocations).
* On Windows, small allocations use siril_malloc/siril_free, large allocations use VirtualAlloc
* to avoid heap fragmentation, with SIMD-aligned pointers (64-byte).
*
* @param size Number of bytes to allocate
* @return Pointer to allocated memory, or NULL on failure
*/
void* siril_malloc(size_t size) {
#ifdef _WIN32
    /* existing Windows code unchanged */
    if (size >= LARGE_ALLOC_THRESHOLD) {
        size_t page_size = 4096;
        size_t total_size = size + sizeof(AllocHeader) + SIMD_ALIGNMENT;
        size_t aligned_size = (total_size + page_size - 1) & ~(page_size - 1);
        LPVOID os_ptr = VirtualAlloc(NULL, aligned_size,
                                     MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
        if (os_ptr == NULL)
            return NULL;
        uintptr_t after_header = (uintptr_t)os_ptr + sizeof(AllocHeader);
        uintptr_t user_ptr = (after_header + SIMD_ALIGNMENT - 1) & ~(SIMD_ALIGNMENT - 1);
        AllocHeader* header = (AllocHeader*)(user_ptr - sizeof(AllocHeader));
        header->original_ptr = os_ptr;
        header->size = aligned_size;
        header->method = 1;
        return (void*)user_ptr;
    } else {
        size_t total_size = size + sizeof(AllocHeader) + SIMD_ALIGNMENT;
        void* raw = _aligned_malloc(total_size, SIMD_ALIGNMENT);
        if (!raw)
            return NULL;
        uintptr_t after_header = (uintptr_t)raw + sizeof(AllocHeader);
        uintptr_t user_ptr = (after_header + SIMD_ALIGNMENT - 1) & ~(SIMD_ALIGNMENT - 1);
        AllocHeader* header = (AllocHeader*)(user_ptr - sizeof(AllocHeader));
        header->original_ptr = raw;
        header->size = total_size;
        header->method = 0;
        return (void*)user_ptr;
    }

#else
    /* POSIX: explicit malloc vs mmap logic */
    if (size >= LARGE_ALLOC_THRESHOLD) {
        size_t page_size = sysconf(_SC_PAGESIZE);
        size_t total_size = size + sizeof(AllocHeader) + SIMD_ALIGNMENT;
        size_t aligned_size = (total_size + page_size - 1) & ~(page_size - 1);
        void* os_ptr = mmap(NULL, aligned_size, PROT_READ | PROT_WRITE,
                            MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (os_ptr == MAP_FAILED)
            return NULL;
        uintptr_t after_header = (uintptr_t)os_ptr + sizeof(AllocHeader);
        uintptr_t user_ptr = (after_header + SIMD_ALIGNMENT - 1) & ~(SIMD_ALIGNMENT - 1);
        AllocHeader* header = (AllocHeader*)(user_ptr - sizeof(AllocHeader));
        header->original_ptr = os_ptr;
        header->size = aligned_size;
        header->method = 1;
        return (void*)user_ptr;
    } else {
        size_t total_size = size + sizeof(AllocHeader) + SIMD_ALIGNMENT;
        void* raw = malloc(total_size);
        if (!raw)
            return NULL;
        uintptr_t after_header = (uintptr_t)raw + sizeof(AllocHeader);
        uintptr_t user_ptr = (after_header + SIMD_ALIGNMENT - 1) & ~(SIMD_ALIGNMENT - 1);
        AllocHeader* header = (AllocHeader*)(user_ptr - sizeof(AllocHeader));
        header->original_ptr = raw;
        header->size = total_size;
        header->method = 0;
        return (void*)user_ptr;
    }
#endif
}

/**
 * Allocate and zero-initialize memory optimally.
 */
void* siril_calloc(size_t nmemb, size_t size) {
#ifdef _WIN32
    if (nmemb != 0 && size > SIZE_MAX / nmemb)
        return NULL;
    size_t total = nmemb * size;
    void* ptr = siril_malloc(total);
    if (ptr)
        memset(ptr, 0, total);
    return ptr;
#else
    if (nmemb != 0 && size > SIZE_MAX / nmemb)
        return NULL;
    size_t total = nmemb * size;
    void* ptr = siril_malloc(total);
    if (ptr)
        memset(ptr, 0, total);
    return ptr;
#endif
}

/**
 * Free memory allocated with siril_malloc/siril_calloc.
 */
void siril_free(void* ptr) {
    if (!ptr)
        return;

    AllocHeader* header = (AllocHeader*)((uintptr_t)ptr - sizeof(AllocHeader));

#ifdef _WIN32
    if (header->method == 1)
        VirtualFree(header->original_ptr, 0, MEM_RELEASE);
    else
        _aligned_free(header->original_ptr);
#else
    if (header->method == 1)
        munmap(header->original_ptr, header->size);
    else
        free(header->original_ptr);
#endif
}

/**
 * Reallocate memory while preserving allocation type logic.
 */
void* siril_realloc(void* ptr, size_t new_size) {
#ifdef _WIN32
    if (!ptr)
        return siril_malloc(new_size);
    if (new_size == 0) {
        siril_free(ptr);
        return NULL;
    }

    AllocHeader* header = (AllocHeader*)((uintptr_t)ptr - sizeof(AllocHeader));
    size_t old_size = header->size - sizeof(AllocHeader) - SIMD_ALIGNMENT;
    int old_large = (old_size >= LARGE_ALLOC_THRESHOLD);
    int new_large = (new_size >= LARGE_ALLOC_THRESHOLD);

    if (old_large == new_large && header->method == 0) {
        void* new_raw = _aligned_realloc(header->original_ptr,
                                         new_size + sizeof(AllocHeader) + SIMD_ALIGNMENT,
                                         SIMD_ALIGNMENT);
        if (!new_raw)
            return NULL;
        uintptr_t after_header = (uintptr_t)new_raw + sizeof(AllocHeader);
        uintptr_t user_ptr = (after_header + SIMD_ALIGNMENT - 1) & ~(SIMD_ALIGNMENT - 1);
        AllocHeader* new_header = (AllocHeader*)(user_ptr - sizeof(AllocHeader));
        new_header->original_ptr = new_raw;
        new_header->size = new_size;
        new_header->method = 0;
        return (void*)user_ptr;
    } else {
        void* new_ptr = siril_malloc(new_size);
        if (!new_ptr)
            return NULL;
        size_t copy_size = (old_size < new_size) ? old_size : new_size;
        memcpy(new_ptr, ptr, copy_size);
        siril_free(ptr);
        return new_ptr;
    }

#else
    if (!ptr)
        return siril_malloc(new_size);
    if (new_size == 0) {
        siril_free(ptr);
        return NULL;
    }

    AllocHeader* header = (AllocHeader*)((uintptr_t)ptr - sizeof(AllocHeader));
    size_t old_size = header->size - sizeof(AllocHeader) - SIMD_ALIGNMENT;
    int old_large = (old_size >= LARGE_ALLOC_THRESHOLD);
    int new_large = (new_size >= LARGE_ALLOC_THRESHOLD);

    if (old_large == new_large && header->method == 0) {
        void* old_raw = header->original_ptr;
        void* new_raw = realloc(old_raw, new_size + sizeof(AllocHeader) + SIMD_ALIGNMENT);
        if (!new_raw)
            return NULL;
        uintptr_t after_header = (uintptr_t)new_raw + sizeof(AllocHeader);
        uintptr_t user_ptr = (after_header + SIMD_ALIGNMENT - 1) & ~(SIMD_ALIGNMENT - 1);
        AllocHeader* new_header = (AllocHeader*)(user_ptr - sizeof(AllocHeader));
        new_header->original_ptr = new_raw;
        new_header->size = new_size;
        new_header->method = 0;
        return (void*)user_ptr;
    } else {
        void* new_ptr = siril_malloc(new_size);
        if (!new_ptr)
            return NULL;
        size_t copy_size = (old_size < new_size) ? old_size : new_size;
        memcpy(new_ptr, ptr, copy_size);
        siril_free(ptr);
        return new_ptr;
    }
#endif
}
