#if defined(__APPLE__) && defined(__MACH__)
#ifndef FILE_SYSTEM_INFO_H
#define FILE_SYSTEM_INFO_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int64_t macos_find_space(const char *path);

#ifdef __cplusplus
}
#endif

#endif // FILE_SYSTEM_INFO_H
#endif // defined(__APPLE__) && defined(__MACH__)
