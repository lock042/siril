#pragma once

#ifdef __cplusplus
#define EXTERNC extern "C" {
#else
#define EXTERNC
#endif
EXTERNC void portable_usleep(size_t usec);
#ifdef __cplusplus
}
#endif
