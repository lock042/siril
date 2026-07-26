#ifndef SRC_IO_AVI_PREVIEW_H_
#define SRC_IO_AVI_PREVIEW_H_

#include <glib.h>

/* Decode just the first frame of a video container (AVI / MP4 / MOV /
 * MKV / ...) via libavformat / libavcodec.  Returns a malloc'd RGB888
 * byte buffer at the file's native frame size, plus the dimensions and
 * a g_strdup'd description string ("WxH pixels[\nN frames]").  Caller
 * owns the buffer (free()) and the descriptor (g_free).
 *
 * Returns NULL on error or when ffmpeg is unavailable at build time.
 *
 * Designed for file-browser previews: avoids the multi-second indexing
 * pass that FFMS2 performs.  Typically only a few hundred KB of the
 * source file are read.  */
guchar *extract_thumbnail_from_avi(const char *filename, gchar **descr,
                                    int *width_out, int *height_out);

#endif /* SRC_IO_AVI_PREVIEW_H_ */
