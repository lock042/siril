/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
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
#ifndef SRC_IO_MY_WINFILE_H_
#define SRC_IO_MY_WINFILE_H_

#ifdef _WIN32
#include <sys/stat.h>
/* my_winfile.c exports, should not be used outside mysys */
extern int     my_win_open(const char *path, int oflag);
extern int      my_win_close(int fd);
extern size_t   my_win_read(int fd, uchar *buffer, size_t  count);
extern size_t   my_win_write(int fd, const uchar *buffer, size_t count);
extern size_t   my_win_pread(int fd, uchar *buffer, size_t count,
                             my_off_t offset);
extern size_t   my_win_pwrite(int fd, const uchar *buffer, size_t count,
                              my_off_t offset);
extern my_off_t my_win_lseek(int fd, my_off_t pos, int whence);
extern int      my_win_chsize(int fd,  my_off_t newlength);
extern FILE*    my_win_fopen(const char *filename, const char *type);
extern int     my_win_fclose(FILE *file);
extern int     my_win_fileno(FILE *file);
extern FILE*    my_win_fdopen(int intdes, const char *type);
extern int      my_win_stat(const char *path, struct _stati64 *buf);
extern int      my_win_fstat(int fd, struct _stati64 *buf);
extern int      my_win_fsync(int fd);
extern int     my_win_dup(int fd);
extern int     my_win_sopen(const char *path, int oflag, int shflag, int perm);
extern int     my_open_osfhandle(HANDLE handle, int oflag);
#endif


#endif /* SRC_IO_MY_WINFILE_H_ */
