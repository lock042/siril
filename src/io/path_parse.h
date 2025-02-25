/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_IO_PATH_PARSE_H_
#define SRC_IO_PATH_PARSE_H_

#include "core/siril.h"

// negative errors are warnings
typedef enum {
	// warning for wilcard read mode
	PATHPARSE_ERR_WILDCARD_SYNTAX_NOFAIL = -12,
	PATHPARSE_ERR_HITS_ALL_NEWER = -11,
	PATHPARSE_ERR_MORE_THAN_ONE_HIT = -10,
	// warnings for write_nofail mode (exact opposite of errors betw 1 and 10)
	PATHPARSE_ERR_HEADER_NULL_NOFAIL = -1,
	PATHPARSE_ERR_WRONG_RESERVED_KEYWORD_NOFAIL = -2,
	PATHPARSE_ERR_LIBRARY_NOTDEFINED_NOFAIL = -3,
	PATHPARSE_ERR_KEY_NOT_FOUND_NOFAIL = -4,
	PATHPARSE_ERR_WRONG_DATE_NOFAIL = -5,
	PATHPARSE_ERR_WRONG_WCS_NOFAIL = -6,
	PATHPARSE_ERR_UNSUPPORTED_FORMAT_NOFAIL = -7,
	PATHPARSE_ERR_NOSEQLOADED_NOFAIL = -8,
	PATHPARSE_ERR_BADSTRING_NOFAIL = -9,
	// no error
	PATHPARSE_ERR_OK = 0,
	// parsing errors for classic read or write modes
	PATHPARSE_ERR_HEADER_NULL = 1,
	PATHPARSE_ERR_WRONG_RESERVED_KEYWORD = 2,
	PATHPARSE_ERR_LIBRARY_NOTDEFINED = 3,
	PATHPARSE_ERR_KEY_NOT_FOUND = 4,
	PATHPARSE_ERR_WRONG_DATE = 5,
	PATHPARSE_ERR_WRONG_WCS = 6,
	PATHPARSE_ERR_UNSUPPORTED_FORMAT = 7,
	PATHPARSE_ERR_NOSEQLOADED = 8,
	PATHPARSE_ERR_BADSTRING = 9,
	// errors for read mode with wildcards
	PATHPARSE_ERR_NO_HIT_FOUND = 10,
	PATHPARSE_ERR_NO_DIR = 11,
	PATHPARSE_ERR_WILDCARD_SYNTAX = 12,
	// internal errors
	PATHPARSE_ERR_WRONG_CALL = 20,
	PATHPARSE_ERR_TMPFIT = 21
} pathparse_errors;

typedef enum {
	PATHPARSE_MODE_READ,
	PATHPARSE_MODE_WRITE,
	PATHPARSE_MODE_WRITE_NOFAIL
} pathparse_mode;

gchar *path_parse(fits *fit, const gchar *expression, pathparse_mode mode, int *status);
gchar *update_header_and_parse(fits *fit, gchar *expression, pathparse_mode mode, gboolean createdir, int *status);
pathparse_errors read_key_from_header_text(gchar **headers, gchar *key, double *numvalue, gchar *strvalue);

#endif /* SRC_IO_PATH_PARSE_H_ */
