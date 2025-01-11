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

#ifndef SRC_CORE_INITFILE_H_
#define SRC_CORE_INITFILE_H_

#include <glib.h>

enum token_index {
	WD = 0,		/* Working Directory */
	BAY = 1,	/* Bayer settings */
	PRE = 2,	/* Preprocessing settings */
	REG = 3,	/* Registration settings */
	STK = 4,	/* Stacking settings */
	AST = 5,	/* Astrometry settings */
	PTM = 6,	/* Photometry settings */
	MISC = 7,	/* Miscellaneous settings */
	CMP = 8,	/* Compression settings */
	NOTOK
};

int writeinitfile();
int checkinitfile();

int readinitfile(char *path);	// to import some settings

int read_keyfile(GKeyFile *kf);

#endif /* SRC_CORE_INITFILE_H_ */
