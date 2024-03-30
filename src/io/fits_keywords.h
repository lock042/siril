/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_IO_FITS_KEYWORDS_H_
#define SRC_IO_FITS_KEYWORDS_H_

enum keywords_type {
	KTYPE_BOOL,
	KTYPE_INT,
	KTYPE_UINT,
	KTYPE_DOUBLE,
	KTYPE_STR,
	KTYPE_DATE
};


struct keywords_access {
	const char *group;	// group name
	const char *key;	// key name
	enum keywords_type type;// type of the keyword
	const char *comment;	// comment
	void *data;		// pointer to the data in keyword struct
};

int save_fits_keywords(fits *fit);


#endif /* SRC_IO_FITS_KEYWORDS_H_ */
