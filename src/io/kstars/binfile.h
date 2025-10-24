/* This file comes from KStars tools to handle star catalogues.
 * It has been modified to remove the need for a fixed HTM level and
 * to be C instead of C++.
 */
/*
    SPDX-FileCopyrightText: 2008 Akarsh Simha <akarshsimha@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

/* NOTE: This file is a NON portable C file that need not be built for
   KStars to be built.  This file is useful only to generate binary
   data files and ensure that they comply with KStars' format, and
   hence need not be cross platform.

   This file does not use Qt or KDE libraries as it is just provided
   here as a tool to build and test data files for KStars

   Hence, we shall hide this file from the watchful eyes of Krazy
*/
//krazy:skip

#ifndef BINFILE_H
#define BINFILE_H

// Siril addition in order to support siril_* alloc functions
#include "core/siril_alloc.h"

#include <stdint.h>
#include <sys/types.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "byteorder.h"

/* NOTE: HTM_LEVEL and other HTM-related stuff must be defined before using this header */

// Bogus Define
#ifndef INDEX_ENTRY_SIZE
#define INDEX_ENTRY_SIZE 12
#endif

/*
 * enum listing out various possible data types
 */

enum dataType
{
	DT_CHAR,   /* Character */
	DT_INT8,   /* 8-bit Integer */
	DT_UINT8,  /* 8-bit Unsigned Integer */
	DT_INT16,  /* 16-bit Integer */
	DT_UINT16, /* 16-bit Unsigned Integer */
	DT_INT32,  /* 32-bit Integer */
	DT_UINT32, /* 32-bit Unsigned Integer */
	DT_CHARV,  /* Fixed-length array of characters */
	DT_STR,    /* Variable length array of characters, either terminated by nullptr or by the limit on field size */
	DT_SPCL = 128 /* Flag indicating that the field requires special treatment (eg: Different bits may mean different things) */
};

/*
 * struct to store the description of a field / data element in the binary files
 */

typedef struct
{
	char name[10];
	int8_t size;
	uint8_t type;
	int32_t scale;
} dataElement;

void charv2str(char *str, char *charv, int n)
{
	int i;
	for (i = 0; i < n; ++i)
	{
		*str = *charv;
		str++;
		charv++;
	}
	*str = '\0';
}

const char *data_element_type_to_string(uint8_t type) {
	switch (type) {
		case DT_CHAR:   return "Character";
		case DT_INT8:   return "8-bit Integer";
		case DT_UINT8:  return "8-bit Unsigned Integer";
		case DT_INT16:  return "16-bit Integer";
		case DT_UINT16: return "16-bit Unsigned Integer";
		case DT_INT32:  return "32-bit Integer";
		case DT_UINT32: return "32-bit Unsigned Integer";
		case DT_CHARV:  return "Fixed-length array of characters";
		case DT_STR:    return "String";
		case DT_SPCL:	return "Special type";
		default:	return "Unknown type";
	}
}

void displayDataElementDescription(dataElement *e)
{
	char str[11];
	charv2str(str, e->name, 10);
	const char *type = data_element_type_to_string(e->type);
	fprintf(stdout, "\tData field %s of type %s, size %d, scale %d\n", str, type, e->size, e->scale);
	/*printf("\nData Field:\n");
	printf("  Name: %s\n", str);
	printf("  Size: %d\n", e->size);
	printf("  Type: %d\n", e->type);
	printf("  Scale: %d\n", e->scale);*/
}

// NOTE: Ineffecient. Not to be used for high-productivity
// applications
void swapbytes(char byteswap, void *ptr, int nbytes)
{
	char *destptr;
	char *i;

	if (!byteswap)
		return;

	destptr = (char *)siril_malloc(nbytes);
	i       = ((char *)ptr + (nbytes - 1));
	while (i >= (char *)ptr)
	{
		*destptr = *i;
		++destptr;
		--i;
	}

	destptr -= nbytes;

	memcpy(ptr, (void *)destptr, nbytes);
	siril_free(destptr);
}

/*
 * Convert a string to an int32_t with a double as an intermediate
 * i : A pointer to the target int32_t
 * str : A pointer to the string that carries the data
 * ndec : Number of decimal places to truncate to
 */

int str2int32(int32_t *i, const char *str, int ndec)
{
	double dbl;

	if (!i)
		return 0;

	dbl = atof(str);

	*i = (int32_t)(round(dbl * pow(10, ndec)));

	return 1;
}

/*
 * Convert a string to an int16_t with a double as an intermediate
 * i : A pointer to the target int16_t
 * str : The string that carries the data
 * ndec : Number of decimal places to truncate to
 */

int str2int16(int16_t *i, const char *str, int ndec)
{
	double dbl;

	if (!i || !str)
		return 0;

	dbl = atof(str);

	*i = (int16_t)(round(dbl * pow(10, ndec)));

	return 1;
}

/*
 * Convert a string into a character array for n characters
 * a : The target array
 * str : The string that carries the data
 * n : Number of characters to convert
 */

int str2charv(char *a, const char *str, int n)
{
	if (!a || !str)
		return 0;

	int ret = 1;

	for (int i = 0; i < n; ++i)
	{
		a[i] = ((ret < 0) ? '\0' : str[i]);
		if (str[i] == '\0') /* We can do this safely because we aren't storing binary data in the DB */
			ret = -1;
	}
	return ret;
}

/*
 * Write one data element description into a binary file header.
 *
 * f    : Handle of file to write into
 * name : Name of the field, as a string. Max as specified in struct dataElement
 * size : Size (in bytes) of the field
 * type : Type of the field, as specified in enum dataType
 * scale : Scale factor used for conversion of fixed-point reals to integers. N/A to DT_CHARV, DT_STR and DT_CHAR
 */

int writeDataElementDescription(FILE *f, char *name, int8_t size, enum dataType type, int32_t scale)
{
	dataElement de;

	if (!f || !name)
		return 0;

	str2charv(de.name, name, 10);
	de.size  = size;
	de.type  = type;
	de.scale = scale;
	fwrite(&de, sizeof(dataElement), 1, f);
	return 1;
}

int writeIndexEntry(FILE *hf, uint32_t trixel_id, uint32_t offset, uint32_t nrec)
{
	if (!hf)
		return 0;

	fwrite(&trixel_id, 4, 1, hf);
	fwrite(&offset, 4, 1, hf);
	fwrite(&nrec, 4, 1, hf);

	/* Put this just for safety, in case we change our mind - we should avoid screwing things up */
	if (4 + 4 + 4 != INDEX_ENTRY_SIZE)
	{
		fprintf(stderr, "CODE ERROR: 4 + 4 + 4 != INDEX_ENTRY_SIZE\n");
	}

	return 1;
}

#endif
