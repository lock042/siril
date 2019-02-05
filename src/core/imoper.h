#ifndef _IMOPER_H_
#define _IMOPER_H_

#include "core/siril.h"

int threshlo(fits *fit, int level);
int threshhi(fits *fit, int level);
int nozero(fits *fit, int level);
int  off(fits *a, int level);

int soper(fits *a, double scalar, char oper);
int imoper(fits *a, fits *b, char oper);

int sub_background(fits* image, fits* background, int layer);
int addmax(fits *a, fits *b);
int siril_fdiv(fits *a, fits *b, float coef);
int siril_ndiv(fits *a, fits *b);

int loglut(fits *fit);
int asinhlut(fits *fit, double beta, double offset, gboolean RGBspace);
int ddp(fits *a, int level, float coeff, float sigma);

#endif
