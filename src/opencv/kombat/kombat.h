#ifndef _KOMBAT_H_
#define _KOMBAT_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * KOMBAT: Kind of Minimal/Basic Alignment Tool
 * https://discuss.pixls.us/t/alignment-over-planetary-captures/30103
 *
 * Author:  Frédéric Trouche
 * e-mail: fred@linuxtribe.org
 *
 * Copyright (2022): Frédéric Trouche
 *
 */

typedef  struct reg_kombat_struct reg_kombat;

struct reg_kombat_struct {
	float dx;
	float dy;
};

int kombat_find_template(int idx, struct registration_args *args, fits *templ, fits *image,  reg_kombat *reg_param,
										   reg_kombat *ref_align, void **vcache);

void kombat_done(void **vcache);

#ifdef __cplusplus
}
#endif

#endif //_KOMBAT_H_
