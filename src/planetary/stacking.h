#ifndef _PLANETARY_STACKING_H_
#define _PLANETARY_STACKING_H_

int the_old_local_multipoint_sum_stacking(struct mpr_args *args);
int the_new_local_multipoint_sum_stacking(struct mpr_args *args);

void add_image_zone_to_stacking_sum(fits *fit,
		const stacking_zone *zone, int frame,
		unsigned long *sum[3], int *count[3]);

#endif
