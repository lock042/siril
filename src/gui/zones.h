#ifndef _ZONES_H_
#define _ZONES_H_

#include "core/siril.h"

int get_number_of_zones();
int get_side(const stacking_zone *zone);
int point_is_inside_zone(int px, int py, const stacking_zone *zone);
gboolean zone_is_too_close(int x, int y, int size);

void add_stacking_zone(double x, double y, double half_side);
void remove_stacking_zones();

#endif
