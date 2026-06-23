#ifndef GUI_USER_POLYGONS_H
#define GUI_USER_POLYGONS_H

#include <stdint.h>
#include "core/siril.h"

// Very liberal limit, purely to avoid unlimited g_malloc0 calls
#define MAX_POLYGON_POINTS 1e6

UserPolygon *find_polygon_by_id(int id);
void free_user_polygon(gpointer data);
int get_unused_polygon_id(void);
/* color is RGBA as double[4], indices 0=R 1=G 2=B 3=A, values in [0, 1] */
int add_user_polygon(point *points, int num_points, const double *color, gboolean fill);
gboolean delete_user_polygon(int id);
void clear_user_polygons(void);
UserPolygon* deserialize_polygon(const uint8_t *data, size_t size);
UserPolygon* create_user_polygon_from_points(GSList *point_list);
int add_existing_polygon(UserPolygon* poly, const double *color, gboolean fill);
uint8_t* serialize_polygon(UserPolygon *polygon, size_t *size);
uint8_t* serialize_polygon_list(GSList *polygons, size_t *out_size);

#endif // GUI_USER_POLYGONS_H
