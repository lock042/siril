#ifndef GUI_USER_POLYGONS_H
#define GUI_USER_POLYGONS_H

typedef struct {
	int id;
	int n_points;
	point *points;
	GdkRGBA color;
	gboolean fill;
	gchar *legend;
} __attribute__((packed)) UserPolygon;

#define MAX_POLYGON_POINTS 100

UserPolygon *find_polygon_by_id(int id);
int get_unused_polygon_id(void);
int add_user_polygon(point *points, int num_points, const GdkRGBA *color, gboolean fill);
gboolean delete_user_polygon(int id);
void clear_user_polygons(void);
UserPolygon* deserialize_polygon(const uint8_t *data, size_t size);
uint8_t* serialize_polygon(UserPolygon *polygon, size_t *size);
uint8_t* serialize_polygon_list(GList *polygons, size_t *out_size);

#endif // GUI_USER_POLYGONS_H
