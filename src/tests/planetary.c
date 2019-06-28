#include <stdio.h>
#include <stdlib.h>
/* test for algos/planetary.c */

struct weighted_AP {
	int zone_index;
	double distance;
};

/* System under test, copied from algos/planetary.c for simpler link */
void check_closest_list(struct weighted_AP *list_for_this_point,
		double distance, int zone_idx, int max_AP) {
	int i, max = 0, max_idx = 0, found_closer = 0;
	for (i = 0; i < max_AP; i++) {
		double ap_dist = list_for_this_point[i].distance;
		if (ap_dist >= 0) {
			// we search the max to replace it by the new
			if (ap_dist > max) {
				max = ap_dist;
				max_idx = i;
			}
			if (distance < ap_dist)
				found_closer = 1;
		}
		else {
			found_closer = 1;
			max_idx = i;
			break;
		}
	}
	if (found_closer) {
		list_for_this_point[max_idx].distance = distance;
		list_for_this_point[max_idx].zone_index = zone_idx;
	}
}

static int test_list_content(struct weighted_AP *list, int *indices, int size) {
	if (!list) return 1;
	int i;
	for (i = 0; i < size; i++) {
		if (indices[i] > 0 && list[i].zone_index != indices[i])
			return 1;
	}
	return 0;
}

int test_closest_list() {
	int nb_closest_AP = 3, ap, ret;
	struct weighted_AP *closest_zones_map = malloc(1 * 1 * nb_closest_AP * sizeof(struct weighted_AP));
	struct weighted_AP *list = closest_zones_map; // (one pixel image)

	for (ap = 0; ap < nb_closest_AP; ap++) {    // init the struct
		list[ap].distance = -1.0;
		list[ap].zone_index = -1.0;
	}

	check_closest_list(list, 0.2, 1, nb_closest_AP);
	ret = test_list_content(list, (int[3]) { 1, -1, -1 }, nb_closest_AP);
	if (ret) return -1;
	check_closest_list(list, 1.0, 2, nb_closest_AP);
	ret = test_list_content(list, (int[3]) { 1, 2, -1 }, nb_closest_AP);
	if (ret) return -2;
	check_closest_list(list, 0.5, 3, nb_closest_AP);
	ret = test_list_content(list, (int[3]) { 1, 2, 3 }, nb_closest_AP);
	if (ret) return -3;
	check_closest_list(list, 1.5, 4, nb_closest_AP);
	ret = test_list_content(list, (int[3]) { 1, 2, 3 }, nb_closest_AP);
	if (ret) return -4;
	check_closest_list(list, 0.4, 5, nb_closest_AP);
	ret = test_list_content(list, (int[3]) { 1, 5, 3 }, nb_closest_AP);
	if (ret) return -5;
	check_closest_list(list, 0.1, 6, nb_closest_AP);
	ret = test_list_content(list, (int[3]) { 1, 5, 6 }, nb_closest_AP);
	if (ret) return -6;

	return 0;
}

int main() {
	int retval = test_closest_list();
	if (retval)
		fprintf(stdout, "failed at %d\n", -retval);
	else fprintf(stdout, "success\n");
	return retval;
}
