//#include <criterion/criterion.h>

#include "../core/siril.h"
#include "stacking/stacking.h"
#include <stdio.h>

cominfo com;	// the main data struct

/* This is the test file for the function that allocates data block to thread
 * for siril median or mean stacking.
 */

#define CHECK(cond, ...) \
	if (cond) { \
		fprintf(stderr, __VA_ARGS__); \
		return 1; \
	}

int check_that_blocks_cover_the_image(long naxes[3], struct _image_block *blocks, int nb_blocks) {
	// assuming that blocks are adjacent 
	long ok_up_to_y = 0L, ok_up_to_z = 0L;
	for (int i = 0; i < nb_blocks; i++) {
		if (blocks[i].start_row != ok_up_to_y)
			return 0;
		ok_up_to_y = blocks[i].end_row + 1L;
		if (ok_up_to_y == naxes[1]) {
			ok_up_to_y = 0L;
			ok_up_to_z++;
		}
	}
	return ok_up_to_y == 0L && ok_up_to_z == naxes[2];
}

int test1() {
	// intputs
	long naxes[] = { 1000L, 1000L, 1L };
	long max_rows = 1001L;
	int nb_threads = 1;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;

	/* case 1: single channel images, enough memory, one thread */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks != 1, "number of blocks returned is %d (expected 1)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	fprintf(stdout, "* test passed *\n");
	return 0;
}

int test2() {
	// intputs
	long naxes[] = { 1000L, 1000L, 1L };
	long max_rows = 1001L;
	int nb_threads = 8;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;

	/* case 2: single channel images, enough memory, 8 threads */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks != 8, "number of blocks returned is %d (expected 8)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	fprintf(stdout, "* test passed *\n");
	return 0;
}

int test3() {
	// intputs
	long naxes[] = { 1000L, 1000L, 1L };
	long max_rows = 999L;
	int nb_threads = 1;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;


	/* case 3: single channel images, not enough memory, 1 thread */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks != 2, "number of blocks returned is %d (expected 2)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	fprintf(stdout, "* test passed *\n");
	return 0;
}

int main() {
	int retval = 0;
	retval |= test1();
	retval |= test2();
	retval |= test3();
	if (retval)
		fprintf(stderr, "TESTS FAILED\n");
	else fprintf(stderr, "ALL TESTS PASSED\n");
	return retval;
}


/* inputs : image 3D dimensions, number of images, channel depth, configured
 *	memory limits, actual memory available, configured thread limit
 * outputs : number of blocks and their size
 */
