/* Temporary diagnostic: stress the FLIS GPU compose tile cache the way the
 * GUI does — small byte budget (mid-render eviction), mip changes from
 * zoom, stretch (lo/hi) changes, layer invalidations and prefetch idles
 * between renders, panning visible rects.  Reproducer hunt for the
 * snapshot-time SIGSEGV in ensure_tile (crash report 2026-06-12). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include <gtk/gtk.h>

#include "core/siril.h"
#include "io/image_format_fits.h"
#include "io/image_format_flis.h"
#include "io/flis_compose.h"
#include "gui-gtk4/gui_state.h"
#include "gui-gtk4/flis_gpu_compose.h"

cominfo com;
fits *gfit;
guiinfo gui;

/* Worker-thread invalidator: mimics generic_image_worker ending with
 * notify_gfit_data_modified() -> flis_invalidate_composite() while the
 * main thread is mid-snapshot. */
static volatile gboolean g_run_invalidator = TRUE;
static gpointer invalidator_thread(gpointer data) {
	guint32 lcg = 99;
	while (g_run_invalidator) {
		lcg = lcg * 1103515245u + 12345u;
		if ((lcg >> 16) & 1)
			flis_gpu_compose_invalidate_all();
		else
			flis_gpu_compose_invalidate_layer(1 + (int)((lcg >> 9) % 4));
		g_usleep(((lcg >> 20) % 400) + 50);
	}
	return NULL;
}

static void render_once(guint W, guint H, float zoom,
                        float vx, float vy, float vw, float vh) {
	GtkSnapshot *snap = gtk_snapshot_new();
	graphene_rect_t dst = GRAPHENE_RECT_INIT(0.f, 0.f, W * zoom, H * zoom);
	graphene_rect_t vis = GRAPHENE_RECT_INIT(vx, vy, vw, vh);
	flis_gpu_compose_render(snap, com.uniq->layers, W, H, &dst, &vis,
	                        GSK_SCALING_FILTER_NEAREST);
	GskRenderNode *node = gtk_snapshot_free_to_node(snap);
	if (node) gsk_render_node_unref(node);
}

int main(int argc, char **argv) {
	const char *path = argc > 1 ? argv[1] : "/workspace/flislrgb.fit";

	memset(&com, 0, sizeof(com));
	memset(&gui, 0, sizeof(gui));
	com.uniq = calloc(1, sizeof(single));
	com.uniq->next_item_id = 1;
	com.max_thread = g_get_num_processors();
	com.pref.swap_dir = g_strdup(g_get_tmp_dir());
	com.pref.gui.flis_tile_budget_mb = 48;   /* force eviction mid-render */
	com.script = TRUE;

	if (!gtk_init_check())
		fprintf(stderr, "note: no display, continuing\n");

	if (load_flis(path)) { fprintf(stderr, "load failed\n"); return 1; }
	const guint W = flis_canvas_rx(), H = flis_canvas_ry();

	gui.lo = 0; gui.hi = USHRT_MAX;
	for (int i = 0; i <= USHRT_MAX; i++)
		gui.remap_index[0][i] = (BYTE)(i >> 8);

	if (!flis_gpu_compose_compatible(com.uniq->layers)) {
		printf("not GPU-compatible, aborting\n");
		return 2;
	}

	GThread *inval = g_thread_new("invalidator", invalidator_thread, NULL);

	GMainContext *ctx = g_main_context_default();
	guint32 lcg = 7;
	for (int round = 0; round < 60; round++) {
		lcg = lcg * 1103515245u + 12345u;
		int kind = (lcg >> 16) % 5;
		switch (kind) {
		case 0:  /* fit-to-window render, whole canvas visible (mip 4/8) */
			render_once(W, H, 0.12f, 0.f, 0.f, (float)W, (float)H);
			break;
		case 1:  /* 100% zoom, small panning viewport (mip 1, eviction) */
		case 2: {
			lcg = lcg * 1103515245u + 12345u;
			float vx = (float)((lcg >> 8) % (W - 1600));
			lcg = lcg * 1103515245u + 12345u;
			float vy = (float)((lcg >> 8) % (H - 1000));
			render_once(W, H, 1.0f, vx, vy, 1600.f, 1000.f);
			break;
		}
		case 3:  /* stretch change → slot rebuild next render */
			gui.lo = (WORD)((lcg >> 4) & 0x0FFF);
			gui.hi = USHRT_MAX - (WORD)((lcg >> 14) & 0x0FFF);
			render_once(W, H, 0.25f, 0.f, 0.f, (float)W, (float)H);
			break;
		case 4: { /* invalidate a random layer mid-sequence */
			lcg = lcg * 1103515245u + 12345u;
			int id = 1 + (int)((lcg >> 9) % 4);
			flis_gpu_compose_invalidate_layer(id);
			render_once(W, H, 1.0f, 100.f, 100.f, 1600.f, 1000.f);
			break;
		}
		}
		/* run pending idles (background prefetch) like the GUI loop */
		int spins = 0;
		while (g_main_context_iteration(ctx, FALSE) && spins++ < 64)
			;
		if ((round % 10) == 0)
			printf("round %d ok\n", round);
	}
	g_run_invalidator = FALSE;
	g_thread_join(inval);
	printf("stress done, no crash\n");
	return 0;
}
