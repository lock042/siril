#ifndef SRC_CORE_NDE_REPLAY_INTERNAL_H_
#define SRC_CORE_NDE_REPLAY_INTERNAL_H_

/*
 * The vocabulary the four NDE replay implementation files share, and nothing
 * else.  nde_replay.h is the contract the rest of Siril sees; if a declaration
 * here becomes useful outside these files, that is the sign it belongs there
 * instead.  The files, and what each owns:
 *
 *   nde_chain.c    can these records be replayed, may this one be edited
 *   nde_replay.c   running records over pixels, and writing the result back
 *   nde_edit.c     amend / delete / reorder, and the cascade that follows
 *   nde_preview.c  the conductor thread, previews, and region-scoped replay
 *
 * Everything below is prefixed nde_ because these are ordinary global symbols
 * in a large C program; a bare commit_pixels() would be a collision waiting to
 * happen.  The prefix also tells a reader of nde_edit.c that commit_pixels is
 * not local to the file they are looking at.
 */

#include "core/siril.h"
#include "core/nde_history.h"
#include "core/nde_replay.h"

/* ---- nde_chain.c -------------------------------------------------------- */

/* @exclude_record_id != 0 builds the TRIAL chain for a pending delete — see
 * the definition.  nde_chain_build() is the ordinary (0) case. */
nde_chain *nde_chain_build_excluding(gint item_id, gint64 exclude_record_id);

/* First reason @rec's Tier-C script cannot be re-run (heap string), or NULL.
 * The chain build asks it to classify the record; the engine asks it again to
 * gate the actual re-run, so the answer cannot differ between the two. */
gchar *nde_tier_c_invalid_reason(const nde_record *rec);

/* ---- nde_replay.c: running records ------------------------------------- */

/* Apply chain members [@from, @upto) to @scratch (consumed), returning a new
 * fits.  @pos_x / @pos_y carry the layer position through geometry members, or
 * are NULL when the chain moves nothing. */
fits *nde_replay_apply_records(fits *scratch, const nde_chain *chain,
                               guint from, guint upto,
                               gint *pos_x, gint *pos_y, gchar **err);

/* Where a replay's layer starts: @restart_id 0 means the baseline, anything
 * else that record's output checkpoint.  FALSE when the chain carries no
 * geometry member, in which case the hooks are handed NULL and move nothing. */
gboolean nde_replay_start_offset(const nde_chain *chain, gint64 restart_id,
                                 gint *pos_x, gint *pos_y);

/* Replay a MASK item's chain up to member @upto, returning a fits carrying the
 * mask (a mask is a mono image with its own item id and its own history). */
fits *nde_mask_chain_replay(const nde_chain *chain, guint upto, gchar **err);

/* The phase-4 restart point: the state the chain's editable tail begins from. */
fits *nde_replay_chain_restart_state(const nde_chain *chain, gchar **err);

/* Install / remove @rec's pinned mask on a scratch fits for one op. */
gboolean nde_mask_pin_install(fits *scratch, const nde_record *rec, gchar **err);
void nde_mask_pin_clear(fits *scratch);

/* Set while an edit to a GEOMETRIC joint record cascades to its participants,
 * which have to re-anchor from the baseline — see the definition for why this
 * is not the general rule.  Owned here because nde_replay_start_offset reads
 * it; set by the cascade in nde_edit.c. */
void nde_replay_set_joint_reanchor(gboolean on);

/* ---- nde_replay.c: writing a value back -------------------------------- */

/* The fits an item's pixels live in: gfit for a plain image, the layer's own
 * fit for FLIS (the same pointer when it is the active layer).  NULL when the
 * item names nothing live. */
fits *nde_edit_target_fits(gint item_id);

/* Quiesced commit window.  Take the lock, swap, release: nde_commit_pixels
 * exchanges @target's and @result's pixels wholesale, so @result comes back
 * holding the old ones and the swap can be undone by repeating it. */
gboolean nde_commit_lock(fits *fit);
void nde_commit_pixels(fits *target, fits *result);
void nde_commit_unlock(fits *fit, gboolean quiesced);

/* Move @item_id's layer to where the replay left it. */
void nde_commit_layer_offset(gint item_id, gint pos_x, gint pos_y);

/* Carry the pre-edit metadata (keywords, header text, FITS HISTORY) back onto
 * @target after a swap, keeping a WCS the replay produced. */
void nde_commit_restore_metadata(fits *target, fits *old);

/* What every successful history edit does once its own work is committed.
 * @target NULL when the edit changed no live pixels of its own. */
void nde_edit_finish(fits *target, const char *done_msg);

/* ---- nde_edit.c --------------------------------------------------------- */

/* Shared core of amend (@new_params != NULL) and delete (NULL).  Conductor
 * context: the caller holds SLOT_REPLAY. */
gboolean nde_edit_execute(gint64 record_id, const gchar *new_params, gchar **err);

/* The state to restart an edited chain from, at or before member @e. */
fits *nde_edit_restart_state(const nde_chain *chain, guint e,
                             gint64 boundary_pre_id,
                             guint *start_idx, gchar **err);

/* The pixels a replay produced, and the log change that describes them.
 * @result is consumed either way. */
typedef struct {
	gint      item_id;
	fits     *target;        /* NULL when @retained */
	fits     *result;        /* consumed */
	gboolean  retained;      /* nde_item_is_retained_input(@item_id) */
	gboolean  carry_offset;  /* the replay moved the layer */
	gint      pos_x, pos_y;
} nde_commit_ctx;

/* Commit a replayed value: pixels first, then the log, and put the pixels back
 * if the log refuses.  @log_commit is what the caller's kind of edit does to
 * the history, or NULL when the log already says what these pixels are. */
gboolean nde_commit_replayed(nde_commit_ctx *c,
                             gboolean (*log_commit)(gpointer, gchar **),
                             gpointer log_user, gchar **err);

/* Everything derived from @item_id's pixels, rebuilt in dependency order.
 * Runs AFTER the pixel and log commits — see the definition. */
void nde_cascade_from(gint item_id, gint64 unchanged_upto,
                      GArray *joint_targets, gboolean joint_reanchor);

/* The participants a joint record at or after @from_record_id derives factors
 * for, collected BEFORE the log is touched. */
GArray *nde_edit_joint_targets(gint item_id, gint64 from_record_id,
                               gboolean include_self);

/* The record before @record_id among @item_id's own, or 0.  Measured before
 * the log is touched, because a delete removes the anchor. */
gint64 nde_log_predecessor(gint item_id, gint64 record_id);

/* ---- nde_preview.c ------------------------------------------------------ */

/* TRUE while a preview holds pre-K pixels in a target fits — every history
 * edit must refuse then, because a commit would swap pixels the preview's
 * restore path is about to overwrite. */
gboolean nde_amend_preview_installed(void);

#endif /* SRC_CORE_NDE_REPLAY_INTERNAL_H_ */
