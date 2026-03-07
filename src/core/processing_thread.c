/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#include "core/siril.h"
#include "core/command_line_processor.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"

/*****************************************************************************
*    I N T E R N A L   T Y P E S   A N D   S T A T E
*****************************************************************************/

/* ── Queue and worker ───────────────────────────────────────────────────── */

static GThread  *worker_thread  = NULL;
static GMutex    queue_mutex;               /* guards everything below       */
static GCond     queue_cond;               /* worker wait + external waiters */
static GQueue    job_queue      = G_QUEUE_INIT;
static gboolean  worker_running = FALSE;

/* Thread-local flag: non-NULL only on the worker thread itself. */
static GPrivate  worker_tls = G_PRIVATE_INIT(NULL);

/* ── Cancellation and job-active ────────────────────────────────────────── */

/*
* cancel_flag — atomic int, 1 when cancellation has been requested.
* Cleared automatically at the start of each new job.
* Polled via get_thread_run() inside worker functions.
*/
static gint cancel_flag = 0;

/*
* job_active_flag — atomic int, 1 while a job is executing or the slot has
* been reserved via reserve_thread().  Distinct from cancel_flag: a job can
* be active without a cancel being requested, and vice versa after the job
* has finished but before the next one starts.
*
* Note on ordering with queue_cond broadcasts: job_active_flag is set to 0
* *before* broadcasting queue_cond after a job finishes (see worker loop),
* so waiters that re-read it after waking will always see the current value.
*/
static gint job_active_flag = 0;

/* ── Python reservation ──────────────────────────────────────────────────
*
* When python_reserved is TRUE, processing_submit_job blocks (on queue_cond)
* until python_releases_thread() clears it.
* Protected by queue_mutex.
*/
static gboolean python_reserved = FALSE;

/* ── Pending-result cache ────────────────────────────────────────────────
*
* When a blocking submission completes (script / python-command context),
* the worker deposits the job's return value here before signalling the
* submitter's job->cond.  The immediately-following waiting_for_thread()
* call consumes it without needing to wait again.
*
* Protected by queue_mutex.
*/
static gpointer pending_result       = NULL;
static gboolean pending_result_valid = FALSE;

/* ── Active-job reference ────────────────────────────────────────────────
*
* Points to the job currently being executed, or NULL between jobs.
* Used by waiting_for_thread() when no cached result is available
* (fire-and-forget GUI path where we need to wait on the job finishing).
*
* Protected by queue_mutex.
*/
static ProcessingJob *active_job = NULL;

/*****************************************************************************
*    I N T E R N A L   H E L P E R S
*****************************************************************************/

static void processing_job_free(ProcessingJob *job) {
	g_mutex_clear(&job->mutex);
	g_cond_clear(&job->cond);
	g_free(job);
}

/*
* Returns TRUE when the calling context requires the submission to block.
* False → fire-and-forget (GUI, arbitrary threads that are not scripts).
*/
static gboolean context_requires_wait(void) {
	return com.script || com.python_command;
}

/*****************************************************************************
*    W O R K E R   T H R E A D
*****************************************************************************/

static gpointer worker_thread_main(gpointer user_data G_GNUC_UNUSED) {
	/* Mark this thread so processing_in_worker_thread() returns TRUE here. */
	g_private_set(&worker_tls, GINT_TO_POINTER(1));

	for (;;) {

		/* ── Wait for a job ──────────────────────────────────────────── */
		g_mutex_lock(&queue_mutex);
		while (worker_running && g_queue_is_empty(&job_queue))
			g_cond_wait(&queue_cond, &queue_mutex);

		if (!worker_running) {
			g_mutex_unlock(&queue_mutex);
			break;
		}

		ProcessingJob *job = g_queue_pop_head(&job_queue);
		g_mutex_unlock(&queue_mutex);

		if (!job)
			break;

		if (job->type == PROCESSING_JOB_QUIT) {
			processing_job_free(job);
			break;
		}

		/* ── Prepare ─────────────────────────────────────────────────── */

		/* Clear any stale cancellation from the previous job. */
		g_atomic_int_set(&cancel_flag, 0);

		/* Mark slot occupied. */
		g_atomic_int_set(&job_active_flag, 1);

		/* Publish active_job so waiting_for_thread() can find us. */
		g_mutex_lock(&queue_mutex);
		active_job = job;
		g_mutex_unlock(&queue_mutex);

		/* ── Execute ─────────────────────────────────────────────────── */

		job->result = job->func(job->data);

		/* ── Publish result ──────────────────────────────────────────── */

		/*
		* Order matters:
		*   1. Clear job_active_flag *before* the broadcast so that any
		*      thread waking on queue_cond and re-reading the flag sees 0.
		*   2. Then lock queue_mutex, update shared result state, broadcast.
		*   3. Then signal the per-job cond (unblocks blocking submitters).
		*   4. Free if fire-and-forget.
		*/
		g_atomic_int_set(&job_active_flag, 0);

		g_mutex_lock(&queue_mutex);
		active_job = NULL;
		/* Only cache the result for blocking (script/python-command) jobs.
		* Fire-and-forget jobs have no caller waiting on pending_result, and
		* overwriting it would corrupt the result of a subsequent blocking
		* job whose waiting_for_thread() call consumes it.               */
		if (!job->fire_and_forget) {
			pending_result       = job->result;
			pending_result_valid = TRUE;
		}
		g_cond_broadcast(&queue_cond);
		g_mutex_unlock(&queue_mutex);

		/* Unblock any thread that called processing_wait_for_job() or the
		* blocking path inside processing_submit_job. */
		g_mutex_lock(&job->mutex);
		job->completed = TRUE;
		g_cond_broadcast(&job->cond);
		g_mutex_unlock(&job->mutex);

		if (job->fire_and_forget)
			processing_job_free(job);
	}

	return NULL;
}

/*****************************************************************************
*    L I F E C Y C L E
*****************************************************************************/

void processing_system_init(void) {
	g_mutex_init(&queue_mutex);
	g_cond_init(&queue_cond);
	worker_running = TRUE;
	worker_thread  = g_thread_new("processing", worker_thread_main, NULL);
}

void processing_system_shutdown(void) {
	if (!worker_running)
		return;

	/* Tell the worker loop to exit after draining. */
	g_mutex_lock(&queue_mutex);
	worker_running = FALSE;
	g_cond_broadcast(&queue_cond);
	g_mutex_unlock(&queue_mutex);

	/* Enqueue a sentinel so the worker wakes even if it is mid-wait.
	* (The loop also exits on worker_running=FALSE, but a QUIT job guarantees
	*  it drains any jobs pushed before the shutdown signal.) */
	ProcessingJob *quit = g_new0(ProcessingJob, 1);
	quit->type = PROCESSING_JOB_QUIT;
	g_mutex_init(&quit->mutex);
	g_cond_init(&quit->cond);

	g_mutex_lock(&queue_mutex);
	g_queue_push_tail(&job_queue, quit);
	g_cond_signal(&queue_cond);
	g_mutex_unlock(&queue_mutex);

	g_thread_join(worker_thread);
	worker_thread = NULL;

	g_mutex_clear(&queue_mutex);
	g_cond_clear(&queue_cond);
}

/*****************************************************************************
*    C O R E   S U B M I S S I O N
*****************************************************************************/

ProcessingJob *processing_submit_job(ProcessingFunc func, gpointer data) {
	g_return_val_if_fail(func != NULL, NULL);

	/* ── Re-entrant call from the worker thread ─────────────────────── */
	/*
	* If we are already on the worker thread a new submission would deadlock
	* (the worker is busy executing us; it can't pick up the new job).
	* Run synchronously instead, then cache the result.
	*/
	if (processing_in_worker_thread()) {
		gpointer result = func(data);
		g_mutex_lock(&queue_mutex);
		pending_result       = result;
		pending_result_valid = TRUE;
		g_mutex_unlock(&queue_mutex);
		return NULL;
	}

	/* ── Python-reservation gate ─────────────────────────────────────── */
	/*
	* If Python has reserved the slot, block until it is released — unless
	* we are on the GTK main thread, where blocking would prevent idles from
	* dispatching and deadlock the worker.  In that case return NULL and
	* let the caller handle the rejection (start_in_new_thread returns FALSE).
	*
	* Non-GTK background threads (script thread, Python thread, etc.) simply
	* wait here; this is safe because none of them run the GTK main loop.
	*/
	if (g_main_context_is_owner(g_main_context_default())) {
		g_mutex_lock(&queue_mutex);
		gboolean reserved = python_reserved;
		g_mutex_unlock(&queue_mutex);
		if (reserved) {
			siril_debug_print("processing_submit_job: Python reservation active, "
							"cannot submit from GTK main thread\n");
			return NULL;
		}
	} else {
		g_mutex_lock(&queue_mutex);
		while (python_reserved)
			g_cond_wait(&queue_cond, &queue_mutex);
		g_mutex_unlock(&queue_mutex);
	}

	/* ── Build and queue the job ─────────────────────────────────────── */

	gboolean must_wait = context_requires_wait();

	ProcessingJob *job = g_new0(ProcessingJob, 1);
	job->type           = PROCESSING_JOB_NORMAL;
	job->func           = func;
	job->data           = data;
	job->fire_and_forget = !must_wait;
	g_mutex_init(&job->mutex);
	g_cond_init(&job->cond);

	g_mutex_lock(&queue_mutex);
	g_queue_push_tail(&job_queue, job);
	g_cond_signal(&queue_cond);
	g_mutex_unlock(&queue_mutex);

	/* ── Blocking path (script / python-command context) ─────────────── */

	if (must_wait) {
		/*
		* Block until the worker signals job->cond.  By that point the
		* worker has already written pending_result (see worker loop), so
		* the next waiting_for_thread() call retrieves it instantly.
		*/
		g_mutex_lock(&job->mutex);
		while (!job->completed)
			g_cond_wait(&job->cond, &job->mutex);
		g_mutex_unlock(&job->mutex);

		processing_job_free(job);
		return NULL;   /* caller gets nothing to wait on — already done */
	}

	/* ── Fire-and-forget path (GUI / arbitrary threads) ─────────────── */

	return job;   /* caller may pass to processing_wait_for_job if needed */
}

void processing_wait_for_job(ProcessingJob *job) {
	if (!job)
		return;
	g_mutex_lock(&job->mutex);
	while (!job->completed)
		g_cond_wait(&job->cond, &job->mutex);
	g_mutex_unlock(&job->mutex);
}

/*****************************************************************************
*    C A N C E L L A T I O N   A N D   S T A T U S
*****************************************************************************/

void processing_request_cancel(void) {
	g_atomic_int_set(&cancel_flag, 1);
}

/*
* get_thread_run — FOR WORKER FUNCTIONS ONLY.
*
* Returns TRUE  → no cancellation requested; keep processing frames.
* Returns FALSE → cancellation was requested; abort the current job.
*
* Do not use this to test whether a job is currently active from outside a
* worker; use processing_is_job_active() for that purpose.
*/
gboolean get_thread_run(void) {
	return g_atomic_int_get(&cancel_flag) == 0;
}

gboolean processing_is_job_active(void) {
	return g_atomic_int_get(&job_active_flag) != 0;
}

gboolean processing_in_worker_thread(void) {
	return g_private_get(&worker_tls) != NULL;
}

/*****************************************************************************
*    P Y T H O N   R E S E R V A T I O N
*****************************************************************************/

int claim_thread_for_python(void) {
	/*
	* In the new design com.thread is never set (the worker thread is an
	* internal implementation detail).  We gate purely on job_active_flag and
	* the queue being empty.
	*/
	if (is_an_image_processing_dialog_opened()) {
		fprintf(stderr, "A Siril image processing dialog is open. "
				"It must be closed before Python can claim the processing thread.\n");
		return 2;
	}

	g_mutex_lock(&queue_mutex);

	if (python_reserved) {
		fprintf(stderr, "Processing thread is already reserved by Python.\n");
		g_mutex_unlock(&queue_mutex);
		return 1;
	}

	/*
	* Wait until any running job finishes and the pending queue drains.
	* The worker broadcasts queue_cond after each job completes (with
	* job_active_flag already cleared), so this wakes reliably.
	*/
	while (g_atomic_int_get(&job_active_flag) ||
		!g_queue_is_empty(&job_queue)) {
		g_cond_wait(&queue_cond, &queue_mutex);
	}

	python_reserved = TRUE;
	g_mutex_unlock(&queue_mutex);
	set_cursor_waiting(TRUE);
	return 0;
}

void python_releases_thread(void) {
	g_mutex_lock(&queue_mutex);
	python_reserved = FALSE;
	g_cond_broadcast(&queue_cond);   /* unblock any threads waiting in
										processing_submit_job              */
	g_mutex_unlock(&queue_mutex);
	set_cursor_waiting(FALSE);
}

/*****************************************************************************
*    C O M P A T I B I L I T Y   W R A P P E R S
*****************************************************************************/

gboolean start_in_new_thread(ProcessingFunc func, gpointer data) {
	/*
	* Re-entrant submissions (from within a running job) bypass all guards:
	* processing_submit_job will run func() synchronously in that case.
	*/
	if (!processing_in_worker_thread()) {

		/* ── Busy guard ──────────────────────────────────────────────────
		*
		* In the old design start_in_new_thread refused to start if
		* com.run_thread was set.  We restore that invariant: only one job
		* may be active at a time from the perspective of the external API.
		* This prevents rapid GUI triggers (e.g. multiple menu clicks) from
		* stacking up jobs that operate on the same gfit with arguments
		* captured at different points in time, producing dimension
		* mismatches and corrupted state in the worker (manifests as
		* "unknown interpolation type" from OpenCV for rotation jobs).
		*
		* We test job_active_flag rather than the queue length because a
		* queued-but-not-yet-started job is just as problematic as one
		* already running.
		*/
		if (processing_is_job_active()) {
			fprintf(stderr, "The processing thread is busy, stop it first.\n");
			return FALSE;
		}

		/* ── Python-reservation gate ─────────────────────────────────── */
		if (g_main_context_is_owner(g_main_context_default())) {
			g_mutex_lock(&queue_mutex);
			gboolean reserved = python_reserved;
			g_mutex_unlock(&queue_mutex);
			if (reserved) {
				fprintf(stderr, "The processing thread is reserved by Python; "
						"cannot start a job from the GTK main thread.\n");
				return FALSE;
			}
		}

		/*
		* Mark the slot active immediately, before the job reaches the
		* worker.  This closes the window between pushing to the queue and
		* the worker's own g_atomic_int_set, during which a second rapid
		* caller could pass the busy guard above.  The worker will set it
		* again before executing; the double-set is harmless.
		*/
		g_atomic_int_set(&job_active_flag, 1);
	}

	if (!com.headless)
		set_cursor_waiting(TRUE);

	if (!add_child((GPid) -2, INT_PROC_THREAD, "Siril processing thread"))
		siril_log_color_message(
			_("Warning: failed to add processing thread to child list\n"),
			"salmon");

	processing_submit_job(func, data);
	return TRUE;
}

gboolean start_in_reserved_thread(ProcessingFunc func, gpointer data) {
	/*
	* In the old design start_in_reserved_thread allowed submission even while
	* com.python_claims_thread was TRUE (it only checked com.thread).  In the
	* new design the single-worker queue serialises everything; we simply
	* delegate to start_in_new_thread.  The "reserved" quality is implicit:
	* there is always exactly one slot and it is always pre-allocated.
	*/
	return start_in_new_thread(func, data);
}

gpointer waiting_for_thread(void) {
	gpointer result = NULL;

	/* ── Guard: never block the GTK main thread ──────────────────────── */
	/*
	* If called from the GTK main thread a running fire-and-forget job may
	* itself need the GTK main thread to dispatch an idle via
	* execute_idle_and_wait_for_it.  Blocking here would deadlock.
	*
	* In practice, the CMD_NOTIFY_GFIT_MODIFIED path in execute_command only
	* calls waiting_for_thread() for commands that modify gfit synchronously
	* (no thread), so no job is in flight and the fast-path below returns 0
	* immediately anyway.  The guard is a safety net for any future caller.
	*/
	if (g_main_context_is_owner(g_main_context_default()))
		return GINT_TO_POINTER(0);

	g_mutex_lock(&queue_mutex);

	/* ── Fast path: result cached by a prior blocking submit ─────────── */
	/*
	* In script / python-command context, start_in_new_thread already blocked
	* until the job finished.  The worker deposited the result in
	* pending_result before signalling; we just retrieve it here.
	*/
	if (pending_result_valid) {
		result               = pending_result;
		pending_result       = NULL;
		pending_result_valid = FALSE;
		g_mutex_unlock(&queue_mutex);
		goto strip_and_return;
	}

	/* ── Slow path: fire-and-forget job may still be running ──────────── */
	/*
	* This path is taken when waiting_for_thread() is called from a non-GTK
	* background thread (e.g. the processcommand wait_for_completion path
	* with com.python_command=FALSE but wait_for_completion=TRUE) and the job
	* was submitted fire-and-forget.  We wait on queue_cond, which the worker
	* broadcasts after every job completion.
	*/
	while (active_job != NULL)
		g_cond_wait(&queue_cond, &queue_mutex);

	if (pending_result_valid) {
		result               = pending_result;
		pending_result       = NULL;
		pending_result_valid = FALSE;
	}

	g_mutex_unlock(&queue_mutex);

strip_and_return:
	/* Strip CMD_NOTIFY_GFIT_MODIFIED from the job's return value, matching
	* the behaviour of the old g_thread_join-based implementation.          */
	return GINT_TO_POINTER(GPOINTER_TO_INT(result) & ~CMD_NOTIFY_GFIT_MODIFIED);
}

void stop_processing_thread(void) {
	/*
	* In the old design this function joined com.thread, which was the source
	* of the intermittent deadlock (the worker would call it while blocked in
	* execute_idle_and_wait_for_it, causing a self-join via the GTK idle).
	*
	* In the new design the persistent worker thread is never stopped between
	* jobs.  This function's only responsibilities are:
	*   1. Signal cancellation so the running job aborts its frame loop.
	*   2. Remove the processing-thread entry from the child list.
	*   3. Reset the wait cursor.
	*
	* All three operations are non-blocking and safe to call from any thread,
	* including the worker itself (end_generic is called from within jobs in
	* the !run_idle / script-mode fallback path).
	*/
	processing_request_cancel();
	remove_child_from_children((GPid) -2);
	if (!com.headless)
		set_cursor_waiting(FALSE);
}

/* Compatibility shim — see header. */
void set_thread_run(gboolean b) {
	siril_debug_print("set_thread_run(%d)\n", b);
	if (!b)
		processing_request_cancel();
	/* b=TRUE is a no-op: job_active_flag is managed by the worker. */
}

gboolean reserve_thread(void) {
	/*
	* Atomically claim the active slot for a short synchronous operation.
	* Returns FALSE if already active (another job or reservation in progress).
	*/
	return g_atomic_int_compare_and_exchange(&job_active_flag, 0, 1);
}

void unreserve_thread(void) {
	g_atomic_int_set(&job_active_flag, 0);
	/* Broadcast so claim_thread_for_python waiters can re-evaluate. */
	g_mutex_lock(&queue_mutex);
	g_cond_broadcast(&queue_cond);
	g_mutex_unlock(&queue_mutex);
}

gboolean end_generic(gpointer arg G_GNUC_UNUSED) {
	stop_processing_thread();
	return FALSE;
}
