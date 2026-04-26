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

/* Internal job type — used by the queue machinery. */
typedef enum {
    PROCESSING_JOB_NORMAL,
    PROCESSING_JOB_QUIT
} ProcessingJobType;

/*
 * Job descriptor.  For blocking jobs (script / python-command context) the
 * mutex and cond are initialised so the submitter can wait on completion.
 * For fire-and-forget jobs they are left uninitialised — no one ever waits
 * on them, and the worker frees the struct immediately after execution.
 */
typedef struct _ProcessingJob {
    ProcessingJobType  type;
    ProcessingFunc     func;
    gpointer           data;
    gpointer           result;

    gboolean           fire_and_forget;

    /* Initialised only when fire_and_forget == FALSE. */
    GMutex             mutex;
    GCond              cond;
    gboolean           completed;
} ProcessingJob;


/* ── Queue and worker ───────────────────────────────────────────────────── */

static GThread  *worker_thread  = NULL;
static GMutex    queue_mutex;               /* guards everything below       */
static GCond     queue_cond;               /* worker wait + external waiters */
static GQueue    job_queue      = G_QUEUE_INIT;
static gboolean  worker_running = FALSE;

/*
 * Guard against double-initialization.  Without this, calling init() twice
 * would create two worker threads with no way to stop the first, leaking
 * its GThread, queue_mutex, and queue_cond.
 * Must use gint for g_atomic_int_* compatibility.
 */
static gint      system_initialized = 0;

/* Thread-local flag: non-NULL only on the worker thread itself. */
static GPrivate  worker_tls = G_PRIVATE_INIT(NULL);

/* ── Cancellation and job-active ────────────────────────────────────────── */

/*
* cancel_flag — atomic int, 1 when cancellation has been requested.
* Cleared automatically at the start of each new job.
* Polled via processing_should_continue() inside worker functions.
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
* pending_result_valid is kept separate from pending_result because NULL
* (0) is a legitimate return value from a ProcessingFunc — it cannot serve
* as a sentinel to indicate "no result cached".
*
* Protected by queue_mutex.
*/
static gpointer pending_result       = NULL;
static gboolean pending_result_valid = FALSE;

/* ── Active-job flag ─────────────────────────────────────────────────────
*
* TRUE while the worker is executing a job.  Used by waiting_for_thread()
* to wait on fire-and-forget job completion via queue_cond.
*
* Protected by queue_mutex.
*/
static gboolean active_job_running = FALSE;

/* ── Pseudo-PID for child list tracking ────────────────────────────────── */
#define PROCESSING_THREAD_PSEUDO_PID ((GPid) -2)

/*****************************************************************************
*    I N T E R N A L   H E L P E R S
*****************************************************************************/

static void processing_job_free(ProcessingJob *job) {
    if (!job->fire_and_forget) {
        g_mutex_clear(&job->mutex);
        g_cond_clear(&job->cond);
    }
    g_free(job);
}

/*
* Returns TRUE when the calling context requires the submission to block.
* False -> fire-and-forget (GUI, arbitrary threads that are not scripts).
*
* NOTE: This captures transient global state (com.script, com.python_command)
* at submission time.  These flags should remain stable for the duration of
* any single command's execution; callers must ensure they are not modified
* concurrently with a call to processing_submit_job or start_in_new_thread.
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

        /* Publish active-job flag so waiting_for_thread() can find us. */
        g_mutex_lock(&queue_mutex);
        active_job_running = TRUE;
        g_mutex_unlock(&queue_mutex);

        /* ── Execute ─────────────────────────────────────────────────── */

        job->result = job->func(job->data);

        /* ── Publish result ──────────────────────────────────────────── */

        /*
        * Order matters:
        *   1. Clear job_active_flag *before* the broadcast so that any
        *      thread waking on queue_cond and re-reading the flag sees 0.
        *   2. Then lock queue_mutex, update shared result state, broadcast.
        *   3. For blocking jobs, signal the per-job cond.
        *   4. Free fire-and-forget jobs.
        */
        g_atomic_int_set(&job_active_flag, 0);

        g_mutex_lock(&queue_mutex);
        active_job_running = FALSE;
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

        if (job->fire_and_forget) {
            /* No waiter exists; free immediately. */
            processing_job_free(job);
        } else {
            /* Unblock the submitter waiting in processing_submit_job. */
            g_mutex_lock(&job->mutex);
            job->completed = TRUE;
            g_cond_broadcast(&job->cond);
            g_mutex_unlock(&job->mutex);
            /* Submitter owns and frees the job. */
        }
    }

    return NULL;
}

/*****************************************************************************
*    L I F E C Y C L E
*****************************************************************************/

void processing_system_init(void) {
    /*
    * Guard against double-initialization.  Without this, calling init()
    * twice would create two worker threads with no way to stop the first,
    * leaking its GThread, queue_mutex, and queue_cond.
    */
    if (!g_atomic_int_compare_and_exchange(&system_initialized, 0, 1))
        return;

    g_mutex_init(&queue_mutex);
    g_cond_init(&queue_cond);
    worker_running = TRUE;
    worker_thread  = g_thread_new("processing", worker_thread_main, NULL);
}

void processing_system_shutdown(void) {
    if (!worker_running)
        return;

    /* Cancel any running job so it can exit its frame loop quickly. */
    processing_request_cancel();

    if (com.headless) {
        /*
         * CLI / headless mode: execute_idle_and_wait_for_it is a no-op, so
         * the worker can never be blocked waiting for an idle callback that
         * needs the GTK main loop.  A clean join is safe here.
         *
         * Push a QUIT sentinel so the worker exits its loop cleanly rather
         * than reading worker_running=FALSE before the queue is drained.
         */
        ProcessingJob *quit = g_new0(ProcessingJob, 1);
        quit->type = PROCESSING_JOB_QUIT;
        quit->fire_and_forget = TRUE;

        g_mutex_lock(&queue_mutex);
        g_queue_push_tail(&job_queue, quit);
        g_cond_signal(&queue_cond);
        g_mutex_unlock(&queue_mutex);

        g_thread_join(worker_thread);
        worker_thread  = NULL;
        worker_running = FALSE;

        g_mutex_clear(&queue_mutex);
        g_cond_clear(&queue_cond);
    } else {
        /*
         * GUI mode: gtk_main_quit() was already called before reaching here,
         * so the GTK main loop is no longer running.  If the worker is blocked
         * inside execute_idle_and_wait_for_it (waiting for an idle that will
         * never be dispatched), g_thread_join() would deadlock.
         *
         * Application close is unconditional — just wake the worker in case
         * it is idle-waiting on queue_cond, then abandon the thread.  The
         * process is about to exit; the OS reclaims all resources.
         */
        g_mutex_lock(&queue_mutex);
        worker_running = FALSE;
        g_cond_signal(&queue_cond);
        g_mutex_unlock(&queue_mutex);

        /* Do NOT join and do NOT clear queue_mutex/queue_cond: the worker
         * thread may still be running and using them.  The OS handles
         * cleanup on process exit. */
        worker_thread = NULL;
    }

    /* Allow re-initialization if init is called again after shutdown. */
    g_atomic_int_set(&system_initialized, 0);
}

/*****************************************************************************
*    C O R E   S U B M I S S I O N
*****************************************************************************/

/*
* processing_submit_job — submit a job to the processing queue.
*
* Returns TRUE if the job was successfully queued (or executed synchronously
* in the re-entrant case), FALSE if submission was rejected (e.g. Python
* reservation active when called from the GTK main thread).
*/
gboolean processing_submit_job(ProcessingFunc func, gpointer data) {
    g_return_val_if_fail(func != NULL, FALSE);

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
        return TRUE;
    }

    /* ── Python-reservation gate ─────────────────────────────────────── */
    /*
    * If Python has reserved the slot, block until it is released — unless
    * we are on the GTK main thread, where blocking would prevent idles from
    * dispatching and deadlock the worker.  In that case return FALSE and
    * let the caller handle the rejection.
    *
    * Non-GTK background threads (script thread, Python thread, etc.) simply
    * wait here; this is safe because none of them run the GTK main loop.
    */
    if (g_main_context_is_owner(g_main_context_default())) {
        g_mutex_lock(&queue_mutex);
        gboolean reserved = python_reserved;
        g_mutex_unlock(&queue_mutex);
        if (reserved) {
            siril_log_color_message(
                _("Processing thread is reserved by Python; "
                  "cannot submit from GTK main thread.\n"),
                "salmon");
            return FALSE;
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

    if (must_wait) {
        g_mutex_init(&job->mutex);
        g_cond_init(&job->cond);
    }

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
    }

    /* Fire-and-forget path: job is now owned by the worker. */
    return TRUE;
}

/*****************************************************************************
*    C A N C E L L A T I O N   A N D   S T A T U S
*****************************************************************************/

void processing_request_cancel(void) {
    g_atomic_int_set(&cancel_flag, 1);
}

/*
* processing_should_continue — FOR WORKER FUNCTIONS ONLY.
*
* Returns TRUE  -> no cancellation requested; keep processing frames.
* Returns FALSE -> cancellation was requested; abort the current job.
*
* Do not use this to test whether a job is currently active from outside a
* worker; use processing_is_job_active() for that purpose.
*/
gboolean processing_should_continue(void) {
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
        siril_log_color_message(
            _("A Siril image processing dialog is open. "
              "It must be closed before Python can claim the processing thread.\n"),
            "salmon");
        return 2;
    }

    g_mutex_lock(&queue_mutex);

    if (python_reserved) {
        siril_log_color_message(
            _("Processing thread is already reserved by Python.\n"),
            "salmon");
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

gboolean processing_is_reserved_for_python(void) {
    g_mutex_lock(&queue_mutex);
    gboolean reserved = python_reserved;
    g_mutex_unlock(&queue_mutex);
    return reserved;
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
            siril_log_color_message(
                _("The processing thread is busy, stop it first.\n"),
                "salmon");
            return FALSE;
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

    if (!add_child(PROCESSING_THREAD_PSEUDO_PID, INT_PROC_THREAD,
                   "Siril processing thread"))
        siril_log_color_message(
            _("Warning: failed to add processing thread to child list\n"),
            "salmon");

    /*
    * Submit the job.  processing_submit_job handles the Python-reservation
    * gate internally — we do NOT duplicate that check here to avoid a
    * TOCTOU race (checking python_reserved, releasing the mutex, then
    * having processing_submit_job re-check under the same mutex).
    *
    * If submission fails, we must clean up the state we just acquired so
    * that the processing system is not left in a permanently broken state.
    */
    if (!processing_submit_job(func, data)) {
        if (!processing_in_worker_thread())
            g_atomic_int_set(&job_active_flag, 0);
        if (!com.headless)
            set_cursor_waiting(FALSE);
        remove_child_from_children(PROCESSING_THREAD_PSEUDO_PID);
        return FALSE;
    }

    return TRUE;
}

gboolean start_in_reserved_thread(ProcessingFunc func, gpointer data) {
    /*
    * The caller has already claimed the active slot via reserve_thread(), which
    * sets job_active_flag = 1.  We must NOT delegate to start_in_new_thread()
    * because that function's busy guard tests processing_is_job_active() and
    * would immediately reject the submission with "The processing thread is
    * busy" — even though the reservation was made legitimately.
    *
    * Instead we reproduce the submission path from start_in_new_thread while
    * skipping only that guard.  The Python-reservation gate is handled
    * inside processing_submit_job, avoiding a TOCTOU race.
    */

    if (!com.headless)
        set_cursor_waiting(TRUE);

    if (!add_child(PROCESSING_THREAD_PSEUDO_PID, INT_PROC_THREAD,
                   "Siril processing thread"))
        siril_log_color_message(
            _("Warning: failed to add processing thread to child list\n"),
            "salmon");

    /* job_active_flag is already 1 from reserve_thread(); no need to set it
     * again.  The worker sets it once more before executing, which is harmless. */
    if (!processing_submit_job(func, data)) {
        if (!com.headless)
            set_cursor_waiting(FALSE);
        remove_child_from_children(PROCESSING_THREAD_PSEUDO_PID);
        return FALSE;
    }

    return TRUE;
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
    *
    * IMPORTANT: If a new job starts between the fast-path check above and
    * entering this loop, we may end up waiting for the wrong job.  The
    * pending_result_valid re-check after the loop mitigates this: if a
    * blocking job completed while we were waiting, its result is returned.
    *
    * NOTE: for fire-and-forget jobs pending_result is never set by the worker
    * (see worker loop comment — overwriting it would corrupt a subsequent
    * blocking job's result).  Therefore the result of the awaited job is NOT
    * recoverable here and this function will return 0.  All command-path
    * callers that need the return value must ensure they submit via a blocking
    * context (com.script = TRUE or com.python_command = TRUE) so that the
    * fast path above is taken instead.
    */
    while (active_job_running)
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
    remove_child_from_children(PROCESSING_THREAD_PSEUDO_PID);
    if (!com.headless)
        set_cursor_waiting(FALSE);
}

gboolean reserve_thread(void) {
    /*
    * Atomically claim the active slot for a short synchronous operation.
    * Returns FALSE if already active (another job or reservation in progress).
    *
    * Also clears cancel_flag, mirroring what the worker does at the start of
    * each queued job.  This matters for synchronous callers (compositing
    * alignment, etc.) that call generic_sequence_worker directly with
    * already_in_a_thread=TRUE: those use processing_should_continue() to poll
    * for abort, so a stale cancel_flag left by a previous
    * stop_processing_thread() call would cause them to terminate immediately
    * on the first frame.
    */
    if (!g_atomic_int_compare_and_exchange(&job_active_flag, 0, 1))
        return FALSE;
    g_atomic_int_set(&cancel_flag, 0);
    return TRUE;
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

/*
 * cancel_and_wait_for_preview  — FOR GUI PREVIEW CANCEL PATHS ONLY.
 *
 * Sets the cancel flag then spins, pumping the GTK main loop, until the
 * active job finishes.  Safe to call from the GTK main thread because
 * preview jobs use siril_add_idle (asynchronous) for their cleanup and
 * never call execute_idle_and_wait_for_it, so they do not need the main
 * thread to complete.
 *
 * Do NOT use for jobs that call execute_idle_and_wait_for_it — those
 * would deadlock if the main thread is blocked here.
 *
 * Memory ordering: g_atomic_int_get() provides acquire semantics, which is
 * sufficient for the boolean job_active_flag.  The worker sets this flag to
 * 0 (release) before broadcasting queue_cond, establishing a happens-before
 * relationship that ensures the worker's side effects are visible once this
 * loop exits.
 */
void cancel_and_wait_for_preview(void) {
    if (!processing_is_job_active())
        return;
    processing_request_cancel();
    /* Pump the GTK main loop (non-blocking iterations) while the worker
     * finishes.  The worker does not need the main thread; we add a brief
     * sleep to avoid spinning 100% CPU. */
    while (processing_is_job_active()) {
        if (g_main_context_pending(g_main_context_default()))
            g_main_context_iteration(g_main_context_default(), FALSE);
        else
            g_usleep(1000); /* 1 ms */
    }
}

/*
 * start_and_wait_from_main_thread  — submit a job and block until done.
 *
 * Intended for the rare case where GUI code (running on the GTK main thread)
 * must submit a synchronous job and read back its return value.  The job
 * MUST NOT call execute_idle_and_wait_for_it internally, or a deadlock will
 * result.  (For generic_sequence_worker called from GUI context, the worker
 * uses siril_add_idle — asynchronous — so blocking here is safe.)
 *
 * Returns FALSE if the job could not be started (same semantics as
 * start_in_new_thread).
 *
 * Memory ordering: see cancel_and_wait_for_preview for the acquire/release
 * relationship between job_active_flag and queue_cond.
 */
gboolean start_and_wait_from_main_thread(ProcessingFunc func, gpointer data) {
    if (!start_in_new_thread(func, data))
        return FALSE;
    /* Pump the GTK main loop while the worker finishes. */
    while (processing_is_job_active()) {
        if (g_main_context_pending(g_main_context_default()))
            g_main_context_iteration(g_main_context_default(), FALSE);
        else
            g_usleep(1000);
    }
    return TRUE;
}
