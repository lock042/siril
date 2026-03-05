/*
 * processing_thread.h — persistent single-worker-thread job queue
 *
 * Design intent
 * -------------
 * 1. One persistent worker thread runs for the lifetime of the application.
 *    Jobs are serialised through a FIFO queue; there is never more than one
 *    job executing at a time.
 *
 * 2. Submission is context-aware:
 *      - script thread (com.script)            → blocks until job completes
 *      - Python-command thread (com.python_command) → blocks until job completes
 *      - all other callers (GUI, arbitrary)    → fire-and-forget
 *    This makes the "one command must finish before the next can start"
 *    guarantee in execute_script automatic, without needing an external join.
 *
 * 3. Cancellation is a separate signal (processing_request_cancel / get_thread_run)
 *    and is distinct from "is a job active?" (processing_is_job_active).
 *    Call sites that previously used get_thread_run() for "is something running?"
 *    should migrate to processing_is_job_active(); get_thread_run() is kept for
 *    the cancel-poll semantics used inside worker functions.
 *
 * 4. stop_processing_thread() only does child-list cleanup and cursor reset.
 *    (Eventually cursor reset should probably go into the generic idle as well, but
 *    stop_processing_thread is a good catch-all to guarantee it gets reset.) It never
 *    joins or blocks.  Calling it from an idle callback (end_generic_*) or directly
 *    from a worker job function is always safe.
 */

#ifndef SIRIL_PROCESSING_THREAD_H
#define SIRIL_PROCESSING_THREAD_H

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif
/* --------------------------------------------------------------------------
 * Types
 * -------------------------------------------------------------------------- */

/*
 * Job function signature.  Identical to GThreadFunc so existing worker
 * functions (generic_image_worker, generic_sequence_worker, …) need no
 * signature change.  The return value is passed back through waiting_for_thread().
 */
typedef gpointer (*ProcessingFunc)(gpointer data);

/* Internal job type — used by the queue machinery. */
typedef enum {
    PROCESSING_JOB_NORMAL,
    PROCESSING_JOB_QUIT
} ProcessingJobType;

/*
 * Full struct definition is public so callers can hold a ProcessingJob on the
 * stack if needed, but treat all fields as private: use only the API below.
 *
 * Only non-NULL job handles are returned for fire-and-forget (GUI) submissions.
 * Blocking submissions (script / python-command) return NULL from
 * processing_submit_job / start_in_new_thread because the job is already
 * complete and freed by the time the call returns.
 */
typedef struct _ProcessingJob {
    ProcessingJobType  type;
    ProcessingFunc     func;
    gpointer           data;
    gpointer           result;

    GMutex             mutex;
    GCond              cond;
    gboolean           completed;
    gboolean           fire_and_forget;
} ProcessingJob;


/* --------------------------------------------------------------------------
 * Lifecycle  (call once at application start / stop)
 * -------------------------------------------------------------------------- */

void processing_system_init     (void);
void processing_system_shutdown (void);


/* --------------------------------------------------------------------------
 * Core submission
 * --------------------------------------------------------------------------
 *
 * If called from the worker thread itself (re-entrant / already_in_a_thread
 * pattern), func is invoked directly and synchronously; NULL is returned.
 *
 * If context_requires_wait() — i.e. com.script or com.python_command is set —
 * the call blocks until the job completes, the result is cached for the next
 * waiting_for_thread() call, and NULL is returned.
 *
 * Otherwise the job is queued and a non-NULL handle is returned; the caller
 * must not free it.  Use processing_wait_for_job() to wait on it explicitly.
 */
ProcessingJob *processing_submit_job (ProcessingFunc func, gpointer data);

/* Block until a specific fire-and-forget job has completed. */
void processing_wait_for_job (ProcessingJob *job);


/* --------------------------------------------------------------------------
 * Cancellation
 * --------------------------------------------------------------------------
 *
 * processing_request_cancel() sets the cancel flag polled by get_thread_run().
 * The flag is automatically cleared when the next job starts.
 *
 * get_thread_run()  — FOR USE INSIDE WORKER FUNCTIONS ONLY.
 *   Returns TRUE → keep running.  Returns FALSE → cancellation was requested.
 *   (This is the old "com.run_thread" semantics for the abort-check loop
 *    inside generic_sequence_worker etc.)
 *
 * processing_is_job_active()  — FOR EXTERNAL CALLERS.
 *   Returns TRUE while a job is executing or the slot has been reserved.
 *   Replaces the "is something running?" use of get_thread_run() in
 *   processcommand's wait loop and GUI button-sensitivity checks.
 */
void     processing_request_cancel (void);
gboolean get_thread_run (void);
gboolean processing_is_job_active (void);


/* --------------------------------------------------------------------------
 * Worker-thread identity
 * --------------------------------------------------------------------------
 *
 * Returns TRUE when called from the persistent worker thread.
 * Replaces the args->already_in_a_thread flag: processing_submit_job uses
 * this internally to run jobs synchronously when called re-entrantly from
 * within another job.  Callers that currently test already_in_a_thread to
 * gate idle dispatch may test this instead; the flag in the args struct can
 * remain for code that sets it explicitly.
 */
gboolean processing_in_worker_thread (void);


/* --------------------------------------------------------------------------
 * Python reservation
 * --------------------------------------------------------------------------
 *
 * claim_thread_for_python() drains the queue, waits for the active job to
 * finish, then gates all further submissions until python_releases_thread().
 *
 * Return values match the old API:
 *   0 → success (Python now owns the slot)
 *   1 → busy or already reserved
 *   2 → an image-processing dialog is open TODO: with this improved thread
 *       control mechanism do we stll need to worry about image processing
 *       dialogs being open when python claims the thread (or multiple image
 *       processing dialogs being open simultaneously)?
 */
int  claim_thread_for_python (void);
void python_releases_thread  (void);


/* --------------------------------------------------------------------------
 * Compatibility wrappers  (same signatures as before)
 * -------------------------------------------------------------------------- */

/*
 * start_in_new_thread / start_in_reserved_thread
 *
 * Both route through processing_submit_job.  The wait/no-wait decision is
 * driven by calling context:
 *   - script or python-command thread → blocks (ordering guarantee preserved)
 *   - any other thread, incl. GTK main (on_command_activate) → fire-and-forget
 *
 * In the new model start_in_reserved_thread is identical to start_in_new_thread
 * because the single-worker queue provides the "only one at a time" guarantee
 * implicitly.
 */
gboolean start_in_new_thread      (ProcessingFunc func, gpointer data);
gboolean start_in_reserved_thread (ProcessingFunc func, gpointer data);

/*
 * waiting_for_thread()
 *
 * Returns the job's return value with CMD_NOTIFY_GFIT_MODIFIED stripped, or 0.
 * TODO: once all image operations are migrated to generic_image_worker, there
 * will no longer be a need for the CMD_NOTIFY_GFIT_MODIFIED status.
 *
 * · Script / python-command context: start_in_new_thread already blocked, so
 *   this is a fast retrieval of the cached result — no additional waiting.
 * · Non-GTK background thread (e.g. processcommand wait_for_completion path):
 *   blocks until the active job finishes.
 * · GTK main thread: returns 0 immediately to avoid deadlocking the idle
 *   dispatcher (the CMD_NOTIFY_GFIT_MODIFIED path only ever calls this when no
 *   async job is in flight anyway).
 */
gpointer waiting_for_thread (void);

/*
 * stop_processing_thread()
 *
 * Requests cancellation of any running job, removes the processing-thread
 * entry from the child list, and resets the wait cursor.  Does NOT join or
 * block.  Safe to call from any thread, including the worker itself.
 */
void stop_processing_thread (void);

/*
 * set_thread_run(b)
 *
 * Compatibility shim:
 *   b = FALSE → processing_request_cancel()
 *   b = TRUE  → no-op (job_active is managed internally)
 */
void set_thread_run (gboolean b);

/*
 * reserve_thread / unreserve_thread
 *
 * Claim / release the "active slot" for short synchronous operations that need
 * to prevent concurrent submissions without running on the worker thread.
 * reserve_thread() is an atomic test-and-set; returns FALSE if already active.
 */
gboolean reserve_thread (void);
void unreserve_thread (void);

gboolean end_generic (gpointer arg);

#ifdef __cplusplus
}
#endif

#endif /* SIRIL_PROCESSING_THREAD_H */
