#ifndef SRC_GUI_PHOTOMETRIC_CC_H_
#define SRC_GUI_PHOTOMETRIC_CC_H_

#include <stdio.h>
#include <glib.h>

void initialize_spcc_mirrors();
void check_gaia_archive_status();
gpointer gaia_check(gpointer user_data);
void get_spectrum_from_ui(xpsampled *spectrum, int chan);
void initialize_photometric_cc_dialog();
void initialize_spectrophotometric_cc_dialog();
int get_photometry_catalog_from_GUI();
gboolean populate_spcc_combos(gpointer user_data);
gpointer populate_spcc_combos_async(gpointer user_data);
void reset_spcc_filters();

/* NDE amend (nde_editors.c): reopen this dialog pre-filled from a recorded
 * PCC/SPCC step — the single-image record or a group joint record whose
 * nested photometric params it edits.  OK commits through the amend path;
 * closing cancels.  FALSE when the record's params cannot be read. */
gboolean pcc_open_amend(gint64 record_id);

#endif /* SRC_GUI_PHOTOMETRIC_CC_H_ */
