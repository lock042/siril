#ifndef SRC_GUI_OPEN_DIALOG_H_
#define SRC_GUI_OPEN_DIALOG_H_

/* cookies for the file chooser */
#define OD_NULL 	0
#define OD_FLAT 	1
#define OD_DARK 	2
#define OD_OFFSET 	3
#define OD_CWD 		4
#define OD_OPEN 	5
#define OD_CONVERT 	6
#define OD_BADPIXEL	7
#define OD_FLATLIB  8
#define OD_DARKLIB  9
#define OD_OFFSETLIB  10
#define OD_DISTOLIB  11

#define FITS_EXTENSIONS "*.fit;*.FIT;*.fits;*.FITS;*.fts;*.FTS;*.fit.fz;*.fits.fz;*.fts.fz;*.fit.gz;*.fits.gz;*.fts.gz"

void header_open_button_clicked();
void cwd_btton_clicked();

#endif /* SRC_GUI_OPEN_DIALOG_H_ */
