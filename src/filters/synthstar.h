#ifndef SRC_SYNTHSTAR_H_
#define SRC_SYNTHSTAR_H_

gpointer do_synthstar();
gpointer fix_saturated_stars();
void makegaussian(float *psf, int size, float fwhm, float lum, float xoffset, float yoffset);
void makemoffat(float *psf, const int size, const float fwhm, const float lum, const float xoff, const float yoff, const float beta);


#define SYNTHESIZE_GAUSSIAN 0
#define SYNTHESIZE_MOFFAT 1

#endif /* SRC_SYNTHSTAR_H_ */
