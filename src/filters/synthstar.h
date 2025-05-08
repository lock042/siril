#ifndef SRC_SYNTHSTAR_H_
#define SRC_SYNTHSTAR_H_

gpointer do_synthstar(gpointer);
gpointer fix_saturated_stars(gpointer);
void makeairy(float *psf, const int size, const float lum, const float xoff, const float yoff, const float wavelength, const float aperture, const float focal_length, const float pixel_scale, const float obstruction);

void makegaussian(float *psf, int size, float fwhm, float lum, float xoffset, float yoffset, float ratio, float angle);
void makemoffat(float *psf, const int size, const float fwhm, const float lum, const float xoff, const float yoff, const float beta, float ratio, float angle);
void makedisc(float *psf, int size, float radius, float lum, float xoffset, float yoffset);
int starcount(psf_star **stars);
int free_psf_starstarstar(psf_star **stars);

#define SYNTHESIZE_GAUSSIAN 0
#define SYNTHESIZE_MOFFAT 1

#endif /* SRC_SYNTHSTAR_H_ */
