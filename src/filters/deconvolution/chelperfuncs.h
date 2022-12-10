// Helpful wrapper functions to call Siril functions without causing clashes with some of the img_expr functionality

extern "C" void updateprogress(const char *text, double percent);
extern "C" void sirillog(const char* text);
extern "C" int is_thread_stopped();
extern "C" int updatenoise(float *array, int nx, int ny, int nchans, double *noise);
