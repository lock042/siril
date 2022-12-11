#ifndef DRAND48_REPLACEMENT
#define DRAND48_REPLACEMENT

#ifdef _WIN32

// Code to replicate drand48() from https://gist.github.com/mortennobel/8665258
#include <math.h>
#define RAND48_SEED_0   (0x330e)
#define RAND48_SEED_1 (0xabcd)
#define RAND48_SEED_2 (0x1234)
#define RAND48_MULT_0 (0xe66d)
#define RAND48_MULT_1 (0xdeec)
#define RAND48_MULT_2 (0x0005)
#define RAND48_ADD (0x000b)

unsigned short _rand48_seed[3] = {
        RAND48_SEED_0,
         RAND48_SEED_1,
         RAND48_SEED_2
};
unsigned short _rand48_mult[3] = {
         RAND48_MULT_0,
         RAND48_MULT_1,
         RAND48_MULT_2
 };
unsigned short _rand48_add = RAND48_ADD;

void
 _dorand48(unsigned short xseed[3]);

double erand48(unsigned short xseed[3]);

double drand48();

void srand48(long seed);
#endif

#endif
