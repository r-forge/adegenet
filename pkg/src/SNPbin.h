#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <R_ext/Utils.h>
#include "adesub.h"

/* EXTERNAL */
void binIntToBytes(int *vecsnp, int *vecsize, unsigned char *vecres, int *ressize);
void bytesToBinInt(unsigned char *vecbytes, int *vecsize, int *vecres);
void bytesToInt(unsigned char *vecbytes, int *veclength, int *nbvec, int *vecres);

/* INTERNAL */
void byteToBinInt(unsigned char in, int out[8]);
void testRaw(unsigned char *a, int *n);
double snpbin_dotprod(struct snpbin *x, struct snpbin *y, double *mean, double*sd);
