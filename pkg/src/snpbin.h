#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>


#define NEARZERO 0.0000000001 */
#define TRUE 1
#define FALSE 0

typedef short bool;


/*
   =========================
   === CLASS DEFINITIONS ===
   =========================
*/

struct snpbin{
	unsigned char *bytevec;
	int *byteveclength, *bytevecnb, *nloc, *nanb, *naposi, *ploidy; /* all but naposi have length 1 */
};




struct snpbin makesnpbin(unsigned char *bytevec, int *byteveclength, int *bytevecnb, int *nloc, int *nanb, int *naposi, int *ploidy) {
	struct snpbin out;
	int i;

	out.bytevec = bytevec;
	out.byteveclength = byteveclength;
	out.bytevecnb = bytevecnb;
	out.nloc = nloc;
	out.nanb = nanb;
	/* need to decrease the indices of NAs by 1, e.g. [1-10]->[0-9] */
	out.naposi = naposi;
	if(*nanb > 0){
		for(i=0;i< *nanb; i++){
			out.naposi[i] = out.naposi[i] - 1;
		}
	}
	out.ploidy = ploidy;
	return out;
};



struct genlightC{
	struct snpbin *x;
	int *nind;
};



/*
   ===========================
   === AUXILIARY FUNCTIONS ===
   ===========================
*/


/* Maps one byte from 0-255 to sequences of 8 (binary) integers values */
void byteToBinInt(unsigned char in, int *out);


/* Maps an array of values from 0-255 to sequences of 8 binary values */
/* Input are unsigned char (hexadecimal), outputs are integers */
void bytesToBinInt(unsigned char *vecbytes, int *vecsize, int *vecres);






/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

void bytesToInt(unsigned char *vecbytes, int *veclength, int *nbvec, int *vecres);
void binIntToBytes(int *vecsnp, int *vecsize, unsigned char *vecres, int *ressize);





/*
   =====================
   === CLASS METHODS ===
   =====================
*/

int nLoc(struct snpbin *x);
int ploidy(struct snpbin *x);
void snpbin2intvec(struct snpbin *x, int *out);
void snpbin2freq(struct snpbin *x, double *out);
void printsnpbin(struct snpbin *x);
short int snpbin_isna(struct snpbin *x, int i);
double snpbin_dotprod(struct snpbin *x, struct snpbin *y, double *mean, double *sd);
struct genlightC genlightTogenlightC(unsigned char *gen, int *nbvecperind, int *byteveclength, int *nbnaperind, int *naposi, int *nind, int *nloc, int *ploidy);








/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/


void testRaw(unsigned char *a, int *n);
