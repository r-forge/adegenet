/*
  Coded by Thibaut Jombart (tjombart@imperial.ac.uk), December 2010.
  Distributed with the adephylo package for the R software.
  Licence: GPL >=2.

  These functions are designed to recode genotypes given as binary integers into new integers which
  map them to unique bytes. One genotype of 8 binary SNPs is mapped uniquely (bijectively) to a
  value between 0 and 255. This is achieved by considering the genotype 'x' in the basis 2^0
  ... 2^7, and summing the values of the vector in this basis. That is, we use the function:

  {0,1}^8 |-> {0,...,255}
  x -> x_1 * 2^0 + ... + x_8 * 2^7 = \sum_i x_i * 2^(i-1)

*/


#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <R_ext/Utils.h>
#include "adesub.h"




/*
   ==========================
   === INTERNAL FUNCTIONS ===
   ==========================
*/


/* Maps one value from 0-255 to sequences of 8 binary values */
void byteToBinInt(unsigned char in, int *out){
	short int rest, i, temp;

	rest = (int)in;

	/* initialize all values to 0*/
	for(i=0;i<=7;i++)
		out[i]=0;

	for(i=7;i>=0;i--){
		temp = pow(2, i);
		if(rest >= temp) {
			out[i] = 1;
			rest = rest- temp;
			if(rest == 0) break;
		}
	}
}












/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/


/* 
   === MAP BINARY SNPS TO 1->256 SCALE ===
   - vecsnp: vector of integers (0/1)
   - vesize: length of vecsnp
   - res: vector of integers valued on 0:255
   - ressize: length of res
*/
void binIntToBytes(int *vecsnp, int *vecsize, unsigned char *vecres, int *ressize){
	/* declarations */
	int i, j, idres, *binBasis; /* must use dynamic allocation */

	/* allocate memory for local variables */
	vecintalloc(&binBasis, 8);

	/* define binary basis */
	for(i=1; i<=8; i++){
		binBasis[i] = pow(2, i-1);
	}

	/* set all values of vecres to 0 */
	for(i=0;i < *ressize;i++){
		vecres[i] = 0x00;
	}



	/* 
	   =============
	   MAIN FUNCTION
	   ============= 
	*/
	
	/* INDICES */
	/* i: idx of snp */
	/* j: idx of binBasis (1:8) */
	/* idres: idx in vector of results */

	idres = 0;
	j = 1;
	for(i=0;i< *vecsize;i++){
		vecres[idres] = vecres[idres] + (unsigned char)(binBasis[j] * vecsnp[i]);
		if(j == 8){
			idres++;
			j = 1;
		} else {
			j++;
		}
	}


	/* free memory */
	freeintvec(binBasis);

} /* end sptips */








/* Maps an array of values from 0-255 to sequences of 8 binary values */
/* Input are unsigned char (hexadecimal), outputs are integers */
void bytesToBinInt(unsigned char *vecbytes, int *vecsize, int *vecres){
	int i, j, idres=0, *temp; /* idres: index in vecres*/
	
	temp = (int *) calloc(8, sizeof(int));

	for(i=0;i<*vecsize;i++){
		byteToBinInt(vecbytes[i], temp);
		for(j=0;j<=7;j++){
			vecres[j+idres] = temp[j];
		}
		idres = idres + 8;
	}

	free(temp);
} /* end binIntToBytes*/








/* Maps an array of values from 0-255 to integers representing counts of alleles */
/* This is done by adding arrays of 0-1 for indiv with ploidy > 1*/
/* Input are unsigned char (hexadecimal), outputs are integers */
/* veclength is the length of one vector of bytes */
/* nbvec is the nb of input vectors*/
/* input 'vecbytes' is actually concatenated, ie of size veclength * nbvec */
void bytesToInt(unsigned char *vecbytes, int *veclength, int *nbvec, int *vecres){
	int i, j, k, idres=0, *temp; /* idres: index in vecres*/

	temp = (int *) calloc(8, sizeof(int));


	for(k=0;k<*nbvec;k++){ /* for all input vector */
		idres = 0;
		for(i=0;i<*veclength;i++){ /* for one input vector */
			byteToBinInt(vecbytes[i+ k* *veclength], temp); /* byte -> 8 int (0/1)*/
			for(j=0;j<=7;j++){ /* fill in the result*/
				vecres[j+idres] = vecres[j+idres] + temp[j];
			}
			idres = idres + 8;
		}
	}
	free(temp);
} /* end bytesToInt */








/* Simple test function */
/* Test: increases for a raw (unsigned char) vector */
void testRaw(unsigned char *a, int *n){
	int i;
	for(i=0; i<*n; i++){
		a[i] = (unsigned char)(i);
	}
}







/* Function to compute all dot products between individuals */
/* No centring, no scaling */
/* a: 2-dim array, dim n x p*/
/* naposi: 2-dim array, dim n x ...*/
/* nbna: array of nb of NAs for each individual*/
void dotProd(unsigned char **a, int *n, int *p, int **naposi, int *nbna){
	/* define variables, allocate memory */


	/* free memory */
	
}





/* TESTING in R */

/*
## test raw conversion
.C("testRaw", raw(256), 256L, PACKAGE="adegenet")

## test raw->int conversion
x <- sample(0:1,800,replace=TRUE)
toto <- .bin2raw(x)$snp
all(.C("bytesToBinInt", toto, length(toto), integer(length(toto)*8))[[3]]==x)

## test raw vec -> binary integers
.C("bytesToBinInt",as.raw(c(12,11)), 2L, integer(16), PACKAGE="adegenet")

## test several raw vec -> int (allele counts, any ploidy)
.C("bytesToInt",as.raw(c(12,11)), 1L, 2L, integer(8), PACKAGE="adegenet")
*/

