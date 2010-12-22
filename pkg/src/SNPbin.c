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
void binIntToBytes(int *vecsnp, int *vecsize, int *vecres, int *ressize){
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
		vecres[i] = 0;
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
		vecres[idres] = vecres[idres] + binBasis[j] * vecsnp[i];
		if(j == 8){
			idres = idres +1;
			j = 1;
		} else {
			j = j+1;
		}
	}
	
	
	/* free memory */
	freeintvec(binBasis);

} /* end sptips */










/* TESTING */
/*


*/
