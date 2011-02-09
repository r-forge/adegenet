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


  # Function named as 'SNPbin...' or 'GL...' are to be called directly from R.
  # The structure 'snpbin' is a C representation of the class 'SNPbin'.
  # Function named as 'snpbin...' are made to be called internally.
*/


#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <R_ext/Utils.h>
#include "adesub.h"

#define NEARZERO 0.0000000001


/*
   ========================
   === CLASS DEFINITION ===
   ========================
*/

/* 'bytevecnb' arrays of bytes concatenated into a single array */
/* of dim 'byteveclength' x 'bytevecnb' */
/* nloc is the number of SNPs - used for recoding to integers */
/* naposi indicates the positions of NAs */
/* nanb is the length of naposi */

struct snpbin{
	unsigned char *bytevec;
	int *byteveclength, *bytevecnb, *nloc, *nanb, *naposi; /* all but naposi have length 1 */
};




struct snpbin makesnpbin(unsigned char *bytevec, int *byteveclength, int *bytevecnb, int *nloc, int *nanb, int *naposi) {
	struct snpbin out;

	out.bytevec = bytevec;
	out.byteveclength = byteveclength;
	out.bytevecnb = bytevecnb;
	out.nloc = nloc;
	out.nanb = nanb;
	out.naposi = naposi;
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










/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/



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
				vecres[j+idres] += temp[j];
			}
			idres = idres + 8;
		}
	}
	free(temp);
} /* end bytesToInt */





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

} /* end binIntToBytes */








/*
   =====================
   === CLASS METHODS ===
   =====================
*/

int nLoc(struct snpbin *x){
	return *(x->nloc);
}



/* transform a snpbin into a vector of integers */
void snpbin2intvec(struct snpbin *x, int *out){
	bytesToInt(x->bytevec, x->byteveclength, x->bytevecnb, out);
/*reminders: 
-void bytesToInt(unsigned char *vecbytes, int *veclength, int *nbvec, int *vecres){
- snpbin: unsigned char *bytevec; int *byteveclength, *bytevecnb, *nloc, *nanb, *naposi; */
}


/* print a snpbin object - used for debugging */
void printsnpbin(struct snpbin *x){
	int i, *temp;
	temp = (int *) calloc(nLoc(x), sizeof(int));
	snpbin2intvec(x, temp);


	for(i=0;i< *(x->byteveclength);i++){
		printf("%i ", (int) (x->bytevec)[i]);
	}
	printf("   ");
	for(i=0;i<nLoc(x);i++){
		printf("%i ", temp[i]);
	}

	free(temp);
}



short int snpbin_isna(struct snpbin *x, int i){
	int j = 0;
	if(*(x->nanb) < 1 || i > nLoc(x)) return 0;

	while(j < *(x->nanb)){
		if( i == (x->naposi)[j]) return 1;
		j++;
	}

	return 0;
}





/* Function to compute one dot products between two individuals */
/* centring and scaling is always used */
/* but need to pass vectors of 0 and 1*/
double snpbin_dotprod(struct snpbin *x, struct snpbin *y, double *mean, double *sd){
	/* define variables, allocate memory */
	int P = nLoc(x), i, *vecx, *vecy;
	short int isna;
	double res = 0.0;
	vecx = (int *) calloc(P, sizeof(int));
	vecy = (int *) calloc(P, sizeof(int));

	/* conversion to integers */
	snpbin2intvec(x, vecx);
	snpbin2intvec(y, vecy);

	/* printf("\nvector x: \n"); */
	/* for(i=0;i<P;i++){ */
	/* 	printf("%i", vecx[i]); */
	/* } */

	/* printf("\nvector y: \n"); */
	/* for(i=0;i<P;i++){ */
	/* 	printf("%i", vecy[i]); */
	/* } */

	/* compute dot product */
	int count=0;
	for(i=0;i<P;i++){
		if(snpbin_isna(x,i) == 0 && snpbin_isna(y,i) == 0){
			/* res += ((vecx[i]-mean[i])/sd[i]) * ((vecy[i]-mean[i])/sd[i]); */
			res += ((vecx[i]-mean[i])/sd[i]) * ((vecy[i]-mean[i])/sd[i]);
			/* printf("\ntemp value of increment: %f", ((vecx[i]-mean[i])/sd[i]) * ((vecy[i]-mean[i])/sd[i])); */
			/* printf("\ntemp value of result: %f", res); */
		}
	}

	/* free memory */
	free(vecx);
	free(vecy);

	return res;
}




/* Function to convert a 'genlight' object (R side) into an array of 'snpbin' (C side) */
/* Each component of the genlight is concatenated into a single vector */
/* and then used to create different 'snpbin' on the C side */
struct genlightC genlightTogenlightC(unsigned char *gen, int *nbvecperind, int *byteveclength, int *nbnaperind, int *naposi, int *nind, int *nloc){
	/* declare variables and allocate memory */
	int i, j, idxByteVec=0, idxNAVec=0;
	struct genlightC out;
	out.x = (struct snpbin *) calloc(*nind, sizeof(struct snpbin));

	/* create the list of snpbin */
	/* printf("\n nind: %d\n", *nind); */
	for(i=0; i < *nind; i++){
		out.x[i] = makesnpbin(&gen[idxByteVec], byteveclength, &nbvecperind[i], nloc, &nbnaperind[i], &naposi[idxNAVec]);
		idxByteVec += *byteveclength * nbvecperind[i]; /* update index in byte array */
		idxNAVec +=  nbnaperind[i]; /* update index in byte array */
		/* printf("\nimported genotype %i: ", i+1); */
		/* printsnpbin(&out.x[i]); */
	}
	
	/* printf("step 3"); */

	out.nind = nind;

	/* printf("step 4"); */
	return out;
}




/* Function to compute all dot products between individuals */
/* centring and scaling is always used */
/* but need to pass vectors of 0 and 1*/
void GLdotProd(unsigned char *gen, int *nbvecperind, int *byteveclength, int *nbnaperind, int *naposi, int *nind, int *nloc, double *mean, double *sd, double *res){
	struct genlightC dat;
	int i, j, k=0;

	/* Check variance vector: do not divide by 0 */
	for(i=0;i< *nloc;i++){
		if(sd[i] < NEARZERO){
			sd[i] = 1;
		}
	}

	dat = genlightTogenlightC(gen, nbvecperind, byteveclength, nbnaperind, naposi, nind, nloc);

	/* Lower triangle - without the diagonal */
	for(i=0; i< (*nind-1); i++){
		for(j=i+1; j< *nind; j++){
			/* printf("\n == pair %i-%i ==\n", i+1,j+1); */
			res[k] = snpbin_dotprod(&dat.x[i], &dat.x[j], mean, sd);
			++k;
		}
	}

	/* add the diagonal to the end of the array */
	for(i=0; i< *nind; i++){
		/* printf("\n == pair %i-%i == \n", i+1,i+1); */
		res[k] = snpbin_dotprod(&dat.x[i], &dat.x[i], mean, sd);
		++k;
	}
}





/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/

/* Simple test function */
/* Test: increases for a raw (unsigned char) vector */
void testRaw(unsigned char *a, int *n){
	int i;
	for(i=0; i< *n; i++){
		a[i] = (unsigned char)(i);
	}
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


## === DOT PRODUCTS === ##

#### NO LONGER NEEDED ####

library(adegenet)
library(ade4)
dat <- rbind("a"=c(1,0,0), "b"=c(1,2,1), "c"=c(1,0,1))
x <- new("genlight",dat)


## NOT CENTRED, NOT SCALED
glDotProd(x)

res2 <- as.matrix(x) %*% t(as.matrix(x))
res2

## DATA > 8 SNPs
dat <- rbind(rep(c(1,0,1), c(8,10,5)))
x <- new("genlight",dat)
glDotProd(x)


## RANDOM DATA
dat <- matrix(sample(0:1, 5*1000, replace=TRUE), nrow=5)
x <- new("genlight",dat)
res1 <- glDotProd(x)

res2 <- as.matrix(x) %*% t(as.matrix(x))

all(res1==res2)


## CENTRED, NOT SCALED
res1 <- glDotProd(x, cent=TRUE)

temp <- as.matrix(x) / ploidy(x)
temp <- scalewt(temp, cent=TRUE, scale=FALSE)
res2 <- temp %*% t(temp)
res2

all(abs(res1-res2)<1e-10)


## CENTRED, SCALED
res1 <- glDotProd(x, cent=TRUE, scale=TRUE)

temp <- as.matrix(x) / ploidy(x)
temp <- scalewt(temp, cent=TRUE, scale=TRUE)
res2 <- temp %*% t(temp)
res2

all(abs(res1-res2)<1e-10)




## TEST WITH NAs
dat <- lapply(1:50, function(i) sample(c(0,1,NA), 1e4, prob=c(.5, .49, .01), replace=TRUE))
names(dat) <- paste("indiv", 1:length(dat))
print(object.size(dat), unit="aut") # size of the original data

x <- new("genlight", dat) # conversion
x

dat <- matrix(sample(0:1, 5*1000, replace=TRUE), nrow=5)
x <- new("genlight",dat)

*/

