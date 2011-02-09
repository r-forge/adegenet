/*
  Coded by Thibaut Jombart (tjombart@imperial.ac.uk), December 2010.
  Distributed with the adephylo package for the R software.
  Licence: GPL >=2.

  Functions based on snpbin and genlightC classes, which mirror the R classes SNPbin and genlight on the C side.
*/


#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "snpbin.h"



/* Function to compute all dot products between individuals */
/* centring and scaling is always used */
/* but need to pass vectors of 0 and 1*/
void GLdotProd(unsigned char *gen, int *nbvecperind, int *byteveclength, int *nbnaperind, int *naposi, int *nind, int *nloc, int *ploidy, double *mean, double *sd, bool *freq, double *res){
	struct genlightC dat;
	int i, j, k=0;

	/* Check variance vector: do not divide by 0 */
	for(i=0;i< *nloc;i++){
		if(sd[i] < NEARZERO){
			sd[i] = 1;
		}
	}

	dat = genlightTogenlightC(gen, nbvecperind, byteveclength, nbnaperind, naposi, nind, nloc, ploidy);

	/* Lower triangle - without the diagonal */
	for(i=0; i< (*nind-1); i++){
		for(j=i+1; j< *nind; j++){
			/* printf("\n == pair %i-%i ==\n", i+1,j+1); */
			res[k] = snpbin_dotprod(&dat.x[i], &dat.x[j], mean, sd, freq);
			++k;
		}
	}

	/* add the diagonal to the end of the array */
	for(i=0; i< *nind; i++){
		/* printf("\n == pair %i-%i == \n", i+1,i+1); */
		res[k] = snpbin_dotprod(&dat.x[i], &dat.x[i], mean, sd, freq);
		++k;
	}
}






/* TESTING in R */

/*

## === DOT PRODUCTS === ##

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
library(adegenet)
library(ade4)

dat <- list(a=c(1,NA,0,0,2), b=c(1,2,3,4,0), c=c(NA,0,1,NA,2))

x <- new("genlight", dat) # conversion
x

res1 <- glDotProd(x)
t(data.frame(dat))
res1


*/

