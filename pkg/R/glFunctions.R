
##########
## glSum
##########
## compute col sums
## removing NAs
##
glSum <- function(x, alleleAsUnit=TRUE){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## DEFAULT, VECTOR-WISE PROCEDURE ##
    ## use ploidy (sum absolute frequencies)
    if(alleleAsUnit){
    res <- integer(nLoc(x))
        for(e in x@gen){
            temp <- as.integer(e)
            temp[is.na(temp)] <- 0L
            res <- res + temp
        }
    } else {
        ## sum relative frequencies
        res <- numeric(nLoc(x))
        myPloidy <- ploidy(x)
        for(i in 1:nInd(x)){
            temp <- as.integer(x@gen[[i]]) / myPloidy[i]
            temp[is.na(temp)] <- 0
            res <- res + temp
        }
    }
    names(res) <- locNames(x)
    return(res)

} # glSum





##########
## glNA
##########
## counts NB of NAs per column
##
## if alleleAsUnit, then effective is the number of alleles sampled (sum(ploidy(x)))
## otherwise, effective is simply the number of individuals (
glNA <- function(x, alleleAsUnit=TRUE){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## DEFAULT, VECTOR-WISE PROCEDURE ##
    res <- integer(nLoc(x))
    temp <- NA.posi(x)

    ## NAs in allele sampling
    if(alleleAsUnit){
        for(i in 1:length(temp)){
            if(length(temp[[i]])>0){
                res[temp[[i]]] <- res[temp[[i]]] + ploidy(x)[i]
            }
        }
    } else { ## NAs amongst individuals
        for(e in temp){
            if(length(e)>0){
                res[e] <- res[e] + 1
            }
        }
    }

    names(res) <- locNames(x)
    return(res)

} # glNA





##########
## glMean
##########
## computes SNPs means
## takes NAs into account
##
glMean <- function(x, alleleAsUnit=TRUE){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## DEFAULT, VECTOR-WISE PROCEDURE ##
    if(alleleAsUnit){ # use alleles
        N <- sum(ploidy(x)) - glNA(x, alleleAsUnit=TRUE)
        res <- glSum(x, alleleAsUnit=TRUE)/N
    } else { # use relative frequencies of individuals
        N <- nInd(x) - glNA(x, alleleAsUnit=FALSE)
        res <- glSum(x, alleleAsUnit=FALSE)/N
    }

    names(res) <- locNames(x)
    return(res)

} # glMean





########
## glVar
########
## computes SNPs variances
## takes NAs into account
##
glVar <- function(x, alleleAsUnit=TRUE){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## DEFAULT, VECTOR-WISE PROCEDURE ##
    res <- numeric(nLoc(x))
    myPloidy <- ploidy(x)

    if(alleleAsUnit){ # use alleles
        N <- sum(ploidy(x)) - glNA(x, alleleAsUnit=TRUE)
        xbar <- glMean(x, alleleAsUnit=TRUE)
        for(i in 1:nInd(x)){
            temp <- (as.integer(x@gen[[i]])/myPloidy[i] - xbar)^2
            temp[is.na(temp)] <- 0
            res <- res + temp*myPloidy[i]
        }
        res <- res/N
    } else { # use relative frequencies of individuals
        N <- nInd(x) - glNA(x, alleleAsUnit=FALSE)
        xbar <- glMean(x, alleleAsUnit=FALSE)

        for(i in 1:nInd(x)){
            temp <- (as.integer(x@gen[[i]])/myPloidy[i] - xbar)^2
            temp[is.na(temp)] <- 0L
            res <- res + temp
        }
        res <- res/N
    }

    names(res) <- locNames(x)
    return(res)

} # glVar








########
## glPca
########
##
## PCA for genlight objects
##
glPca <- function(x, center=TRUE, scale=FALSE, scannf=TRUE, nf=2, loadings=TRUE){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## COMPUTE MEANS AND VARIANCES ##
    if(center) {
       vecMeans <- glMean(x, alleleAsUnit=FALSE)
        if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
    }

    if(scale){
        vecSd <- sqrt(glVar(x, alleleAsUnit=FALSE))
        if(any(is.na(vecSd))) stop("NAs detected in the vector of variances")
    }

    ## COMPUTE DOT PRODUCTS BETWEEN GENOTYPES ##
    ## to be fast, a particular function is defined for each case of centring/scaling

    myPloidy <- ploidy(x)

    ## NO CENTRING / NO SCALING
      if(!centre & !scale){
        dotProd <- function(a,b, ploid.a, ploid.b){ # a and b are two SNPbin objects
            a <- as.integer(a) / ploid.a
            a[is.na(a)] <- 0
            b <- as.integer(b) / ploidy.b
            b[is.na(b)] <- 0
            res <- sum( a*b, na.rm=TRUE )
            return(res)
        }
    }

    ## CENTRING / NO SCALING
    if(centre & !scale){
        dotProd <- function(a,b, ploid.a, ploid.b){ # a and b are two SNPbin objects
            a <- as.integer(a) / ploid.a
            a[is.na(a)] <- vecMeans[is.na(a)]
            b <- as.integer(b) / ploid.b
            b[is.na(b)] <- vecMeans[is.na(b)]
            res <- sum( (a-vecMeans) * (b-vecMeans), na.rm=TRUE )
            return(res)
        }
    }


    ## NO CENTRING / SCALING (odd option...)
    if(!centre & scale){
        dotProd <- function(a,b, ploid.a, ploid.b){ # a and b are two SNPbin objects
            a <- as.integer(a) / ploid.a
            a[is.na(a)] <- 0
            b <- as.integer(b) / ploid.b
            b[is.na(b)] <- 0
            res <- sum( a/vecSd * b/vecSd, na.rm=TRUE )
            return(res)
        }
    }


    ## CENTRING/ SCALING
    if(centre & scale){
        dotProd <- function(a,b, ploid.a, ploid.b){ # a and b are two SNPbin objects
            a <- as.integer(a) / ploid.a
            a[is.na(a)] <- vecMeans[is.na(a)]
            b <- as.integer(b) / ploid.b
            a[is.na(a)] <- vecMeans[is.na(a)]
            res <- sum( (a-vecMeans)/vecSd * (b-vecMeans)/vecSd, na.rm=TRUE )
            return(res)
        }
    }


    ## COMPUTE ALL POSSIBLE DOT PRODUCTS
    allComb <- combn(1:nInd(x), 2)
    allProd <- unlist(lapply(1:ncol(allComb), function(i) dotProd(x@gen[[allComb[1,i]]], x@gen[[allComb[2,i]]], myPloidy[allComb[1,i]], myPloidy[allComb[2,i]]) ))
    allProd <- allProd / nInd(x) # assume uniform weights

    ## shape result as a matrix
    attr(allProd,"Size") <- nInd(x)
    attr(allProd,"Diag") <- FALSE
    attr(allProd,"Upper") <- FALSE
    class(allProd) <- "dist"
    allProd <- as.matrix(allProd)


    ## PERFORM EIGENANALYSIS ##
    ## eigenanalysis
    eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
    rank <- sum(eigRes > 1e-10)
    eigRes$values <- eigRes$value[1:rank]
    eigRes$vectors <- eigRes$vectors[1:rank]

    ## scan nb of axes retained
    if(scannf){
        barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
        cat("Select the number of axes: ")
        nf <- as.integer(readLines(n = 1))
    }

    ## rescale PCs
    res <- list()
    ## use: li = XQU = V\Lambda^(1/2)
    res$scores <- sweep(eigRes$vectors[, 1:nf],2, sqrt(eigRes$values[1:nf]), FUN="*")

    ## get loadings
    if(loadings){
        res$loadings <- matrix(0
        ## use: c1 = X^TDV
    ##     ##### TO FINISH !!!
    ##     for(i in 1:nInd(x)){
    ##         temp <- as.integer(x@gen[[i]]) / myPloidy[i]
    ##         if(centre) {
    ##             temp[is.na(temp)] <- vecMeans[is.na(temp)]
    ##             temp <- temp - vecMeans
    ##         } else {
    ##             temp[is.na(temp)] <- 0
    ##         }
    ##         if(scale){
    ##             temp <- temp/vecSd
    ##         }

    ##         res$loadings <- t(temp) %*% res$scores / N
    ##     }
    ## }
    ## ## compute loadings

    names(res) <- locNames(x)
    return(res)

} # glPca




## TESTING ##
## x <- new("genlight", list(c(0,0,1,1,0), c(1,1,1,0,0,1), c(2,1,1,1,1,NA)))
## as.matrix(x)
## glNA(x)
## glSum(x)
## glMean(x)
##
## same ploidy everywhere
## x <- new("genlight", list(c(0,0,1,1,0), c(1,1,1,0,0,1), c(0,0,0,1,1,1)))
## f1 <- function(e) {return(mean((e-mean(e, na.rm=TRUE))^2, na.rm=TRUE))}
## all.equal(glVar(x), apply(as.matrix(x), 2, f1 )) # MUST BE TRUE
## all.equal(glVar(x,FALSE), apply(as.matrix(x), 2, f1 )) # MUST BE TRUE

## ## differences in ploidy
## x <- new("genlight", list(c(0,0,1,1,0), c(1,1,1,0,0,1), c(2,1,1,1,1,NA)))
## temp <- sweep(as.matrix(x), 1, c(1,1,2), "/")
## f2 <- function(e,w) {
##     mu <- weighted.mean(e, w, na.rm=TRUE)
##     res <- weighted.mean((e-mu)^2, w, na.rm=TRUE)
##     return(res)
## }

## all.equal(glMean(x), apply(temp,2,weighted.mean, w=c(1,1,2), na.rm=TRUE)) # MUST BE TRUE
## all.equal(glVar(x), apply(temp, 2, f2,w=c(1,1,2) )) # MUST BE TRUE

## all.equal(glMean(x,FALSE), apply(temp,2,mean,na.rm=TRUE)) # MUST BE TRUE
## all.equal(glVar(x,FALSE), apply(temp,2,f1)) # MUST BE TRUE
