
##########
## glSum
##########
## compute col sums
## removing NAs
##
glSum <- function(x){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## DEFAULT, VECTOR-WISE PROCEDURE ##
    res <- integer(nLoc(x))
    for(e in x@gen){
        temp <- as.integer(e)
        temp[is.na(temp)] <- 0L
        res <- res + temp
    }

    names(res) <- locNames(x)
    return(res)

} # glSum





##########
## glNA
##########
## counts NB of NAs per column
##
glNA <- function(x){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## DEFAULT, VECTOR-WISE PROCEDURE ##
    res <- integer(nLoc(x))
    temp <- NA.posi(x)
    for(e in temp){
        if(length(e)>0){
            res[e] <- res[e] + 1
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
glMean <- function(x){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## DEFAULT, VECTOR-WISE PROCEDURE ##
    N <- nInd(x) - glNA(x)
    res <- glSum(x)/N
    names(res) <- locNames(x)
    return(res)

} # glMean





########
## glVar
########
## computes SNPs variances
## takes NAs into account
##
glVar <- function(x){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## DEFAULT, VECTOR-WISE PROCEDURE ##
    N <- nInd(x) - glNA(x)
    xbar <- glMean(x)

    res <- numeric(nLoc(x))
    for(e in x@gen){
        temp <- (as.integer(e) - xbar)^2
        temp[is.na(temp)] <- 0L
        res <- res + temp
    }

    res <- res/N
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
        vecMeans <- glMean(x)
        if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
    } else {
        vecMeans <- 0
    }
    if(scale){
        if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
        vecSd <- sqrt(glVar(x))
    } else {
        vecSd <- 1
    }


    ## COMPUTE DOT PRODUCTS BETWEEN GENOTYPES ##
    dotProd <- function(a,b){ # a and b are two SNPbin objects
        res <- sum( ((as.integer(a)-vecMeans)/vecSd) * ((as.integer(b)-vecMeans)/vecSd), na.rm=TRUE )
        return(res)
    }

    allComb <- combn(1:nInd(x), 2)
    allProd <- unlist(lapply(1:ncol(allComb), function(i) dotProd(x@gen[[allComb[1,i]]]), x@gen[[allComb[2,i]]]))
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

    ## scan  nb of axes retained
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
        ## use: c1 = X^TDV
        ##### TO FINISH !!!
        res$loadings <- matLoad
    }
    ## compute loadings

    names(res) <- locNames(x)
    return(res)

} # glPca




## TESTING ##
## all.equal(glVar(x), apply(as.matrix(x), 2, function(e) mean((e-mean(e, na.rm=TRUE))^2, na.rm=TRUE)))
