
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



## TESTING ##
## all.equal(glVar(x), apply(as.matrix(x), 2, function(e) mean((e-mean(e, na.rm=TRUE))^2, na.rm=TRUE)))
