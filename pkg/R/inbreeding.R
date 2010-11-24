#############
## inbreeding
#############
inbreeding <- function(x, pop=NULL, truenames=TRUE, res.type=c("mean","byloc"), plot=TRUE, ...){
    ## CHECKS ##
    if(!is.genind(x)) stop("x is not a valid genind object")
    checkType(x)
    res.type <- match.arg(res.type)

    if(x$ploidy != 2) stop("this inbreeding coefficient is designed for diploid genotypes only")

    if(!is.null(pop)) pop(x) <- pop
    if(is.null(x@pop) && is.null(pop)) {
        pop(x) <- factor(rep(1, nrow(x@tab)))
    }


    ## COMPUTATIONS ##

    ## get allele frequencies and \sum p_i^2 by pop and loc ##
    tabfreq2 <- (makefreq(x = genind2genpop(x, quiet = TRUE), quiet=TRUE, truenames=truenames)$tab) ^2
    sumpi2 <- t(apply(tabfreq2, 1, tapply, x$loc.fac, sum))

    ## function to check a 1-locus genotype for homozigosity
    ## returns 1 if homoz, 0 otherwise
    ## !!! NOTE : reverse the values returned by f1 to obtain a strange thing !!!
    f1 <- function(gen){
        if(any(is.na(gen))) return(NA)
        if(any(round(gen, 10)==1)) return(1)
        return(0)
    }

    ## get the table of binary hetero/homo data
    if(truenames) {
        X <- truenames(x)$tab
    } else
    X <- x$tab

    homotab <- t(apply(X, 1, tapply, x@loc.fac, f1))


    ## get pi2 for the appropriate pop
    if(truenames){
    popx <- pop(x)
    } else {
        popx <- x$pop
    }

    popx <- as.character(popx)
    tabpi2 <- sumpi2[popx, , drop=FALSE]


    ## COMPUTE FINAL RESULT ##
    num <- homotab - tabpi2
    denom <- tabpi2 * (1 - tabpi2)
    res <- num / denom
    if(res.type=="byloc") return(res)

    res <- apply(res, 1, mean, na.rm=TRUE)
    if(plot){
        par(bg="grey")
        nPop <- length(unique(popx))
        myCol <- rainbow(nPop)[as.integer(pop(x))]
        plot(res, col=myCol, type="h", ylab="Inbreeding", xlab="Individuals", ...)
    }

    return(res)

} # end inbreeding
