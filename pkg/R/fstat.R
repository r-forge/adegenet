#################
# fstat function
#################
#
# Wrapper for fst estimator from hierfstat package
#
fstat <- function(x, pop=NULL, fstonly=FALSE){
    ## misc checks
    if(!is.genind(x)) stop("x is not a valid genind object")
    if(!require(hierfstat)) stop("hierfstat package is required. Please install it.")
    if(x@ploidy != as.integer(2)) stop("not implemented for non-diploid genotypes")
    checkType(x)

    if(is.null(pop)) pop <- x@pop
    if(is.null(pop)) stop("no pop factor provided")
    if(length(pop)!=nrow(x@tab)) stop("pop has a wrong length.")

    ## computations
    dat <- genind2hierfstat(x)[,-1]
    res <- varcomp.glob(levels=data.frame(pop), loci=dat)$F

    if(fstonly) {res <- res[1,1]}
    return(res)
}






###############
# fst function
###############
#
# classical fst sensu Nei
#
pairwise.fst <- function(x, pop=NULL){
    ## MISC CHECKS ##
    if(!is.genind(x)) stop("x is not a valid genind object")

    if(is.null(pop)) pop <- pop(x)
    if(is.null(pop)) stop("no pop factor provided")
    if(length(pop)!=nrow(x@tab)) stop("pop has a wrong length.")


    ## COMPUTATIONS ##



    return(res)
} # end of pairwise.fst
