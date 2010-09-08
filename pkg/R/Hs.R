############################
# Hs (expected heterozygosity)
############################
Hs <- function(x, truenames=TRUE) {

    ## CHECKS
    if(!is.genpop(x)) stop("x is not a valid genpop object")
    if(x@type=="PA") stop("not implemented for presence/absence markers")


    ## MAIN COMPUTATIONS

    ## have to handle loci with no polymorphism
    x.byloc <- seploc(x, truenames=truenames)

    toRemove <- which(x@loc.nall==1)
      if(length(toRemove)>0){
          x.byloc <- x.byloc[-toRemove]

    }

    lX <- lapply(x.byloc, function(e) makefreq(e, quiet=TRUE, truenames=truenames)$tab)
    lres <- lapply(lX, function(X) 1- apply(X^2,1,sum))

    res <- apply(as.matrix(data.frame(lres)),1,mean)

    return(res)
} # end Hs
