seqTrack2 <- function(x, x.names, x.dates, ...(){
    ## CHECKS ##
    if(!require(ape)) stop("package ape is required")

    ## FIND NUMBER OF MUTATIONS / SEQUENCE LENGTH ##
    L <- ifelse(is.matrix(x), ncol(x), length(x[[1]])) # haplo length
    D <- as.matrix(dist.dna(x, model="raw")) * L # nb of mutations
    N <- nrow(D)

    ## AUXILLIARY FUNCTIONS ##
    ## function generating infection durations
    aug.infdur <- function(){

    }

    ## function generating infection dates

    ## function moving nu1, nu2

    ## function moving alpha_i

    ## function moving pi


    ## ESTIMATION ##
    ## initialisation

    ## run mcmc

    ##

}







                      ## NUCL <- c('a','t','g','c')
    ## dist.nbmut <- function(a,b){
    ##     a[!a %in% NUCL] <- NA
    ##     b[!b %in% NUCL] <- NA
    ##     isNA <- is.na(a) | is.na(b)
    ##     a <- a[!isNA]
    ##     b <- b[!isNA]
    ##     return(sum(a!=b))
    ## }

    ## nbNuclCommon <- function(a,b){
    ##     NUCL <- c('a','t','g','c')
    ##     a[!a %in% NUCL] <- NA
    ##     b[!b %in% NUCL] <- NA
    ##     isNA <- is.na(a) | is.na(b)

    ##     return(sum(!isNA))
    ## }
