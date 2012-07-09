#############
## GENERIC ##
#############
gengraph <-  function (x, ...) UseMethod("gengraph")





#############
## DEFAULT ##
#############
gengraph.default <- function(x, ...){
    stop(paste("No method for objects of class",class(x)))
} # end gengraph.default





#############
## DEFAULT ##
#############
gengraph.genind <- function(x, ...){
    ## CHECKS ##
    if(!require("igraph")) stop("igraph is required")

    ## COMPUTE DISTANCES ##
    temp <- 1-propShared(x)

}
