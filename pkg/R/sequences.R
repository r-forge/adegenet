######################################
##
## The code below implements import
## from alignement data.
##
######################################



################
# DNAbin2genind
################
DNAbin2genind <- function(x, pop=NULL, na.char=c("n","-","?"), polyThres=1/100){

    ## misc checks
    if(!inherits(x,"DNAbin")) stop("x is not a DNAbin object")
    if(!require(ape)) stop("The package ape is required.")

    ## DNA bin to matrix of characters
    x <- as.character(x) # should output a matrix

    if(is.list(x)) { # if this is a list
        temp <- unique(sapply(x,length)) # check lengths of sequences
        if(length(temp)>1) stop("Sequences have different length - please use alignements only.")
        else{ # if sequences have same length, build the matrix
            temp <- names(x)
            x <- t(as.data.frame(x))
            rownames(x) <- temp
        }
    }

    if(is.null(colnames(x))) {
        colnames(x) <- 1:ncol(x)
    }

    ## replace NAs
    x[x %in% na.char] <- NA

    ## keep only columns with polymorphism (i.e., SNPs)
    isPoly <- function(vec){
        N <- sum(!is.na(vec)) # N: number of sequences
        temp <- table(vec)/N
        if(sum(temp > polyThres) >= 2) return(TRUE)
        return(FALSE)
    }

    toKeep <- apply(x, 2, isPoly)
    x <- x[,toKeep]

    ## build output
    res <- df2genind(x, pop=pop, ploidy=1, ncode=1)
    res$call <- match.call()

    return(res)
} # end DNAbin2genind
