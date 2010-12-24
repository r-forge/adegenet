

###############
##
##   CLASSES
##
###############

###############
## SNPbin class
###############
setClass("SNPbin", representation(snp = "list",
                                  n.loc = "integer",
                                  NA.posi = "integer",
                                  label = "charOrNULL",
                                  ploidy = "integer"),
         prototype(snp = list(0), n.loc = integer(0), label = NULL, ploidy = 1L))




###############
## genlight class
###############
setClass("genlight", representation(tab = "list",
                                    n.loc = "integer",
                                    ind.names = "charOrNULL",
                                    loc.names = "charOrNULL",
                                    ploidy = "integer"),
         prototype(tab = list(0), n.loc = integer(0), ind.names = NULL, loc.names = NULL, ploidy = 1L))








#####################
##
##   CONSTRUCTORS
##
#####################

####################
## SNPbin constructor
####################
setMethod("initialize", "SNPbin", function(.Object, ...) {
    x <- .Object
    input <- list(...)
    if(length(input)==1) names(input) <- "snp"


    ## handle snp data ##
    if(!is.null(input$snp)){
        ## a vector of raw is provided
        if(is.raw(input$snp)){
            x@snp <-list(input$snp)
        }

        ## a list of raw vectors is provided
        if(is.list(input$snp)){
            if(all(sapply(input$snp, class)=="raw")){
                x@snp <- input$snp
            }
        }

        ## a numeric/integer vector is provided
        ## conversion from a vector of 0/1 (integers)
        if(is.numeric(input$snp) | is.integer(input$snp)){
            input$snp <- as.integer(input$snp)
            ## determine ploidy
            if(is.null(input$ploidy)){
                input$ploidy <- max(input$snp, na.rm=TRUE)
            }
            input$ploidy <- as.integer(input$ploidy)
            if(input$ploidy<1) stop("Ploidy is less than 1")

            ## check values in the vector
            if(any(input$snp<0 | input$snp>input$ploidy, na.rm=TRUE)){
                stop("Values of SNPs < 0 or > ploidy")
            }


            ## handle ploidy (may have to split info into binary vectors)
            x@snp <- list()
            i <- max(input$snp, na.rm=TRUE) # determine max nb of alleles
            if(i > 1){ # haplotype can be 0/1/2/...
                j <- 0 # index for the length of the list @snp
                while(i > 0){
                    j <- j+1 # update length of the result
                    temp <- as.integer(input$snp==i)
                    x@snp[[j]] <- .bin2raw(temp)$snp # make a vector of 1
                    input$snp <- input$snp-temp # deflate data (subtract the recoded alleles)
                    i <- max(input$snp, na.rm=TRUE) # update the max nb of alleles
                }
            } else { # haplotype is only 0/1/NA
                x@snp[[1]] <- .bin2raw(input$snp)$snp
            }
            x@n.loc <- length(input$snp)
            x@NA.posi <- which(is.na(input$snp))
            x@ploidy <- input$ploidy
            return(x)
        }
    }


    ## handle n.loc ##
    if(!is.null(input$n.loc)){
        x@n.loc <- as.integer(input$n.loc)
    } else {
        x@n.loc <- as.integer(length(x@snp)*8)
    }


    ## handle NA.posi ##
    if(!is.null(input$NA.posi)){
        x@NA.posi <- as.integer(input$NA.posi)
    }


    ## handle ploidy ##
    if(!is.null(input$ploidy)){
        x@ploidy <- as.integer(input$ploidy)
    }

    return(x)
}) # end SNPbin constructor














################################
##
##   METHODS AND ACCESSORS
##
################################

###############
## show SNPbin
###############
setMethod ("show", "SNPbin", function(object){
    cat(" === S4 class SNPbin ===")
    if(!is.null(object@label)) {
        cat("\n", object@label)
    }
    cat("\n", nLoc(object), "SNPs coded as bits")
    cat("\n Ploidy:", object@ploidy)
    temp <- round(length(object@NA.posi)/nLoc(object) *100,2)
    cat("\n ", length(object@NA.posi), " (", temp," %) missing data\n", sep="")
}) # end show method




############
## accessors
############
setMethod("nLoc","SNPbin", function(x,...){
    return(x@n.loc)
})


setMethod("nLoc","genlight", function(x,...){
    return(x@n.loc)
})


setMethod("$","SNPbin",function(x,name) {
    return(slot(x,name))
})


setMethod("$","genlight",function(x,name) {
    return(slot(x,name))
})


setMethod("names", signature(x = "SNPbin"), function(x){
    return(slotNames(x))
})


setMethod("names", signature(x = "genlight"), function(x){
    return(slotNames(x))
})





###############
## '[' operators
###############
## SNPbin
setMethod("[", signature(x="SNPbin", i="ANY"), function(x, i) {
    if (missing(i)) i <- TRUE
    temp <- .SNPbin2int(x) # data as integers with NAs
    x <- new("SNPbin", snp=temp[i], label=x@label, ploidy=x@ploidy)
    return(x)
}) # end [] for SNPbin







###################
##
##   CONVERSIONS
##
###################

############
## .bin2raw
###########
## each byte takes a value on [0,255]

## function to code multiple SNPs on a byte
## 8 combinations of SNPs can be coded onto a single byte (0->255)
.bin2raw <- function(vecSnp){
    ## was required in pure-R version
    ## SNPCOMB <- as.matrix(expand.grid(rep(list(c(0,1)), 8)))
    ## colnames(SNPCOMB) <- NULL

    ## handle missing data
    NAposi <- which(is.na(vecSnp))
    if(length(NAposi)>0){
        vecSnp[is.na(vecSnp)] <- 0L
    }


    nbBytes <- 1+ length(vecSnp) %/% 8
    ori.length <- length(vecSnp)
    new.length <- 8*nbBytes
    vecSnp <- c(vecSnp, rep(0, new.length-ori.length)) # fill the end with 0 of necessary


    ## map info to bytes (0:255)
    vecSnp <- as.integer(vecSnp)
    vecRaw <- integer(nbBytes)

    vecRaw <- .C("binIntToBytes", vecSnp, length(vecSnp), vecRaw, nbBytes, PACKAGE="adegenet")[[3]]
    ## vecraw <- sapply(seq(1, by=8, length=nbBytes), function(i) which(apply(SNPCOMB,1, function(e) all(temp[i:(i+7)]==e))) ) # old R version

    ## code information as raw and add missing data
    vecRaw <- as.raw(vecRaw)
    res <- list(snp=vecRaw, n.loc=as.integer(ori.length), NA.posi=as.integer(NAposi))
    return(res)
} # end .bin2raw






###########
## .raw2bin
###########
## convert vector of raw to 0/1 integers
.raw2bin <- function(x){
    SNPCOMB <- as.matrix(expand.grid(rep(list(c(0,1)), 8)))
    colnames(SNPCOMB) <- NULL
    res <- unlist(lapply(as.integer(x), function(i) SNPCOMB[i+1,]))
} # end .raw2bin





#############
## .SNPbin2int
#############
## convert SNPbin to integers (0/1)
.SNPbin2int <- function(x){
    res <- lapply(x@snp, .raw2bin)
    res <- lapply(res, function(e) e[1:x@n.loc])
    res <- as.integer(Reduce("+", res))
    if(length(x@NA.posi)>0){
        res[x@NA.posi] <- NA
    }
    return(res)
} # end .SNPbin2int






#############
## as methods
#############
setAs("SNPbin", "integer", def=function(from){
    res <- .SNPbin2int(from)
    return(res)
})


as.integer.SNPbin <- function(x, ...){
    return(as(x, "integer"))
}












################################
## testing :
##
## HAPLOID DATA
##
## library(adegenet)

## dat <- sample(c(0L,1L,NA), 1e6, prob=c(.5, .495, .005), replace=TRUE)
## x <- new("SNPbin", dat)

## identical(as(x, "integer"),dat) # SHOULD NORMALLY BE TRUE
## all(as(x, "integer") == dat, na.rm=TRUE) # MUST BE TRUE

## object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION



 ## DIPLOID DATA
## dat <- sample(c(0:2,NA), 1e6, prob=c(.4, .4, .195 ,.005), replace=TRUE)
## x <- new("SNPbin", dat)

## identical(as(x, "integer"),dat) # MUST BE TRUE

## object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION



 ## POLYLOID DATA
## dat <- sample(c(0:5,NA), 1e6, prob=c(rep(.995/6,6), 0.005), replace=TRUE)
##  x <- new("SNPbin", dat)

## identical(as(x, "integer"),dat) # MUST BE TRUE

## object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION
