

###############
##
##   CLASSES
##
###############

###############
## SNPbin class
###############
setClass("SNPbin", representation(snp = "raw",
                                  n.loc = "integer",
                                  NA.posi = "integer",
                                  label = "charOrNULL"),
         prototype(snp = raw(0), n.loc = integer(0), label = NULL))




###############
## genlight class
###############
setClass("genlight", representation(tab = "list",
                                    n.loc = "integer",
                                    ind.names = "charOrNULL",
                                    loc.names = "charOrNULL"),
         prototype(tab = list(0), n.loc = integer(0), ind.names = NULL, loc.names = NULL))








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
        if(is.raw(input$snp)){
            x@snp <-input$snp
        }

        ## conversion from a vector of 0/1 (integers)
        if(is.numeric(input$snp) | is.integer(input$snp)){
            temp <- na.omit(unique(input$snp))
            if(!all(temp %in% c(0,1))){
                stop("Values of SNPs are not all 0, 1, or NA")
            }

            temp <- .bin2raw(input$snp)
            x@snp <- temp$snp
            x@n.loc <- temp$n.loc
            x@NA.posi <- temp$NA.posi
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
    if(!is.null(x@label)) {
        cat("\n", x@label)
    }
    cat("\n", nLoc(object), "SNPs coded as bits")
    temp <- round(length(object@NA.posi)/nLoc(object) *100,2)
    cat("\n", length(object@NA.posi), " (", temp," %) missing data\n", sep="")
}) # end show method




############
## accessors
############
setMethod("nLoc","SNPbin", function(x,...){
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
    x <- new("SNPbin", snp=temp[i], label=x@label)
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
## 8 combinations of SNPs can be coded on a single byte
## we use bytes values from [1,128]
## 200 is a missing value
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






#############
## .SNPbin2int
#############
## convert SNPbin to integers (0/1)
.SNPbin2int <- function(x){
    SNPCOMB <- as.matrix(expand.grid(rep(list(c(0,1)), 8)))
    colnames(SNPCOMB) <- NULL
    res <- unlist(lapply(as.integer(x@snp), function(i) SNPCOMB[i+1,]))
    res <- res[1:x@n.loc]
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
## dat <- sample(c(0,1,NA), 1e6, prob=c(.5, .495, .005), replace=TRUE)
## x <- new("SNPbin", dat)

## identical(as(x, "integer"),dat) # MUST BE TRUE

## object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION
