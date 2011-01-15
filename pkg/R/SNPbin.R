

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
         prototype(snp = list(), n.loc = integer(0), label = NULL, ploidy = 1L))




###############
## genlight class
###############
setClass("genlight", representation(gen = "list",
                                    n.loc = "integer",
                                    ind.names = "charOrNULL",
                                    loc.names = "charOrNULL",
                                    loc.all = "charOrNULL",
                                    ploidy = "intOrNULL"),
         prototype(gen = list(), n.loc = integer(0), ind.names = NULL, loc.names = NULL, loc.all = NULL, ploidy=NULL))








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
    if(!is.null(input$snp) && length(input$snp)>0){
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
                if(input$ploidy==0) input$ploidy <- 1
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







####################
## genlight constructor
####################
setMethod("initialize", "genlight", function(.Object, ...) {
    x <- .Object
    input <- list(...)
    if(length(input)==1) names(input) <- "gen"


    ## HANDLE INPUT$GEN ##
    if(!is.null(input$gen)){
        ## input$gen is a list of SNPbin ##
        if(is.list(input$gen) && all(sapply(input$gen, class)=="SNPbin")){
            ## check nb of loci in each SNPbin
            if(length(unique(sapply(input$gen, nLoc)))>1) {
                warning("SNPbin objects have different numbers of loci")
                input$gen <- lapply(input$gen, as.integer)
            } else { # all seems fine
                x@gen <- input$gen
                if(is.null(input$ind.names)){
                    input$ind.names <- names(input$gen)
                }
            }
        }


        ## input$gen is a matrix or a data.frame
        if(is.matrix(input$gen) | is.data.frame(input$gen)){
            if(is.null(input$ind.names)){
                input$ind.names <- rownames(input$gen)
            }
            if(is.null(input$loc.names)){
                input$loc.names <- colnames(input$gen)
                if(is.data.frame(input$gen)){ # do not use names if these are the default names of a data.frame
                    if(identical(colnames(input$gen), paste("V", 1:ncol(input$gen), sep=""))){
                        input$loc.names <- NULL
                    }
                }
            }
            ##input$gen <- lapply(1:nrow(input$gen), function(i) as.integer(input$gen[i,]))
            x@gen <- lapply(1:nrow(input$gen), function(i) new("SNPbin", as.integer(input$gen[i,])) )
        }


        ## input$gen is a list of integers/numeric ##
        if(is.list(input$gen) && all(sapply(input$gen, class) %in% c("integer","numeric"))){
            ## check length consistency
            lengthvec <- sapply(input$gen, length)

            ## complete with NA is necessary
            if(length(unique(lengthvec))>1) {
                warning("Genotypes have variable length; completing shorter ones with NAs.")
                for(i in 1:length(input$gen)){
                    input$gen[[i]] <- c(input$gen[[i]], rep(NA, max(lengthvec)-length(input$gen[[i]])))
                }
            }

            ## name individuals if needed
            if(is.null(input$ind.names)){
                input$ind.names <- names(input$gen)
            }

            ## create SNPbin list
            x@gen <- lapply(input$gen, function(e) new("SNPbin",e))
        }
    }


    if(length(x@gen) > 0) { # if non-emtpy object
        ## HANDLE INPUT$IND.NAMES ##
        if(!is.null(input$ind.names)){
            input$ind.names <- as.character(input$ind.names)

            ## check length consistency
            if(length(input$ind.names) != nInd(x)){
                stop("Inconsistent length for ind.names.")
            } else {
                ## assign value to the output object
                x@ind.names <- input$ind.names
                ## ## name list and each SNPbin ## THIS DUPLICATES THE INFORMATION
                ## names(x@gen) <- input$ind.names
                ## for(i in 1:length(x@gen)){
                ##     x@gen[[i]]@label <- input$ind.names[i]
                ## }

            }


        }


        ## HANDLE INPUT$N.LOC ##
        if(!is.null(input$n.loc)){ # n.loc is provided
            input$n.loc <- as.integer(input$n.loc)

            ## check length consistency
            if(input$n.loc != nLoc(x@gen[[1]])) {
                warning("Inconsistent number of loci (n.loc) - ignoring this argument.")
            } else {
                x@n.loc <- input$n.loc
            }
        } else { # n.loc is not provided
            x@n.loc <- nLoc(x@gen[[1]])
        }


        ## HANDLE INPUT$PLOIDY ##
        ## note: if not provided, @ploidy is NULL (saves some space)
        if(!is.null(input$ploidy)){ # ploidy is provided
            input$ploidy <- as.integer(input$ploidy)
            input$ploidy <- rep(input$ploidy, length=length(x@gen))
            x@ploidy <- input$ploidy
        }


        ## HANDLE INPUT$LOC.NAMES ##
        if(!is.null(input$loc.names) && length(input$loc.names)>0){ # ploidy is provided
            input$loc.names <- as.character(input$loc.names)

            ## check length consistency
            if(length(input$loc.names) != x@n.loc){
                warning("Inconsistent length for loc.names - ignoring this argument.")
            } else {
                x@loc.names <- input$loc.names
            }
        }


        ## HANDLE INPUT$LOC.ALL ##
        if(!is.null(input$loc.all) && length(input$loc.all)>0){ # ploidy is provided
            input$loc.all <- as.character(input$loc.all)

            ## check length consistency
            if(length(input$loc.all) != x@n.loc){
                warning("Inconsistent length for loc.all - ignoring this argument.")
            } else {
                ## check string consistency (format is e.g. "a/t")
                if(any(grep("^[[:alpha:]]{1}/[[:alpha:]]{1}$", input$loc.all) != 1:length(x@gen))){
                    input$loc.all <- gsub("[[:space:]]","", input$loc.all)
                    warning("Miss-formed strings in loc.all (must be e.g. 'c/g') - ignoring this argument.")
                } else {
                    x@loc.all <- input$loc.all
                }
            }
        }

    } # end if non-empty @gen


    ## RETURN OBJECT ##
    names(x@gen) <- NULL # do not store ind.names twice
    return(x)
}) # end genlight constructor










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





###############
## show genlight
###############
setMethod ("show", "genlight", function(object){
    cat(" === S4 class genlight ===")
    cat("\n", nInd(object), "genotypes with", nLoc(object),  "binary SNPs")
    temp <- unique(ploidy(object))
    if(length(temp)==1){
        cat("\n Ploidy:", temp)
    } else {
        temp <- summary(ploidy(object))
        cat("\n Ploidy statistics (min/median/max):", temp[1], "/", temp[3], "/", temp[6])
    }
    temp <- sapply(object@gen, function(e) length(e@NA.posi))
    cat("\n ", sum(temp), " (", round(sum(temp)/(nInd(object)*nLoc(object)),2)," %) missing data\n", sep="")
}) # end show method





############
## accessors
############

## nLoc
setMethod("nLoc","SNPbin", function(x,...){
    return(x@n.loc)
})

setMethod("nLoc","genlight", function(x,...){
    return(x@n.loc)
})


## nInd
setMethod("nInd","genlight", function(x,...){
    return(length(x@gen))
})


## $
setMethod("$","SNPbin",function(x,name) {
    return(slot(x,name))
})

setMethod("$","genlight",function(x,name) {
    return(slot(x,name))
})

setMethod("$<-","SNPbin",function(x,name,value) {
  slot(x,name,check=TRUE) <- value
  return(x)
})

setMethod("$<-","genlight",function(x,name,value) {
  slot(x,name,check=TRUE) <- value
  return(x)
})


## names
setMethod("names", signature(x = "SNPbin"), function(x){
    return(slotNames(x))
})

setMethod("names", signature(x = "genlight"), function(x){
    return(slotNames(x))
})


## ploidy
setMethod("ploidy","SNPbin", function(x,...){
    return(x@ploidy)
})

setMethod("ploidy","genlight", function(x,...){
    if(!is.null(x@ploidy)){
        res <- x@ploidy
    } else {
        res <- sapply(x@gen, function(e) e@ploidy)
    }
    names(res) <- x@ind.names
    return(res)
})


setMethod("ploidy<-","SNPbin",function(x,value, ...) {
    value <- as.integer(value)
    if(any(value)<1) stop("Negative or null values provided")
    if(any(is.na(value))) stop("NA values provided")
    if(length(value)>1) warning("Several ploidy numbers provided; using only the first integer")
    slot(x,"ploidy",check=TRUE) <- value[1]
    return(x)
})

setMethod("ploidy<-","genlight",function(x,value, ...) {
    value <- as.integer(value)
    if(any(value)<1) stop("Negative or null values provided")
    if(any(is.na(value))) stop("NA values provided")
    if(length(value) != nInd(x)) stop("Length of the provided vector does not match nInd(x)")
    slot(x,"ploidy",check=TRUE) <- value
    return(x)
})


## locNames
setMethod("locNames","genlight", function(x,...){
    return(x@loc.names)
})


setMethod("locNames<-","genlight",function(x,value, ...) {
    value <- as.character(value)
    if(length(value) != nLoc(x)) stop("Vector length does no match number of loci")
    slot(x,"loc.names",check=TRUE) <- value
    return(x)
})


## indNames
setMethod("indNames","genlight", function(x,...){
    return(x@ind.names)
})


setMethod("indNames<-","genlight",function(x,value, ...) {
    value <- as.character(value)
    if(length(value) != nInd(x)) stop("Vector length does no match number of individuals")
    slot(x,"ind.names",check=TRUE) <- value
    return(x)
})


## allNames
setMethod("allNames","genlight", function(x,...){
    return(x@loc.all)
})


## NA.posi
setGeneric("NA.posi", function(x, ...) standardGeneric("NA.posi"))

setMethod("NA.posi","SNPbin", function(x,...){
    return(x@NA.posi)
})

setMethod("NA.posi","genlight", function(x,...){
    res <- lapply(x@gen, function(e) e@NA.posi)
    names(res) <- indNames(x)
    return(res)
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





## lightgen
setMethod("[", signature(x="genlight", i="ANY", j="ANY", drop="ANY"), function(x, i, j, ...) {
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE

    ## subset individuals
    x@gen <- x@gen[i]
    x@ind.names <- x@ind.names[i]
    if(!is.null(x@ploidy)) {
        ori.ploidy <- ploidy(x)[i]
    } else {
        ori.ploidy <- NULL
    }

    ## subset loci
    x <- as.matrix(x)[, j, drop=FALSE]
    x <- x[!apply(x, 1, function(e) all(is.na(e))), , drop=FALSE] # remove indiv that are all NAs
    x <- x[, !apply(x, 2, function(e) all(is.na(e))), drop=FALSE] # remove loci that are all NAs

    x <- new("genlight", gen=x, ploidy=ori.ploidy)

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
    ## handle missing data
    NAposi <- which(is.na(vecSnp))
    if(length(NAposi)>0){
        vecSnp[is.na(vecSnp)] <- 0L
    }


    nbBytes <- length(vecSnp) %/% 8
    if(length(vecSnp) %% 8 > 0) {nbBytes <- nbBytes +1}
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


setAs("genlight", "matrix", def=function(from){
    res <- unlist(lapply(from@gen, as.integer))
    res <- matrix(res, ncol=nLoc(from), nrow=nInd(from), byrow=TRUE)
    colnames(res) <- locNames(from)
    rownames(res) <- indNames(from)
    return(res)
})


## KLUDGE - needed for as.matrix.genlight to be dispatched correctly (R-2.12.1)
setGeneric("as.matrix")

as.matrix.genlight <- function(x, ...){
    return(as(x, "matrix"))
}


setAs("genlight", "data.frame", def=function(from){
    return(as.data.frame(as.matrix(from)))
})


as.data.frame.genlight <- function(x, ...){
    return(as(x, "data.frame"))
}


setAs("genlight", "list", def=function(from){
    res <- lapply(from@gen, as.integer)
    names(res) <- indNames(from)
    return(res)
})


as.list.genlight <- function(x, ...){
    return(as(x, "list"))
}













################################
## testing SNPbin
##
##
## library(adegenet)

## HAPLOID DATA - NO NA
## dat <- sample(c(0L,1L), 1e6, replace=TRUE)
## x <- new("SNPbin", dat)
## identical(as(x, "integer"),dat) # SHOULD NORMALLY BE TRUE
## all(as(x, "integer") == dat, na.rm=TRUE) # MUST BE TRUE
## object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION


## HAPLOID DATA - WITH NAs
## dat <- sample(c(0,1,NA), 1e6, prob=c(.5, .49, .01), replace=TRUE)
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




################################
## testing genlight
##
##


## SIMPLE TEST
## library(adegenet)
## dat <- list(toto=c(1,1,0,0), titi=c(NA,1,1,0), tata=c(NA,0,3, NA))
## x <- new("genlight", dat)
## x
## as.list(x)
## as.matrix(x)

## identical(x, new("genlight", as.list(x))) # round trip - list - MUST BE TRUE
## identical(x, new("genlight", as.matrix(x))) # round trip - matrix - MUST BE TRUE
## identical(x, new("genlight", as.data.frame(x))) # round trip - data.frame - MUST BE TRUE

## ## test subsetting
## identical(as.list(x[c(1,3)]), as.list(x)[c(1,3)]) # MUST BE TRUE
## identical(x, x[]) # MUST BE TRUE
## all.equal(t(as.matrix(as.data.frame(dat)))[,1:3], as.matrix(x[,1:3])) # MUST BE TRUE




## ## ## BIG SCALE TEST - HAPLOID DATA WITH NA
## library(adegenet)
## dat <- lapply(1:50, function(i) sample(c(0,1,NA), 1e6, prob=c(.5, .49, .01), replace=TRUE))
## names(dat) <- paste("indiv", 1:length(dat))
## print(object.size(dat), unit="aut")

## system.time(x <- new("genlight", dat)) # conversion + time taken
## print(object.size(x), unit="au")
## object.size(dat)/object.size(x) # conversion efficiency


## ## time taken by subsetting (quite long, +- 35sec)
## system.time(y <- x[1:10, 1:5e5])
