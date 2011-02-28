
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




## genlight
setMethod("[", signature(x="genlight", i="ANY", j="ANY", drop="ANY"), function(x, i, j, ..., treatOther=TRUE, quiet=TRUE, drop=FALSE) {
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE

    ori.n <- nInd(x)


    ## SUBSET INDIVIDUALS ##
    x@gen <- x@gen[i]
    x@ind.names <- x@ind.names[i]
    if(!is.null(x@ploidy)) {
        ori.ploidy <- ploidy(x)[i]
    } else {
        ori.ploidy <- NULL
    }

    ## HANDLE 'OTHER' SLOT ##
    nOther <- length(other(x))
    namesOther <- names(other(x))
    counter <- 0
    if(treatOther & !(is.logical(i) && all(i))){
        f1 <- function(obj,n=ori.n){
            counter <<- counter+1
            if(!is.null(dim(obj)) && nrow(obj)==ori.n) { # if the element is a matrix-like obj
                obj <- obj[i,,drop=FALSE]
            } else if(length(obj) == ori.n) { # if the element is not a matrix but has a length == n
                obj <- obj[i]
                if(is.factor(obj)) {obj <- factor(obj)}
            } else {if(!quiet) warning(paste("cannot treat the object",namesOther[counter]))}

            return(obj)
        } # end f1

        other(x) <- lapply(x@other, f1) # treat all elements

    } # end treatOther


    ## SUBSET LOCI ##
    if(length(j)==1 && is.logical(j) && j){ # no need to subset SNPs
        return(x)
    } else { # need to subset SNPs
        old.other <- other(x)
        x <- as.matrix(x)[, j, drop=FALSE] # maybe need to process one row at a time
        x <- new("genlight", gen=x, ploidy=ori.ploidy, other=old.other)
    }

    return(x)
}) # end [] for genlight







######################
##
## c, cbind, rbind...
##
######################

################
## cbind SNPbin
################
##setMethod("cbind", signature("SNPbin"), function(..., deparse.level = 1) {
cbind.SNPbin <- function(..., checkPloidy=TRUE){
    myList <- list(...)
    if(!all(sapply(myList, class)=="SNPbin")) stop("some objects are not SNPbin objects")
    ## remove empty objects
    myList <- myList[sapply(myList,nLoc)>0]
    if(length(myList)==0) {
        warning("All objects are empty")
        return(NULL)
    }


    if(checkPloidy && length(unique(sapply(myList, ploidy))) !=1 ) stop("objects have different ploidy levels")
    x <- new("SNPbin", unlist(lapply(myList, as.integer)))
    return(x)
} # end cbind.SNPbin
##})



c.SNPbin <- function(...){
    return(cbind(...))
}




##################
## cbind genlight
##################
##setMethod("cbind", signature(x="genlight"), function(..., deparse.level = 1) {
cbind.genlight <- function(...){
    myList <- list(...)
    if(!all(sapply(myList, class)=="genlight")) stop("some objects are not genlight objects")
    ## remove empty objects
    myList <- myList[sapply(myList,nLoc)>0 & sapply(myList,nInd)>0]
    if(length(myList)==0) {
        warning("All objects are empty")
        return(NULL)
    }

    ## different checks
    if(length(unique(sapply(myList, nInd))) > 1 ) stop("objects have different numbers of individuals")
    n.obj <- length(myList)
    n.ind <- nInd(myList[[1]])
    if(n.ind==0){
        warning("All objects are empty")
        return(NULL)
    }
    temp <- as.matrix(as.data.frame(lapply(myList, ploidy)))
    if(any(apply(temp,1,function(r) length(unique(r)))>1)) stop("non-consistent ploidy across datasets")


    ## merge one individual at a time ##
    res <- list()
    for(i in 1:n.ind){
        res[[i]] <- Reduce(function(a,b) {cbind(a,b,checkPloidy=FALSE)}, lapply(myList, function(e) e@gen[[i]]) )
    }

    res <- new("genlight",res)

    ## handle loc.names, alleles, etc. ##
    indNames(res) <- indNames(myList[[1]])
    locNames(res) <- unlist(lapply(myList, locNames))
    alleles(res) <- unlist(lapply(myList, alleles))
    pop(res) <- pop(myList[[1]])

    ## return object ##
    return(res)
} # end cbind.genlight
##})






##################
## rbind genlight
##################
##setMethod("cbind", signature(x="genlight"), function(..., deparse.level = 1) {
rbind.genlight <- function(...){
    myList <- list(...)
    if(!all(sapply(myList, class)=="genlight")) stop("some objects are not genlight objects")
    ## remove empty objects
    myList <- myList[sapply(myList,nLoc)>0 & sapply(myList,nInd)>0]
    if(length(myList)==0) {
        warning("All objects are empty")
        return(NULL)
    }

    if(length(unique(sapply(myList, nLoc))) !=1 ) stop("objects have different numbers of SNPs")


    ## build output
    res <- new("genlight", Reduce(c, lapply(myList, function(e) e@gen)))
    locNames(res) <- locNames(myList[[1]])
    alleles(res) <- alleles(myList[[1]])
    indNames(res) <- unlist(lapply(myList, indNames))
    pop(res) <- factor(unlist(lapply(myList, pop)))

    ## return object ##
    return(res)

} # end rbind.genlight






##########
## seppop
##########
setMethod("seppop", signature(x="genlight"), function(x, pop=NULL, treatOther=TRUE, quiet=TRUE){
    ## HANDLE POP ARGUMENT ##
    if(!is.null(pop)) {
        pop(x) <- pop
    }

    if(is.null(pop(x))) stop("pop not provided and pop(x) is NULL")

    ## PERFORM SUBSETTING ##
    kObj <- lapply(levels(pop(x)), function(lev) x[pop(x)==lev, , treatOther=treatOther, quiet=quiet])
    names(kObj) <- levels(pop(x))

    return(kObj)
})








##########
## glSim
##########
glSim <- function(n.ind, n.snp.nonstruc, n.snp.struc=0, grp.size=round(n.ind/2), ploidy=1, alpha=0,
                  block.size=NULL){

    ## BASIC CHECKS ##
    if( any(c(n.ind, n.snp.nonstruc+n.snp.struc) <1)) stop("null numbers of individuals and/or SNPs requested")
    ## alpha parameter
    if(alpha>0.5){
        alpha <- 0.5
        warning("alpha cannot exceed 0.5 - changing the value to 0.5 (total forced asymmetry)")
    }
    if(alpha<0){
        alpha <- 0
        warning("alpha cannot be lower than 0 - changing the value to 0 (no forced asymmetry)")
    }

    ## handle block size
    if(is.null(block.size)){
        block.size <- n.snp.nonstruc + n.snp.struc
    }

    ## handle group sizes
    grpA.size <- grp.size
    if(grpA.size >= n.ind) stop("grpA.size is >= n.ind")
    grpB.size <- n.ind - grpA.size


    ## AUXIL FUNCTIONS ##
    ## draw p snp for i indiv and convert into a genlight - no structure
    f1 <- function(n,p){
        temp <- sapply(1:p, function(i) rbinom(n, ploidy, runif(1)))
        if(n==1) {temp <- matrix(temp,nrow=1)}
        return(new("genlight", temp, ploidy=ploidy))
    }

    ## draw p snp for i indiv and convert into a genlight - differences between 2 groups
    if(n.snp.struc > 0){
        f2 <- function(n,p){
            probA <- runif(p, min=0, max=0.5-alpha)
            probB <- 1-probA
            tempA <- sapply(probA, function(i) rbinom(grpA.size, ploidy, i) )
            if(grpA.size==1) {tempA <- matrix(tempA,nrow=1)}
            tempB <- sapply(probB, function(i) rbinom(grpB.size, ploidy, i) )
            if(grpB.size==1) {tempB <- matrix(tempB,nrow=1)}
            return(new("genlight", rbind(tempA,tempB), ploidy=ploidy))
        }
    }

    ## NON-STRUCTURED DATA ##
    ## generate data
    if(n.snp.nonstruc <= block.size){ # no need to use blocks
        res.ns <- f1(n.ind, n.snp.nonstruc)
    } else { # proceed by blocks of SNPs
        res.ns <- f1(n.ind, block.size)
        snp.to.go <- n.snp.nonstruc - nLoc(res.ns) # number of SNPs left to simulate
        while(snp.to.go > 0){
            nb.snp.tmp <- min(snp.to.go, block.size)
            res.ns <- cbind(res.ns, f1(n.ind, nb.snp.tmp))
            snp.to.go <- n.snp.nonstruc - nLoc(res.ns) # number of SNPs left to simulate
        }
    } # end non-structured data


    ## STRUCTURED DATA ##
    if(n.snp.struc > 0){
        if(n.snp.struc <= block.size){ # no need to use blocks
            res.s <- f2(n.ind, n.snp.struc)
        } else { # proceed by blocks of SNPs
            res.s <- f2(n.ind, block.size)
            snp.to.go <- n.snp.struc - nLoc(res.s) # number of SNPs left to simulate
            while(snp.to.go > 0){
                nb.snp <- min(snp.to.go, block.size)
                res.s <- cbind(res.s, f2(n.ind, nb.snp))
                snp.to.go <- n.snp.struc - nLoc(res.s) # number of SNPs left to simulate
            }

        }
    } # end structured data

    ## RETURN RESULTS ##
    res <- res.ns
    if(n.snp.struc > 0){
        res <- cbind(res,res.s)
    }
    if(!is.null(grp.size)){
        pop(res) <- rep(c("A","B"), c(grpA.size, grpB.size))
    }
    indNames(res) <- paste("ind", 1:n.ind)

    return(res)

} # end glSim






###################
### TESTING
###################


## c, cbind, rbind ##
## a <- new("genlight", list(c(1,0,1), c(0,0,1,0)) )
## b <- new("genlight", list(c(1,0,1,1,1,1), c(1,0)) )
## locNames(a) <- letters[1:4]
## locNames(b) <- 1:6
## c <- cbind(a,b)
## identical(as.matrix(c),cbind(as.matrix(a), as.matrix(b))) # MUST BE TRUE
## identical(as.matrix(rbind(a,a)),rbind(as.matrix(a),as.matrix(a)))




## test subsetting with/without @other ##
## x <- new("genlight", list(a=1,b=0,c=1), other=list(1:3, letters, data.frame(2:4)))
## pop(x) <- c("pop1","pop1", "pop2")
