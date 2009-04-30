#############
## seqTrack
#############
seqTrack <- function(seq.names, seq.dates, D, k=5, lag=3, ...){

    ## CHECKS ##
    if(length(seq.names) != length(seq.dates)){
        stop("inconsistent length for seq.dates")
    }

    D <- as.matrix(D)

    if(length(seq.names) != nrow(D)){
        stop("inconsistent dimension for D")
    }


    ## ASSIGNEMENTS ##
    N <- length(seq.names)
    STARTPOINTS <- 1:N # global variable, modified thoughout
    id <- 1:N
    INPATH <- NULL
    CURRENTPATH <- list()
    CURRENTPATHDIST <- list()
    rownames(D) <- id
    colnames(D) <- id


    ## UTILITARY FUNCTIONS ##
    chooseStartPoint <- function(){
        return(STARTPOINTS[which.max(seq.dates[STARTPOINTS])])
    }


    findAncestors <- function(currentPoint){
        candidates <- id[seq.dates < seq.dates[currentPoint]]
        if(length(candidates)==0) return(NULL) # this will indicate the end of the path
        res <- D[currentPoint,candidates]
        if(length(res) <= k) return(list(id=as.integer(names(res)), d=res)) # if less than k candidates
        res <- sort(res)[1:k] # if more than k, take k closest candidates
        return(list(id=as.integer(names(res)), d=res))
    }



    discardPath <- function(step){
        tempDist <- sapply(CURRENTPATHDIST, sum)
        pathFac <- sapply(CURRENTPATH, function(e) e[step])
        pathLengths <- tapply(tempDist, pathFac, sum)
        toKeep <- names(pathLengths)[which.min(temp)]
        toKeep <- as.integer(toKeep)
        CURRENTPATH <- CURRENTPATH[pathFac==toKeep]
        CURRENTPATHDIST <- CURRENTPATHDIST[pathFac==toKeep]
        return()
    }


} # end seqTrack
