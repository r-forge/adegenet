#############
## seqTrack
#############
seqTrack <- function(seq.names, seq.dates, W, optim=c("min","max"), ...){

    ## CHECKS ##
    optim <- match.arg(optim)
    if(optim=="min"){
        which.optim <- which.min
    } else {
        which.optim <- which.max
    }

    if(length(seq.names) != length(seq.dates)){
        stop("inconsistent length for seq.dates")
    }

    if(is.character(seq.dates)){
        msg <- paste("seq.dates is a character vector; " ,
                     "please convert it as dates using 'as.POSIXct'" ,
                     "\n(making sure dates are given as 'YYYY/MM/DD' or 'YYYY-MM-DD').", sep="")
        stop(msg)
    }

    N <- length(seq.names)
    id <- 1:N

    W <- as.matrix(W)
    ## rename dimensions using id
    colnames(W) <- rownames(W) <- id

    if(length(seq.names) != nrow(W)){
        stop("inconsistent dimension for W")
    }


    ## findAncestor
    findAncestor <- function(idx){ # returns the index of one seq's ancestor
        candid <- which(seq.dates < seq.dates[idx])
        if(length(candid)==0) return(list(ances=NA, weight=NA))
        if(length(candid)==1) return(list(ances=candid, weight=W[idx, candid]))
        ancesId <- as.numeric(names(which.optim(W[idx, candid])))
        return(list(ances=ancesId, weight=W[idx, ancesId])) # Id of the ancestor
    }


    ## BUILD THE OUTPUT ##
    res <- sapply(id, findAncestor)
    res <- rbind(id,res)
    colnames(res) <- seq.names

    return(res)
} # end seqTrack






################
## plotSeqTrack
################
plotSeqTrack <- function(x, xy, useArrows=TRUE, annot=TRUE, dateRange=NULL, dates=NULL,
                         col=NULL, bg="grey", add=FALSE, quiet=TRUE,...){

    ## CHECKS ##
    if(class(x) != "matrix") stop("x is not a matrix")
    if(nrow(x) != 3) stop("x does not have two rows")
    if(ncol(xy) != 2) stop("xy does not have two columns")
    if(nrow(xy) != ncol(x)) stop("x and xy have inconsistent dimensions")


    ## FIND SEGMENTS COORDS ##
    isNA <- is.na(x[2,])
    from <- unlist(x[2,!isNA])
    to <- unlist(x[1,!isNA])

    x.from <- xy[from,1]
    y.from <- xy[from,2]
    x.to <- xy[to,1]
    y.to <- xy[to,2]

    if(useArrows){
        plotFn <- arrows
    } else {
        plotFn <- segments
    }


    ## handle segments/arrows with length 0 ##
    nullLength <- (x.from==x.to) & (y.from==y.to)

    ## FIND THE COLOR FOR EDGES ##
    if(is.null(col)){
        w <- as.numeric(x[3,!isNA])
        w <- max(w) - w
        w <- w-min(w)
        w <- 1+ w/max(w) * 99

        opalette <- palette()
        on.exit(palette(opalette))
        palette(heat.colors(100))

        col <- w
    }

    ## recycle col
    col <- rep(col,length=length(x.from))


    ## HANDLE RANGE OF DATES ##
    if(!is.null(dateRange)){
        if(is.null(dates)){
            stop("dateRange is require without providing dates.")
        }

        if(is.character(dateRange)){
            msg <- paste("dateRange is a character vector; " ,
                     "please convert it as dates using 'as.POSIXct'" ,
                     "\n(making sure dates are given as 'YYYY/MM/DD' or 'YYYY-MM-DD').", sep="")
            stop(msg)
        }

        if(is.character(dates)){
            msg <- paste("dates is a character vector; " ,
                         "please convert it as dates using 'as.POSIXct'" ,
                         "\n(making sure dates are given as 'YYYY/MM/DD' or 'YYYY-MM-DD').", sep="")
            stop(msg)
        }

        if(length(dates) != ncol(x)) stop("length of 'dates' does not match number of rows in 'x'")

        toKeep <- (dates > min(dateRange)) & (dates < max(dateRange))
        if(sum(toKeep)==0) {
            if(!quiet) cat("\nNo item in the specified date range.\n")
            return(NULL)
        }

        ## do the subsetting
        x.from <- x.from[toKeep]
        y.from <- y.from[toKeep]
        x.to <- x.to[toKeep]
        y.to <- y.to[toKeep]
        col <- col[toKeep]
        xy <- xy[toKeep,,drop=FALSE]
        x <- x[,toKeep,drop=FALSE]
    }

    ## DO THE PLOTTING ##
    obg <- par("bg")
    on.exit(par(bg=obg))
    if(!add){
        par(bg=bg)
        plot(xy, type="n")
    }

    suppressWarnings(plotFn(x.from, y.from, x.to, y.to, col=col,...)) # for arrows with length 0
    if(annot) text(xy,lab=colnames(x), font=2)
    if(any(nullLength)) {
        points(x.from[nullLength], y.from[nullLength], cex=2, col=col[nullLength],...)
    }

    ## RESULT ##
    res <- data.frame(x.from, y.from, x.to, y.to)
    return(invisible(res))
} # end plotSeqTrack







##############
## .pAbeforeB
##############
##
## proba that a sequence A, sampled a time Ta,
## actually preceeded B, sampled at time Tb.
##
## - Ta, Tb: sampling times (Ta < Tb)
## - mu0: mutation rate / site / year

.pAbeforeB <- function(Ta, Tb, mu0){
    mu <- mu0/365 # mutation rate / site / day
    t <- 1:1000 # in days
    Pt <- (1-mu)^(t*seq.length)

    ## define the range of surrounding days
    timeSep <- Tb-Ta
    rangeSize <- timeSep + 2000 -1 + 2
    t <- 1:rangeSize # idx of Ta = 1001; idx of Tb = 1001 + timeSep

    ## define probabilities to observe A (resp, B) by days
    Pa <- rep(0, rangeSize)
    Pa[1:2001] <- c(sort(Pt),1,Pt) # proba to observe A for surrounding days
    Pb  <- rep(0, rangeSize)
    Pb[timeSep+(1:2001)] <-  c(sort(Pt),1,Pt)  # proba to observe B for surrounding days

    ## compute proba A preceeded B for a given day
    f1 <- function(day){
        L <- rangeSize
        if(day==L) return(0)
        idSum <- seq(day+1, L, by=1)
        res <- Pa[day] * sum(Pb[idSum])
        return(res)
    } # end f1
}







###############
## .ambigDates
###############
## x: result of seqTrack
## mu0: mutation rate / site / year
## p: threshold probability for a sequence not to mutate
.ambigDates <- function(x, mu0, seq.dates, seq.length, p=0.99){
    mu <- mu0/365 # mutation rate / site / day
    t <- 0:1000 # in days
    Pt <- (1-mu)^(t*seq.length)
    ##plot(Pt,ylim=c(0,1))
    nbDays <- which.min(Pt>p)-1

    isNA <- is.na(unlist(x[2,]))
    date.to <- seq.dates[unlist(x[1,])]
    date.from <- seq.dates[unlist(x[2,])]

    res <- (date.to-date.from) < nbDays*2
    res[is.na(res)] <- FALSE
    return(res)

} # end checkTime






#####################
## optimize.seqTrack
#####################
optimize.seqTrack <- function(seq.names, seq.dates, W, optim=c("min","max"),
                              mu0, seq.length, p=0.99, ...){

    ## CHECKS ##
    optim <- match.arg(optim)
    if(optim=="min"){
        which.optim <- which.min
    } else {
        which.optim <- which.max
    }

    if(length(seq.names) != length(seq.dates)){
        stop("inconsistent length for seq.dates")
    }

    if(is.character(seq.dates)){
        msg <- paste("seq.dates is a character vector; " ,
                     "please convert it as dates using 'as.POSIXct'" ,
                     "\n(making sure dates are given as 'YYYY/MM/DD' or 'YYYY-MM-DD').", sep="")
        stop(msg)
    }

    N <- length(seq.names)
    id <- 1:N

    W <- as.matrix(W)
    ## rename dimensions using id
    colnames(W) <- rownames(W) <- id

    if(length(seq.names) != nrow(W)){
        stop("inconsistent dimension for W")
    }


    ## AUXILIARY FUNCTIONS ##
    valRes <- function(res){
        return(sum(unlist(res[3,]), na.rm=TRUE))
    }

    useNewRes <- function(oldRes,newRes, optim){
        if(optim=="min"){
            res <- valRes(oldRes) > valRes(newRes)
        } else {
            res <- valRes(oldRes) < valRes(newRes)
        }

        return(res)
    }


    ## FIND INITIAL SEQTRACK RESULT ##
    res.ini <- seqTrack(seq.names, seq.dates, W, optim=c("min","max"))


    ## LOOK FOR AMBIGUOUS DATES ##
    isAmbig <- .ambigDates(res.ini, mu0, seq.dates, seq.length, p=0.99)
    if(!any(isAmbig)){
        cat("\nNo ambiguity in dates was found; unique solution returned.\n")
        return(res)
    }


    ## DEFINE NEW TEMPORAL ORDER ##
    permutDate <- function(res, id, newDates){ # id: segment whose vertices are permuted
        idAnc <- unlist(res[2,id])
        idDes <- unlist(res[1,id])
        temp <- newDates[idAnc]
        newDates[idAnc] <- newDates[idDes]
        newDates[idDes] <- temp
        return(newDates)
    }


    ## DO THE OPTIMISATION ##
    res.cur <- res.ini # initialization
    idToPermut <- which(isAmbig) # indices of segments to permut
    dates.cur <- order(seq.dates) # order matching initial one
    permuted <- integer(0) # keeps track of permutated segments (by the index of the descendent)

    for(i in idToPermut){
        dates.new <- permutDate(res.cur, i, dates.cur) # new temporal order
        res.new <- seqTrack(seq.names, newDates, W, optim=c("min","max")) # new result
        if(useNewRes(res.cur, res.new, optim)){ # keep the best result
            res.cur <- res.new
            dates.cur <- dates.new
            permuted <- c(permuted, i)
        }
    }


    ## RESULT ##
    res <- list(res=res.cur, permutID=permuted, newDates=dates.cur)
    nperm <- length(permuted)
    cat("\nBest result found after",nperm,"permutations of ambiguous dates\n")
    return(res)

} # optimize.seqTrack







###############
## seqTrack.ml
###############
seqTrack.ml <- function(aln, seq.dates, , optim=c("min","max"), ...){

} # end seqTrack.ml






#########################
## OLD ADD-HOC VERSION
#########################

## #############
## ## seqTrack
## #############
## seqTrack <- function(seq.names, seq.dates, D, k=5, lag=3, ...){

##     ## CHECKS ##
##     if(length(seq.names) != length(seq.dates)){
##         stop("inconsistent length for seq.dates")
##     }

##     D <- as.matrix(D)

##     if(length(seq.names) != nrow(D)){
##         stop("inconsistent dimension for D")
##     }


##     ## ASSIGNEMENTS ##
##     N <- length(seq.names)
##     STARTPOINTS <- 1:N # global variable, modified thoughout
##     id <- 1:N
##     INPATH <- NULL
##     CURRENTPATH <- list() # current path
##     CURRENTPATHDIST <- list() # current path distances
##     listPaths <- list() # pre-final output
##     listPathsDist <- list() # pre-final output
##     rownames(D) <- id
##     colnames(D) <- id


##     ## AUXILIARY FUNCTIONS ##
##     ## choose a starting sequence, as recent as possible
##     chooseStartPoint <- function(){
##         return(STARTPOINTS[which.max(seq.dates[STARTPOINTS])])
##     }


##     ## find k closest ancestors: returns a list(id=[id of the ancestors], d=[distances to ancestors])
##     findAncestors <- function(currentPoint){
##         candidates <- id[seq.dates < seq.dates[currentPoint]]
##         if(length(candidates)==0) return(NULL) # this will indicate the end of the path
##         res <- D[currentPoint,candidates]
##         if(length(res) <= k) return(list(id=as.integer(names(res)), d=res)) # if less than k candidates
##         res <- sort(res)[1:k] # if more than k, take k closest candidates
##         return(list(id=as.integer(names(res)), d=res))
##     }


##     ## discard paths stemming from the worst ancestors a step ...
##     discardPath <- function(step){
##         tempDist <- sapply(CURRENTPATHDIST, sum)
##         pathFac <- sapply(CURRENTPATH, function(e) e[step])
##         pathLengths <- tapply(tempDist, pathFac, sum)
##         toKeep <- names(pathLengths)[which.min(temp)]
##         toKeep <- as.integer(toKeep)
##         CURRENTPATH <<- CURRENTPATH[pathFac==toKeep]
##         CURRENTPATHDIST <<- CURRENTPATHDIST[pathFac==toKeep]
##         return()
##     } # end discardPath


##     ## update id in INPATH
##     INPATH.up <- function(){
##         temp <- unique(unlist(listPaths))
##         INPATH <<- union(INPATH,temp)
##         return()
##     } # end INPATH.up


##     ## retrieve an already known path
##     findExistingPath <- function(id){
##         if(id %in% INPATH){
##             temp <- sapply(listPaths, function(e) id %in% e) # find in which path
##             res <- listPaths[[which(temp)[1]]] # retrieve the path
##             res <- res[which(res==id):length(res)] # cut the path
##             return(res)
##         } else return(NULL)
##     } # end findExistingPath


##     ## expand one path
##     expandOnePath <- function(onePath){ # onePath has a single path (a vector)
##         newPoints <- findAncestors(onePath[length(onePath)]) # get new points
##         if(is.null(newPoints)) return(onePath)
##         res <- lapply(1:length(newPoints), function(i) c(onePath, newPoints)) # duplicate path and add new pts
##         return(res)
##     } # end expandOnePath


##     ## check if CURRENTPATH should stop (TRUE if yes, FALSE otherwise)
##     checkEndCurrentPath <- function(){
##         isEnded <- function(onePath){
##             oldest <- onePath[length(onePath)]
##             temp <- setdiff(id,oldest)
##             if(all(seq.dates[temp] < seq.dates[oldest])) return(TRUE)
##             return(FALSE)
##         }

##         temp <- sapply(CURRENTPATH, checkEndCurrentPath)
##         if(all(temp)) return(TRUE) # if all paths are ended
##         return(FALSE)
##     } # end checkEndCurrentPath



##     ## FIND ONE PATH ##
##     findPath <- function(id){
##         temp <- findExistingPath(id) # search for an already existing path
##         if(!is.null(temp)) return(temp) # return existing path if needed

##         ## WHILE LOOP ##
##         keepSearching <- TRUE
##         while(keepSearching){
##             for(i in length(CURRENTPATH)){
##                 temp <- CURRENPATH[[1]]

##             }



##             temp <- findExistingPath(id) # search for an already existing path
##             if(!is.null(temp))
##             keepSearching <- !checkEndCurrentPath() # stop searching if all paths are ended





##         }

##         INPATH.up() # update id with known paths

##         path.lengths <- sapply(CURRENTPATHDIST, sum)# compute path length
##         toKeep <- which.min(path.lengths) # keep best (shortest) one
##         res <- list(id=CURRENTPATH[[toKeep]], d=CURRENTPATHDIST[[toKeep]]) # return best path
##         return(res)
##     }


## } # end seqTrack
