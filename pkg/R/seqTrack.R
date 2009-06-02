#############
## seqTrack
#############
seqTrack <- function(seq.names, seq.dates, W, optim=c("min","max"),
                     proxMat=NULL, ...){

    ## CHECKS ##
    optim <- match.arg(optim)
    if(optim=="min"){
        optim <- min
        which.optim <- which.min
    } else {
        optim <- max
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

    W <- as.matrix(W)

    if(!is.null(proxMat) && !identical(dim(proxMat),dim(W))) {
        stop("proxMat is provided but its dimensions are inconsistent with that of W")
    }

    N <- length(seq.names)
    id <- 1:N

    seq.dates <- as.POSIXct(round.POSIXt(seq.dates,units="days")) # round dates to the day

    W <- as.matrix(W)
    ## rename dimensions using id
    colnames(W) <- rownames(W) <- id
    if(!is.null(proxMat)){
        colnames(proxMat) <- rownames(proxMat) <- id
    }

    if(length(seq.names) != nrow(W)){
        stop("inconsistent dimension for W")
    }



    ## AUXILIARY FUNCTIONS ##
    ## test equality in floats
    test.equal <- function(val,vec){
        return(abs(val-vec) < 1e-12)
    }


    ## return the names of optimal value(s) in a named vector
    which.is.optim <- function(vec){
        res <- names(vec)[test.equal(optim(vec), vec)]
        return(res)
    }


    ## select among different possible ancestors
    selAmongAncestors <- function(idx,ances){
        ## Choose the most otherwise connected ancestor, given proxMat
        if(!is.null(proxMat)){ # if we've got no other info
            toKeep <- test.equal(max(proxMat[idx,ances]), proxMat[idx,ances])
            ances <- ances[toKeep]
        }

        ## If several ancestors remain, take the oldest.
        if(length(ances)>1){
            ances <- ances[which.min(seq.dates[ances])]
        }

        return(ances)
    }


    ## findAncestor
    findAncestor <- function(idx){ # returns the index of one seq's ancestor
        candid <- which(seq.dates < seq.dates[idx])
        if(length(candid)==0) return(list(ances=NA, weight=NA))
        if(length(candid)==1) return(list(ances=candid, weight=W[idx, candid]))
        ancesId <- as.numeric(which.is.optim(W[idx, candid]))
        if(length(ancesId)>1) {
            ancesId <- selAmongAncestors(idx,ancesId) # handle several 'best' ancestors
        }
        return(list(ances=ancesId, weight=W[idx, ancesId])) # Id of the ancestor
    }


    ## BUILD THE OUTPUT ##
    res <- sapply(id, findAncestor)
    res <- data.frame(ances=unlist(res[1,]), weight=unlist(res[2,]))
    ances.date <- seq.dates[res[,1]]
    res <- cbind.data.frame(id,res, date=seq.dates, ances.date)
    rownames(res) <- seq.names

    return(res)
} # end seqTrack






################
## plotSeqTrack
################
plotSeqTrack <- function(x, xy, useArrows=TRUE, annot=TRUE, dateRange=NULL,
                         col=NULL, bg="grey", add=FALSE, quiet=TRUE,
                         showAmbiguous=TRUE, mu0=NULL, seq.length=NULL, prob=0.75,
                         plot=TRUE,...){

    ## CHECKS ##
    if(class(x) != "data.frame") stop("x is not a data.frame")
    if(ncol(x) != 5) stop("x does not have five columns")
    if(ncol(xy) != 2) stop("xy does not have two columns")
    if(nrow(xy) != nrow(x)) stop("x and xy have inconsistent dimensions")
    if(showAmbiguous & (is.null(mu0) | is.null(seq.length)) ){
        stop("showAmbiguous is TRUE, but mu0 and seq.length are not all provided.")
    }

    isAmbig <- NULL


    ## SUBSET DATA (REMOVE NAs) ##
    isNA <- is.na(x[,2])
    x <- x[!isNA,,drop=FALSE]
    xy.all <- xy ## used to retrieve all coordinates
    xy <- xy[!isNA,,drop=FALSE]


    ## FIND AMBIGUOUS TEMPORAL ORDERING ##
    if(showAmbiguous){
        temp <- .pAbeforeB(x$ances.date, x$date, mu0, seq.length)
        isAmbig <- temp < prob
    }


    ## FIND SEGMENTS COORDS ##
    from <- unlist(x[,2])
    to <- unlist(x[,1])

    x.from <- xy.all[from,1]
    y.from <- xy.all[from,2]
    x.to <- xy.all[to,1]
    y.to <- xy.all[to,2]


    ## FIND THE COLOR FOR EDGES ##
    if(is.null(col)){
        if(is.null(mu0) & is.null(seq.length)) {
            col <- "black"
        } else {
            ##  w <- .pAbeforeB(x$ances.date, x$date, mu0, seq.length, 200) # beware, lots of steps take time
            ##             isAmbig <-  w < prob
            ##             w[w<.5] <- .5
            ##             w <- (1 - w)
            ##             w <- w - min(w) # rescale to 0-1
            ##             w <- 100*w/max(w)  # rescale to 0-100
            ##             w[w<1] <- 1

            opalette <- palette()
            on.exit(palette(opalette))

            w <- as.integer(round(x$weight))
            col <- rep("yellow", length(w))
            col[w <= 1] <- "orange"
            col[w < 1] <- "red"
        }
    }

    ## THIS WAS USED WHEN COLOR REPRESENTED THE NUMBER OF MUTATIONS ##
    ##  if(is.null(col)){
    ##         w <- as.numeric(x[,3])
    ##         w <- max(w) - w
    ##         w <- w-min(w)
    ##         w <- 1+ w/max(w) * 99

    ##         opalette <- palette()
    ##         on.exit(palette(opalette))
    ##         palette(heat.colors(100))

    ##         col <- w
    ##     }


    ## recycle col
    col <- rep(col,length=length(x.from))


    ## HANDLE RANGE OF DATES ##
    if(!is.null(dateRange)){

        if(is.character(dateRange)){
            msg <- paste("dateRange is a character vector; " ,
                         "please convert it as dates using 'as.POSIXct'" ,
                         "\n(making sure dates are given as 'YYYY/MM/DD' or 'YYYY-MM-DD').", sep="")
            stop(msg)
        }

        dates <- x$date
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
        x <- x[toKeep,,drop=FALSE]
        if(!is.null(isAmbig)) {
            isAmbig <- isAmbig[toKeep]
        }
    }

    ## DO THE PLOTTING ##
    if(plot){
        obg <- par("bg")
        on.exit(par(bg=obg))
        if(!add){
            par(bg=bg)
            plot(xy, type="n")
        }
    }

    ## ARROWS
    if(useArrows & plot){
        if(showAmbiguous & any(isAmbig)){ # plot arrows & segments
            suppressWarnings(arrows(x.from[!isAmbig], y.from[!isAmbig],
                                    x.to[!isAmbig], y.to[!isAmbig], col=col[!isAmbig], angle=15, ...))
            segments(x.from[isAmbig], y.from[isAmbig],
                     x.to[isAmbig], y.to[isAmbig], col=col,...)
        } else{ # plot all arrows
            suppressWarnings(arrows(x.from, y.from, x.to, y.to, col=col, angle=15, ...))
        }
    } else{
        ## SEGMENTS
        if(plot) segments(x.from, y.from, x.to, y.to, col=col,...)
    }


    if(annot & plot) text(xy,lab=rownames(x), font=2)


    ## handle segments/arrows with length 0 ##
    nullLength <- (x.from==x.to) & (y.from==y.to)

    if(any(nullLength) & plot) {
        points(x.from[nullLength], y.from[nullLength], col=col[nullLength], cex=2, pch=20, ...)
        sunflowerplot(x.from[nullLength], y.from[nullLength], seg.lwd=2, size=1/6,
                      col=col[nullLength], seg.col=col[nullLength], add=TRUE, ...)
    }


    ## RESULT ##
    res <- data.frame(x.from, y.from, x.to, y.to, col=col)
    if(!is.null(isAmbig)) {
        res <- cbind.data.frame(res, isAmbig)
    }
    return(invisible(res))
} # end plotSeqTrack






#############
## .dTimeSeq
#############
.dTimeSeq <- function(mu0, L, maxNbDays=100){
    mu <- mu0/365 # mutation rate / site / day
    t <- 0:maxNbDays # in days added / substracted
    Pt <- (1-mu)^(t*L)
    t <- c(-rev(t[-1]), t)
    Pt <- c(rev(Pt[-1]), Pt)
    return(list(t, Pt))
}


#############
## .rTimeSeq
#############
.rTimeSeq <- function(n, mu0, L, maxNbDays=100){
    temp <- .dTimeSeq(mu0, L, maxNbDays)
    p <- temp[[2]]/sum(temp[[2]])
    res <- sample(temp[[1]], size=n, replace=TRUE, prob=p)
    return(res)
}



#################
## .rUnifTimeSeq
#################
.rUnifTimeSeq <- function(n, dateMin, dateMax){
    rangeSize <-  as.integer(difftime(dateMax,dateMin, units="days"))
    nbDays <- round(runif(n, min=0, max=rangeSize))
    res <- dateMin + nbDays*3600*24
    res <- as.POSIXct(round(res, units="days"))
    return(res)
}





###############
## .ambigDates
###############
## x: result of seqTrack
## mu0: mutation rate / site / year
## ## p: threshold probability for a sequence not to mutate
## .ambigDates <- function(x, mu0, seq.length, p=0.99){
##     mu <- mu0/365 # mutation rate / site / day
##     t <- 0:1000 # in days
##     Pt <- (1-mu)^(t*seq.length)
##     ##plot(Pt,ylim=c(0,1))
##     nbDays <- which.min(Pt>p)-1

##     isNA <- is.na(x[,2])
##     date.to <- x$date
##     date.from <- x$ances.date

##     res <- (date.to-date.from) < nbDays*2
##     res[is.na(res)] <- FALSE
##     return(res)

## } # end checkTime






##############
## .pAbeforeB
##############
##
## TODO: replace mu0 and L by muA, muB, LA, LB
##
.pAbeforeB <- function(dateA, dateB, mu0, L, maxNbDays=100){
    temp <- .dTimeSeq(mu0, L, maxNbDays)
    days <- temp[[1]]
    p <- temp[[2]]/sum(temp[[2]]) # proba distribution

    nbDaysDiff <- as.integer(round(difftime(dateA,dateB,units="days"))) # dateA - dateB, in days
    daysA <- days
    daysB <- days - nbDaysDiff

    f1 <- function(i){ # proba A before B for one day
        idx <- daysB>daysA[i]
        res <- p[i] * sum(p[idx])
        return(res)
    }

    res <- sapply(1:length(p), f1) # proba for all days
    res <- sum(res) # sum
    return(res)
}

.pAbeforeB <- Vectorize(.pAbeforeB, vectorize.args=c("dateA","dateB")) ## end .pAbeforeB








#####################
## optimize.seqTrack
#####################
##
## TODO:
## 1) Change the output to retain xxx simulations
## 2) VECTORIZE mu0 and seq.length, recycle if needed with a warning
## 3) uncomment, adapt, and test code for missing data
##
optimize.seqTrack <- function(nsim, seq.names, seq.dates, W, thres, optim=c("min","max"),
                              proxMat=NULL, mu0, seq.length, rMissDate=.rUnifTimeSeq, ...){

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

    isMissDate <- is.na(seq.dates)


    N <- length(seq.names)
    id <- 1:N

    seq.dates <- as.POSIXct(round.POSIXt(seq.dates,units="days")) # round dates to the day

    W <- as.matrix(W)

    if(!is.null(proxMat) && !identical(dim(proxMat),dim(W))) {
        stop("proxMat is provided but its dimensions are inconsistent with that of W")
    }

    ## rename dimensions using id
    colnames(W) <- rownames(W) <- id

    if(length(seq.names) != nrow(W)){
        stop("inconsistent dimension for W")
    }


    ## AUXILIARY FUNCTIONS ##

    ## to compare results -> returns a list of length two: logical, and the value of the res
    val.res <- function(res){
        return(sum(res$weight, na.rm=TRUE))
    }


    ## DO THE OPTIMISATION ##
    RANGE.DATES <- as.integer(round(diff(range(seq.dates, na.rm=TRUE)))) # time window of the sample, in days
    NB.DATES.TO.SIM <- sum(!isMissDate)


    ## for loop is not to slow for < 1e6 rep
    ## and allows not to handle huge objects
    ## (which would grow exponentially)

    ##    res.best <- res.ini # initialization


    ## DEFINE OUTPUTS ##
    ances <- integer(0)
    date <- character(0)
    ances.date <- character(0)
    valRes <- numeric(0)


    ## DEFAULT CASE: NO MISSING DATES
    if(!any(isMissDate)){
        for(i in 1:nsim){
            myDates <- seq.dates +
                .rTimeSeq(n=NB.DATES.TO.SIM, mu0=mu0, L=seq.length, maxNbDays=RANGE.DATES)*24*3600
            myDates <- as.POSIXct(round(myDates, units="days"))
            res.new <- seqTrack(seq.names=seq.names, seq.dates=myDates, W=W, optim=optim, proxMat=proxMat, ...)
            temp <- val.res(res.new)
            if(ifelse(optim=="min", temp < thres, temp > thres)){
                ances <- cbind(ances, res.new$ances)
                date <- cbind(date, res.new$date)
                ances.date <- cbind(ances.date, res.new$ances.date)
                valRes <- c(valRes, temp)
            }
        }
    }


    ##  ## OTHER CASE: HANDLE MISSING DATES
    ##     if(any(isMissDate)){

    ##         ## Handle distribution and its parameters ##
    ##         argList <- list(...)

    ##         if(is.null(argList$dateMin) & identical(rMissDate, .rUnifTimeSeq)){ # earliest date
    ##             argList$dateMin <- min(seq.dates,na.rm=TRUE)
    ##         } else {
    ##             argList$dateMin[is.na(argList$dateMin)] <- min(seq.dates,na.rm=TRUE)
    ##         }
    ##         if(is.null(argList$dateMax) & identical(rMissDate, .rUnifTimeSeq)){ # latest date
    ##             argList$dateMax <- max(seq.dates,na.rm=TRUE)
    ##         } else {
    ##             argList$dateMax[is.na(argList$dateMax)] <- max(seq.dates,na.rm=TRUE)
    ##         }

    ##         argList$n <- sum(isMissDate)

    ##         ## Make simulations ##
    ##         for(i in 1:nsim){
    ##             myDates <- seq.dates
    ##             ## distribution for available dates
    ##             myDates[!isMissDate] <- myDates[!isMissDate] +
    ##                 .rTimeSeq(n=NB.DATES.TO.SIM, mu0=mu0, L=seq.length, maxNbDays=2*RANGE.DATES)
    ##             ## distribution for missing dates
    ##             myDates[isMissDate] <- do.call(rMissDate, argList)

    ##             res.new <- seqTrack(seq.names=seq.names, seq.dates=myDates, W=W, optim=optim, proxMat=proxMat, ...)

    ##             valRes[i] <- sum(res.new$weight,na.rm=TRUE)
    ##             if(use.new.res(res.best, res.new)){
    ##                 res.best <- res.new
    ##             }
    ##         }
    ##     }


    ## RESULT ##
    res <- list(best=res.best, valsim=valRes)
    return(res)

} # optimize.seqTrack







###############
## seqTrack.ml
###############
seqTrack.ml <- function(aln, seq.dates, optim=c("min","max"), ...){

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
