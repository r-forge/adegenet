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

    W <- as.matrix(W)

    if(length(seq.names) != nrow(W)){
        stop("inconsistent dimension for W")
    }

    N <- length(seq.names)
    id <- 1:N


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
plotSeqTrack <- function(x, xy, useArrows=TRUE, col=NULL,bg="grey", add=FALSE,...){

    ## CHECKS ##
    if(class(x) != "matrix") stop("x is not a matrix")
    if(nrow(x) != 3) stop("x does not have two rows")
    if(ncol(xy) != 2) stop("xy does not have two columns")
    if(nrow(xy) != ncol(x)) stop("x and xy have inconsistent dimensions")


    ## FIND SEGMENTS COORDS ##
    NA.posi <- which(is.na(x[2,]))
    from <- unlist(x[2,-NA.posi])
    to <- unlist(x[1,-NA.posi])

    x.from <- xy[from,1]
    y.from <- xy[from,2]
    x.to <- xy[to,1]
    y.to <- xy[to,2]

    if(useArrows){
        plotFn <- arrows
    } else {
        plotFn <- segments
    }


    ## FIND THE COLOR FOR EDGES ##
    if(is.null(col)){
        w <- as.numeric(x[3,-NA.posi])
        w <- max(w) - w
        w <- w-min(w)
        w <- 1+ w/max(w) * 99

        opalette <- palette()
        on.exit(palette(opalette))
        palette(heat.colors(100))

        col <- w
    }

    ## DO THE PLOTTING ##
    obg <- par("bg")
    on.exit(par(bg=obg))
    if(!add){
        par(bg=bg)
        plot(xy, type="n")
    }
    plotFn(x.from, y.from, x.to, y.to, col=col,...)
    text(xy,lab=colnames(x), font=2)

    ## RESULT ##
    res <- data.frame(x.from, y.from, x.to, y.to)
    return(invisible(res))
} # end plotSeqTrack








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
