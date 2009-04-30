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
    CURRENTPATH <- list() # current path
    CURRENTPATHDIST <- list() # current path distances
    listPaths <- list() # pre-final output
    listPathsDist <- list() # pre-final output
    rownames(D) <- id
    colnames(D) <- id


    ## AUXILIARY FUNCTIONS ##
    ## choose a starting sequence, as recent as possible
    chooseStartPoint <- function(){
        return(STARTPOINTS[which.max(seq.dates[STARTPOINTS])])
    }


    ## find k closest ancestors: returns a list(id=[id of the ancestors], d=[distances to ancestors])
    findAncestors <- function(currentPoint){
        candidates <- id[seq.dates < seq.dates[currentPoint]]
        if(length(candidates)==0) return(NULL) # this will indicate the end of the path
        res <- D[currentPoint,candidates]
        if(length(res) <= k) return(list(id=as.integer(names(res)), d=res)) # if less than k candidates
        res <- sort(res)[1:k] # if more than k, take k closest candidates
        return(list(id=as.integer(names(res)), d=res))
    }


    ## discard paths stemming from the worst ancestors a step ...
    discardPath <- function(step){
        tempDist <- sapply(CURRENTPATHDIST, sum)
        pathFac <- sapply(CURRENTPATH, function(e) e[step])
        pathLengths <- tapply(tempDist, pathFac, sum)
        toKeep <- names(pathLengths)[which.min(temp)]
        toKeep <- as.integer(toKeep)
        CURRENTPATH <<- CURRENTPATH[pathFac==toKeep]
        CURRENTPATHDIST <<- CURRENTPATHDIST[pathFac==toKeep]
        return()
    } # end discardPath


    ## update id in INPATH
    INPATH.up <- function(){
        temp <- unique(unlist(listPaths))
        INPATH <<- union(INPATH,temp)
        return()
    } # end INPATH.up


    ## retrieve an already known path
    findExistingPath <- function(id){
        if(id %in% INPATH){
            temp <- sapply(listPaths, function(e) id %in% e) # find in which path
            res <- listPaths[[which(temp)[1]]] # retrieve the path
            res <- res[which(res==id):length(res)] # cut the path
            return(res)
        } else return(NULL)
    } # end findExistingPath


    ## expand one path
    expandOnePath <- function(onePath){ # onePath has a single path (a vector)
        newPoints <- findAncestors(onePath[length(onePath)]) # get new points
        if(is.null(newPoints)) return(onePath)
        res <- lapply(1:length(newPoints), function(i) c(onePath, newPoints)) # duplicate path and add new pts
        return(res)
    } # end expandOnePath


    ## check if CURRENTPATH should stop (TRUE if yes, FALSE otherwise)
    checkEndCurrentPath <- function(){
        isEnded <- function(onePath){
            oldest <- onePath[length(onePath)]
            temp <- setdiff(id,oldest)
            if(all(seq.dates[temp] < seq.dates[oldest])) return(TRUE)
            return(FALSE)
        }

        temp <- sapply(CURRENTPATH, checkEndCurrentPath)
        if(all(temp)) return(TRUE) # if all paths are ended
        return(FALSE)
    } # end checkEndCurrentPath



    ## FIND ONE PATH ##
    findPath <- function(id){
        temp <- findExistingPath(id) # search for an already existing path
        if(!is.null(temp)) return(temp) # return existing path if needed

        ## WHILE LOOP ##
        keepSearching <- TRUE
        while(keepSearching){
            for(i in length(CURRENTPATH)){
                temp <- CURRENPATH[[1]]

            }



            temp <- findExistingPath(id) # search for an already existing path
            if(!is.null(temp))
            keepSearching <- !checkEndCurrentPath() # stop searching if all paths are ended





        }

        INPATH.up() # update id with known paths

        path.lengths <- sapply(CURRENTPATHDIST, sum)# compute path length
        toKeep <- which.min(path.lengths) # keep best (shortest) one
        res <- list(id=CURRENTPATH[[toKeep]], d=CURRENTPATHDIST[[toKeep]]) # return best path
        return(res)
    }


} # end seqTrack
