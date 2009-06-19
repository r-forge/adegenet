###########
# generics
###########
seqTrack <- function(...){
    UseMethod("seqTrack")
}


optimize.seqTrack <- function(...){
    UseMethod("optimize.seqTrack")
}


get.likelihood <- function(...){
    UseMethod("get.likelihood")
}

get.likelihood.seqTrack.default <- function(...){
    cat("Method not implemented.")
    return()
}



#############
## seqTrack
#############
seqTrack.default <- function(x, x.names, x.dates, optim=c("min","max"),
                     prox.mat=NULL, ...){

    ## CHECKS ##
    optim <- match.arg(optim)
    if(optim=="min"){
        optim <- min
        which.optim <- which.min
    } else {
        optim <- max
        which.optim <- which.max
    }

    if(length(x.names) != length(x.dates)){
        stop("inconsistent length for x.dates")
    }

    if(is.character(x.dates)){
        msg <- paste("x.dates is a character vector; " ,
                     "please convert it as dates using 'as.POSIXct'" ,
                     "\n(making sure dates are given as 'YYYY/MM/DD' or 'YYYY-MM-DD').", sep="")
        stop(msg)
    }

    x <- as.matrix(x)

    if(!is.null(prox.mat) && !identical(dim(prox.mat),dim(x))) {
        stop("prox.mat is provided but its dimensions are inconsistent with that of x")
    }

    N <- length(x.names)
    id <- 1:N

    x.dates <- as.POSIXct(round.POSIXt(x.dates,units="days")) # round dates to the day

    x <- as.matrix(x)
    ## rename dimensions using id
    colnames(x) <- rownames(x) <- id
    if(!is.null(prox.mat)){
        colnames(prox.mat) <- rownames(prox.mat) <- id
    }

    if(length(x.names) != nrow(x)){
        stop("inconsistent dimension for x")
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
        ## Choose the most otherwise connected ancestor, given prox.mat
        if(!is.null(prox.mat)){ # if we've got no other info
            toKeep <- test.equal(max(prox.mat[idx,ances]), prox.mat[idx,ances])
            ances <- ances[toKeep]
        }

        ## If several ancestors remain, take the oldest.
        if(length(ances)>1){
            ances <- ances[which.min(x.dates[ances])]
        }

        return(ances)
    }


    ## findAncestor
    findAncestor <- function(idx){ # returns the index of one seq's ancestor
        candid <- which(x.dates < x.dates[idx])
        if(length(candid)==0) return(list(ances=NA, weight=NA))
        if(length(candid)==1) return(list(ances=candid, weight=x[idx, candid]))
        ancesId <- as.numeric(which.is.optim(x[idx, candid]))
        if(length(ancesId)>1) {
            ancesId <- selAmongAncestors(idx,ancesId) # handle several 'best' ancestors
        }
        return(list(ances=ancesId, weight=x[idx, ancesId])) # Id of the ancestor
    }


    ## BUILD THE OUTPUT ##
    res <- sapply(id, findAncestor)
    res <- data.frame(ances=unlist(res[1,]), weight=unlist(res[2,]))
    ances.date <- x.dates[res[,1]]
    res <- cbind.data.frame(id,res, date=x.dates, ances.date)
    rownames(res) <- x.names

    class(res) <- c("seqTrack","data.frame")

    return(res)
} # end seqTrack






################
## plotSeqTrack
################
plotSeqTrack <- function(x, xy, useArrows=TRUE, annot=TRUE, dateRange=NULL,
                         col=NULL, bg="grey", add=FALSE, quiet=FALSE,
                         showAmbiguous=TRUE, mu0=NULL, chr.length=NULL, prob=0.75,
                         plot=TRUE,...){

    ## CHECKS ##
    if(!inherits(x,"seqTrack")) stop("x is not a seqTrack object")
    ##if(ncol(x) != 5) stop("x does not have five columns")
    if(ncol(xy) != 2) stop("xy does not have two columns")
    if(nrow(xy) != nrow(x)) stop("x and xy have inconsistent dimensions")
    if(showAmbiguous & (is.null(mu0) | is.null(chr.length)) ){
        stop("showAmbiguous is TRUE, but mu0 and chr.length are not all provided.")
    }

    isAmbig <- NULL


    ## SUBSET DATA (REMOVE NAs) ##
    isNA <- is.na(x[,2])
    x <- x[!isNA,,drop=FALSE]
    xy.all <- xy ## used to retrieve all coordinates
    xy <- xy[!isNA,,drop=FALSE]


    ## FIND AMBIGUOUS TEMPORAL ORDERING ##
    if(showAmbiguous){
        temp <- .pAbeforeB(x$ances.date, x$date, mu0, chr.length)
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
        if(is.null(mu0) & is.null(chr.length)) {
            col <- "black"
        } else {
            ##  w <- .pAbeforeB(x$ances.date, x$date, mu0, chr.length, 200) # beware, lots of steps take time
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
##
## mu0 and L are vectors, having one value per segment/chromosome
.dTimeSeq <- function(mu0, L, maxNbDays=100){
    mu <- mu0/365 # mutation rates / site / day
    t <- 0:maxNbDays # in days added / substracted
    temp <- sapply((1-mu)^L, function(x) x^t  )
    Pt <- apply(temp,1,prod)
    t <- c(-rev(t[-1]), t)
    Pt <- c(rev(Pt[-1]), Pt)
    return(list(t, Pt))
}


#############
## .rTimeSeq
#############
##
## mu0 and L are vectors, having one value per segment/chromosome
##
## this returns nb days
.rTimeSeq <- function(n, mu0, L, maxNbDays=100){
    temp <- .dTimeSeq(mu0, L, maxNbDays)
    res <- sample(temp[[1]], size=n, replace=TRUE, prob= temp[[2]]/sum(temp[[2]]))
    return(res)
}



#################
## .rUnifDate
#################
##
## this returns random uniform dates in a given range
##
.rUnifDate <- function(n, dateMin, dateMax, ...){
    rangeSize <-  as.integer(difftime(dateMax,dateMin, units="days"))
    nbDays <- round(runif(n, min=0, max=rangeSize))
    res <- dateMin + nbDays*3600*24
    res <- as.POSIXct(round(res, units="days"))
    return(res)
}



#################
## .rNormTimeSeq
#################
##
## this returns nb of days
.rNormTimeSeq <- function(n, mean, sd, ...){
    res <- round(rnorm(n, mean=mean, sd=sd))
    return(res)
}



#################
## .rSampTimeSeq
#################
##
## this returns nb of days
.rSampTime <- function(n,...){
    res <- round(rnorm(n*2, -2))
    res <- res[res < 0 & res > -7][1:n]
    return(res)
}




###############
## .ambigDates
###############
## x: result of seqTrack
## mu0: mutation rate / site / year
## ## p: threshold probability for a sequence not to mutate
## .ambigDates <- function(x, mu0, chr.length, p=0.99){
##     mu <- mu0/365 # mutation rate / site / day
##     t <- 0:1000 # in days
##     Pt <- (1-mu)^(t*chr.length)
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
.pAbeforeB <- function(dateA, dateB, muA, muB, LA, LB, maxNbDays=100){
    ## proba dist for A
    tempA <- .dTimeSeq(muA, LA, maxNbDays)
    days <- tempA[[1]]
    pA <- tempA[[2]]/sum(tempA[[2]]) # proba distribution

    ## proba dist for B
    tempB <- .dTimeSeq(muB, LB, maxNbDays)
    pB <- tempB[[2]]/sum(tempB[[2]]) # proba distribution

    ## days for A and B
    nbDaysDiff <- as.integer(round(difftime(dateA,dateB,units="days"))) # dateA - dateB, in days
    daysA <- days
    daysB <- days - nbDaysDiff

    f1 <- function(i){ # proba A before B for one day
        idx <- daysB > daysA[i]
        return(pA[i] * sum(pB[idx]))
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
## 1) Change the output to retain xxx simulations | ok. -- done.
## 2) VECTORIZE mu0 and chr.length, recycle if needed with a warning
## 3) uncomment, adapt, and test code for missing data
##
optimize.seqTrack.default <- function(x, x.names, x.dates, typed.chr=NULL, mu0=NULL, chr.length=NULL,
                                      thres=0.2, optim=c("min","max"), prox.mat=NULL, nstep=10, step.size=1e3,
                                      rDate=.rTimeSeq, arg.rDate=NULL, rMissDate=.rUnifDate, ...){


    ## CHECKS ##
    optim <- match.arg(optim)
    if(optim=="min"){
        which.optim <- which.min
    } else {
        which.optim <- which.max
    }

    if(length(x.names) != length(x.dates)){
        stop("inconsistent length for x.dates")
    }

    if(is.character(x.dates)){
        msg <- paste("x.dates is a character vector; " ,
                     "please convert it as dates using 'as.POSIXct'" ,
                     "\n(making sure dates are given as 'YYYY/MM/DD' or 'YYYY-MM-DD').", sep="")
        stop(msg)
    }

    isMissDate <- is.na(x.dates)

    if(!identical(rDate, .rTimeSeq)){
        if(is.null(arg.rDate)){
            warning("Specific time distribution specified without arguments.")
            arg.rDate <- list(n=step.size)
        } else {
            if(!is.list(arg.rDate)) stop("If provided, arg.rDate must be a list.")
            if(!is.null(arg.rDate$n)) {
                warning("arg.rDate$n is provided, but will be replaced by step.size.")
            }
            arg.rDate$n <- step.size
        }
    }


    N <- length(x.names)
    id <- 1:N
    ## if(length(mu0) < N) { # recycle mu0
    ##         mu0 <- rep(mu0, length=N)
    ##     }
    ##     if(length(chr.length) < N) {# recycle chr.length
    ##         chr.length <- rep(chr.length, length=N)
    ##     }


    ## handle typed.chr, mu0, chr.length
    if(identical(rDate, .rTimeSeq)){
        if(is.null(typed.chr)|is.null(mu0)|is.null(chr.length)){
            stop("typed.chr, mu0, and chr.length must be provided if rDate is .rTimeSeq")
        }

        if(!is.list(typed.chr)) {
            stop("typed.chr must be a list")
        }
        if(length(typed.chr)!=N) {
            stop("typed.chr has an inconsistent length")
        }

        if(is.null(names(mu0))) stop("mu0 has no names")
        if(is.null(names(chr.length))) stop("chr.length has no names")
        if(any(mu0 > 1)) stop("mu0 has values > 1")
        if(any(mu0 < 0)) stop("mu0 has negative values")

        if(!identical(names(mu0) , names(chr.length))) stop("Names of mu0 and chr.length differ.")
        if(any(!unique(unlist(typed.chr)) %in% names(mu0))) {
            stop("Some chromosomes indicated in typed.chr are not in mu0.")
        }

        list.mu0 <- lapply(typed.chr, function(e) mu0[e])
        list.chr.length <- lapply(typed.chr, function(e) chr.length[e])
    }

    x.dates <- as.POSIXct(round.POSIXt(x.dates,units="days")) # round dates to the day

    x <- as.matrix(x)

    if(!is.null(prox.mat) && !identical(dim(prox.mat),dim(x))) {
        stop("prox.mat is provided but its dimensions are inconsistent with that of x")
    }


    ## rename dimensions using id
    colnames(x) <- rownames(x) <- id

    if(length(x.names) != nrow(x)){
        stop("inconsistent dimension for x")
    }


    ## SET THRESHOLD IF NEEDED ## ## NO LONGER USED
    ##  if(is.null(thres)){
    ##         thres <- sum(seqTrack(x.names=x.names, x.dates=x.dates, W=W,
    ##                                     optim=optim, prox.mat=prox.mat, ...)$weight, na.rm=TRUE)
    ##     }


    ## AUXILIARY FUNCTIONS ##

    ## to compare results -> returns a list of length two: logical, and the value of the res
    val.res <- function(res){
        return(sum(res$weight, na.rm=TRUE))
    }


    ## DO THE OPTIMISATION ##
    RANGE.DATES <- as.integer(round(diff(range(x.dates, na.rm=TRUE)))) # time window of the sample, in days
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
        ## dates initialisation, taken from initial prior
        ## If dates distrib is .rTimeSeq
        if(identical(rDate, .rTimeSeq)){
            newDates <- sapply(1:N, function(i)
                               rDate(n=step.size, mu0=list.mu0[[i]], L=list.chr.length[[i]],
                                     maxNbDays=RANGE.DATES))
        } else { ## Else, any other distrib with free arguements
            newDates <- sapply(1:N, function(i) do.call(rDate, arg.rDate))
        }

        newDates <- t(newDates)*24*3600 + x.dates

        ## >> one step of 'step.size' simulations, all with same prior << ##
        for(i in 1:nstep){
            ## >> each step contains 'step.size' iterations << ##
            for(j in 1:step.size){
                myDates <- as.POSIXct(newDates[,j])

                res.new <- seqTrack(x, x.names=x.names, x.dates=myDates,
                                    optim=optim, prox.mat=prox.mat, ...)

                ##ances <- cbind(ances, res.new$ances) # not needed now
                date <- cbind(date, as.character(res.new$date))
                ##ances.date <- cbind(ances.date, as.character(res.new$ances.date)) # not needed now
                valRes <- c(valRes, val.res(res.new))
                ##}
            } # end for j

            ## retain a given % (thres) of the dates ##
            toKeep <- valRes <= quantile(valRes, thres) ## NOT WORKING FOR optim==max !!!
            valRes <- valRes[toKeep]

            date <- date[,toKeep,drop=FALSE] # retained posterior

            ## DEBUGING ##
            ## cat("\ntoKeep:\n")
            ##             print(toKeep)
            ##             cat("\nhead date (posterior):\n")
            ##             print(head(date))
            ## END DEBUGING ##

            newDates <- apply(date, 1, function(vec)
                              sample(vec, size=step.size, replace=TRUE)) # new prior
            newDates <- t(newDates)

            ## stop if all dates are fixed
            if(all(apply(newDates, 1, function(r) length(unique(r))==1))){
                cat("\nConvergence reached at step",i,"\n")
                break # stop the algorithm
            }

            ## re-initialize posterior distributions
            if(i<nstep){
                ## ances <- integer(0) # not needed now
                date <- character(0)
                ## ances.date <- character(0) # not needed now
                valRes <- numeric(0)
            } # end if

        } # end for i

        ##  ## dates: new prior taken from obtained posterior
        ##             if(length(valRes)==0) { # if no simul are retained
        ##                 warning(paste("No simulation was retained at the given threshold at step",i))
        ##             } else {
        ##  if(optim=="min"){ # define weights for further samplings
        ##                     w <- max(valRes,na.rm=TRUE) - valRes
        ##                     w <- w/sum(w)
        ##                 } else {
        ##                     w <- valRes
        ##                     w <- w/sum(w)
        ##                 }

        ## newDates <- apply(date, 1, function(vec) #  used a weighted sampling
        ##                                  sample(vec, size=step.size, replace=TRUE, prob=w))


    } # end if(!any(isMissDate))


    ##  ## OTHER CASE: HANDLE MISSING DATES
    ##     if(any(isMissDate)){

    ##         ## Handle distribution and its parameters ##
    ##         argList <- list(...)

    ##         if(is.null(argList$dateMin) & identical(rMissDate, .rUnifDate)){ # earliest date
    ##             argList$dateMin <- min(x.dates,na.rm=TRUE)
    ##         } else {
    ##             argList$dateMin[is.na(argList$dateMin)] <- min(x.dates,na.rm=TRUE)
    ##         }
    ##         if(is.null(argList$dateMax) & identical(rMissDate, .rUnifDate)){ # latest date
    ##             argList$dateMax <- max(x.dates,na.rm=TRUE)
    ##         } else {
    ##             argList$dateMax[is.na(argList$dateMax)] <- max(x.dates,na.rm=TRUE)
    ##         }

    ##         argList$n <- sum(isMissDate)

    ##         ## Make simulations ##
    ##         for(i in 1:nstep){
    ##             myDates <- x.dates
    ##             ## distribution for available dates
    ##             myDates[!isMissDate] <- myDates[!isMissDate] +
    ##                 rDate(n=NB.DATES.TO.SIM, mu0=mu0, L=chr.length, maxNbDays=2*RANGE.DATES)
    ##             ## distribution for missing dates
    ##             myDates[isMissDate] <- do.call(rMissDate, argList)

    ##             res.new <- seqTrack(x.names=x.names, x.dates=myDates, W=W, optim=optim, prox.mat=prox.mat, ...)

    ##             valRes[i] <- sum(res.new$weight,na.rm=TRUE)
    ##             if(use.new.res(res.best, res.new)){
    ##                 res.best <- res.new
    ##             }
    ##         }
    ##     }


    ## RESULT ##

    ## reconstruct the result with new dates
    res <- lapply(1:ncol(date), function(i)
                   seqTrack(x=x, x.names=x.names, x.dates=as.POSIXct(date[,i]),
                                    optim=optim, prox.mat=prox.mat, ...))
    ances <- data.frame(lapply(res, function(e) e$ances))
    ances <- matrix(as.integer(unlist(ances)), nrow=nrow(ances))

    ances.date <- data.frame(lapply(res, function(e) as.character(e$ances.date)))
    ances.date <- matrix(as.character(unlist(ances.date)), nrow=nrow(ances.date))

    res <- list(ances=ances, date=date, ances.date=ances.date, valsim=valRes)
    return(res)

} # optimize.seqTrack






#################
## get.result.by
#################
get.result.by <- function(x, ...){
    dat <- list(...)
    if(length(dat)==0) return(x)


    ## DEFINE NEW VALUES ##

    convertElem <- function(e){
        if(class(e)=="DNAbin") {
            e <- as.matrix(e)
        }
        ori.dim <- dim(e)
        e <- as.character(e)
        dim(e) <- ori.dim
        return(e)
    }


    dat <- lapply(dat,convertElem)
    dat <- as.matrix(data.frame(dat))

    newval <- apply(dat, 1, function(vec) paste(vec, collapse=""))
    newval <- unclass(factor(newval))
    newlev <- levels(newval)


    ## if x is a single output of seqTrack
    if(is.vector(x$ances)){
        newId <- newval # new values
        newAnces <- newval[x$ances] # new values
        ## make output
        res <- x
        res$id <- newId
        res$ances <- newAnces
        attr(res$ances, "levels") <- newlev
    }


    ## if x is an optimize.seqTrack output
    if(is.matrix(x$ances)){
        res <- x
        ori.ncol <- ncol(res$ances)
        res$ances <- matrix(newval[res$ances], ncol=ori.ncol)
        attr(res$ances, "levels") <- newlev
    }

    ## method for haploSim
    if(inherits(x,"haploSim")){
        res <- x
        ances.id <- match(x$ances, labels(x))
    }

    return(res)

} # end get.result.by






#################
## get.consensus
#################
get.consensus <- function(orires, listres, mode=c("majority","best")){
    ## handle mode
    mode <- match.arg(mode)

    res <- orires

    if(mode=="majority"){
        nbDraws <- 0

        ## tables of occurences of ancestors
        temp <- apply(listres$ances, 1, table)

        ## compute compromise
        if(!is.list(temp)){
            newances <- temp
            ances.support <- rep(1,length(temp))
        } else {
            f1 <- function(tab){
                if(length(tab)==0) return(NA)

                res <- names(tab)[tab==max(tab)]
                ## if(length(res)==1) return(res)
                ##             return(NA)
                if(length(res)>1) {
                    nbDraws <- nbDraws+1
                }
                return(res[1])
            }

            newances <- sapply(temp, f1)
            ances.support <- sapply(temp, function(e) max(e, na.rm=TRUE)/sum(e, na.rm=TRUE))
            ances.support[is.na(newances)] <- NA
        }

        ## form the output
        olev <- levels(orires$ances)
        res$ances <- newances
        levels(res$ances) <- olev
        res$support <- ances.support
        res$weight <- rep(1, length(res$date))

        if(is.numeric(listres$ances)){
            res$ances <- as.numeric(res$ances)
        }
        cat("\nThere were\n",nbDraws, "draws.\n")
    } # end majority


    if(mode=="best"){
        toKeep <- which.min(listres$valsim)
        nbDraws <- sum(listres$valsim < (min(listres$valsim) + 1e-10 )) -1
        cat("\nThere were\n",nbDraws, "draws.\n")

        res$ances <- listres$ances[,toKeep]
        res$inf.date <- listres$date[,toKeep]
        res$ances.date <- listres$ances.date[,toKeep]
        res$weight <- rep(-1, length(res$date))
    }

    return(res)
} # end get.consensus







###########################
## get.likelihood.seqTrack
###########################
get.likelihood.seqTrack <-function(x, method=("genetic"), mu0=NULL, seq.length=NULL,...){
    method <- match.arg(method)
    if(method=="genetic"){ # p(nb mutations occur in the time interval)
        if(any(na.omit(res$weight - round(res$weight) > 1e-10))){
            warning("Non-integer weights: number of mutations expected in x$weight.")
        }

        if(is.null(mu0)) stop("mu0 is required.")
        if(is.null(seq.length)) stop("seq.length is required.")

        dates <- as.POSIXct(x$date)
        anc.dates <- as.POSIXct(x$ances.date)
        nb.days <- abs(as.integer(anc.dates-dates))
        nb.mut <- x$weight
        mu <- mu0/365
        mu <- mu*nb.days

        res <- dbinom(nb.mut, size=seq.length, prob=mu)
    } else{
        cat("Method not implemented.")
    }

    return(res)
} # end get.likelihood.seqTrack






###############
## seqTrack.ml
###############
seqTrack.ml <- function(aln, x.dates, optim=c("min","max"), ...){

} # end seqTrack.ml























































#########################
## OLD ADD-HOC VERSION
#########################

## #############
## ## seqTrack
## #############
## seqTrack <- function(x.names, x.dates, D, k=5, lag=3, ...){

##     ## CHECKS ##
##     if(length(x.names) != length(x.dates)){
##         stop("inconsistent length for x.dates")
##     }

##     D <- as.matrix(D)

##     if(length(x.names) != nrow(D)){
##         stop("inconsistent dimension for D")
##     }


##     ## ASSIGNEMENTS ##
##     N <- length(x.names)
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
##         return(STARTPOINTS[which.max(x.dates[STARTPOINTS])])
##     }


##     ## find k closest ancestors: returns a list(id=[id of the ancestors], d=[distances to ancestors])
##     findAncestors <- function(currentPoint){
##         candidates <- id[x.dates < x.dates[currentPoint]]
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
##             if(all(x.dates[temp] < x.dates[oldest])) return(TRUE)
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
