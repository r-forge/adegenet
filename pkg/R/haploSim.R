############
## haploSim
############
##
## N: number of sequences to simulate
## mu: mutation rate per nucleotid per year
## Tmax: periode of time to simulate
## mean.gen.time, sd.gen.time: average time for transmission and its standard deviation (normal dist)
## mean.repro, sd.repro: average number of transmissions and its standard deviation (normal dist)
##
haploSim <- function(seq.length=1000, mu=0.0001,
                     Tmax=50, mean.gen.time=5, sd.gen.time=1,
                     mean.repro=2, sd.repro=1, max.nb.haplo=1e3,
                     geo.sim=TRUE, grid.size=10, lambda.xy=0.5, lambdaMat.xy=NULL){

    ## CHECKS ##
    if(!require(ape)) stop("The ape package is required.")


    ## GENERAL VARIABLES ##
    NUCL <- as.DNAbin(c("a","t","c","g"))
    res <- list(seq=as.matrix(as.DNAbin(character(0))), dates=integer(), ances=character())
    toExpand <- logical()
    mu <- mu/365 # mutation rate by day


    ## AUXILIARY FUNCTIONS ##
    ## generate sequence from scratch
    seq.gen <- function(){
        res <- list(sample(NUCL, size=seq.length, replace=TRUE))
        class(res) <- "DNAbin"
        return(res)
    }

    ## create substitutions for defined SNPs
    substi <- function(snp){
        res <- sapply(snp, function(e) sample(setdiff(NUCL,e),1))
        class(res) <- "DNAbin"
        return(res)
    }

    ## duplicate a sequence (including possible mutations)
    seq.dupli <- function(seq){
        toChange <- as.logical(rbinom(n=seq.length, size=1, prob=mu))
        res <- seq
        if(sum(toChange)>0) {
            res[toChange] <- substi(res[toChange])
        }
        return(res)
    }

    ## what is the name of the new sequences?
    seqname.gen <- function(nb.new.seq){
        res <- max(as.integer(rownames(res$seq)), 0) + 1:nb.new.seq
        return(as.character(res))
    }

    ## how many days before duplication occurs ?
    time.dupli <- function(){
        res <- round(rnorm(1, mean=mean.gen.time, sd=sd.gen.time))
        res[res<0] <- 0
        return(res)
    }

    ## when duplication occurs?
    date.dupli <- function(curDate){
        res <- curDate + time.dupli()
        return(res)
    }

    ## how many duplication/transmission occur?
    nb.desc <- function(){
        res <- round(rnorm(1, mean=mean.repro, sd=sd.repro))
        res[res<0] <- 0
        return(res)
    }

      ## where does an haplotype emerges in the first place?
    xy.gen <- function(){
        return(sample(1:grid.size, size=2, replace=TRUE))
    }

    ## where does a transmission occur (destination)?
    if(is.null(lambdaMat.xy)){ # use universal lambda param
        xy.dupli <- function(cur.xy, nbLoc){
            mvt <- rpois(2*nbLoc, lambda.xy) * sample(c(-1,1), size=2*nbLoc, replace=TRUE)
            res <- t(matrix(mvt, nrow=2) + as.vector(cur.xy))
            res[res < 1] <- 1
            res[res > grid.size] <- grid.size
            return(res)
        }
    } else { # use location-dependent lambda param
        xy.dupli <- function(cur.xy){
            lambda.xy <- lambdaMat.xy[cur.xy[1] , cur.xy[2]]
            mvt <- rpois(2, lambda.xy) * sample(c(-1,1), size=2, replace=TRUE)
            res <- cur.xy + mvt
            res[res < 1] <- 1
            res[res > grid.size] <- grid.size
            return(res)
        }
    }


    ## check result size and resize it if needed
    resize.result <- function(){
        curSize <- length(res$dates)
        if(curSize <= max.nb.haplo) return(NULL)
        toKeep <- rep(FALSE, curSize)
        toKeep[sample(1:curSize, size=max.nb.haplo, replace=FALSE)] <- TRUE
        removed.strains <- rownames(res$seq)[!toKeep]
        res$seq <<- res$seq[toKeep,]
        res$dates <<- res$dates[toKeep]
        res$ances <<- res$ances[toKeep]
        toExpand <<- toExpand[toKeep]
        temp <- as.character(res$ances) %in% removed.strains
        if(any(temp)) {
            res$ances[temp] <<- NA
        }

        return(NULL)
    }

    ## check result size and resize it if needed - spatial version
    resize.result.xy <- function(){
        curSize <- length(res$dates)
        if(curSize <= max.nb.haplo) return(NULL)
        toKeep <- rep(FALSE, curSize)
        toKeep[sample(1:curSize, size=max.nb.haplo, replace=FALSE)] <- TRUE
        removed.strains <- rownames(res$seq)[!toKeep]
        res$seq <<- res$seq[toKeep,]
        res$dates <<- res$dates[toKeep]
        res$ances <<- res$ances[toKeep]
        res$xy <<- res$xy[toKeep,,drop=FALSE]
        toExpand <<- toExpand[toKeep]
        temp <- as.character(res$ances) %in% removed.strains
        if(any(temp)) {
            res$ances[temp] <<- NA
        }

        return(NULL)
    }



    ## MAIN SUB-FUNCTION: EXPANDING FROM ONE SEQUENCE - NON SPATIAL ##
    expand.one.strain <- function(seq, date, idx){
        toExpand[idx] <<- FALSE # this one is no longer to expand
        nbDes <- nb.desc()
        if(nbDes==0) return(NULL) # stop if no descendant
        newDates <- sapply(1:nbDes, function(i) date.dupli(date)) # find dates for descendants
        newDates <- newDates[newDates <= Tmax] # don't store future sequences
        nbDes <- length(newDates)
        if(nbDes==0) return(NULL) # stop if no suitable date
        newSeq <- lapply(1:nbDes, function(i) seq.dupli(seq)) # generate new sequences
        class(newSeq) <- "DNAbin" # lists of DNAbin vectors must also have class "DNAbin"
        newSeq <- as.matrix(newSeq) # list DNAbin -> matrix DNAbin with nbDes rows
        rownames(newSeq) <- seqname.gen(nbDes) # find new labels for these new sequences
        res$seq <<- rbind(res$seq, newSeq) # append to general output
        res$dates <<- c(res$dates, newDates) # append to general output
        res$ances <<- c(res$ances, rep(rownames(res$seq)[idx], nbDes)) # append to general output
        toExpand <<- c(toExpand, rep(TRUE, nbDes))
        return(NULL)
    }


    ## 2nd MAIN SUB-FUNCTION: EXPANDING FROM ONE SEQUENCE - SPATIAL ##
    expand.one.strain.xy <- function(seq, date, idx, cur.xy){
        toExpand[idx] <<- FALSE # this one is no longer to expand
        nbDes <- nb.desc()
        if(nbDes==0) return(NULL) # stop if no descendant
        newDates <- sapply(1:nbDes, function(i) date.dupli(date)) # find dates for descendants
        newDates <- newDates[newDates <= Tmax] # don't store future sequences
        nbDes <- length(newDates)
        if(nbDes==0) return(NULL) # stop if no suitable date
        newSeq <- lapply(1:nbDes, function(i) seq.dupli(seq)) # generate new sequences
        class(newSeq) <- "DNAbin" # lists of DNAbin vectors must also have class "DNAbin"
        newSeq <- as.matrix(newSeq) # list DNAbin -> matrix DNAbin with nbDes rows
        rownames(newSeq) <- seqname.gen(nbDes) # find new labels for these new sequences
        res$seq <<- rbind(res$seq, newSeq) # append to general output
        res$dates <<- c(res$dates, newDates) # append to general output
        res$ances <<- c(res$ances, rep(rownames(res$seq)[idx], nbDes)) # append to general output
        res$xy <<- rbind(res$xy, xy.dupli(cur.xy, nbDes))
        toExpand <<- c(toExpand, rep(TRUE, nbDes))
        return(NULL)
    }





    ## PERFORM SIMULATIONS - NON SPATIAL CASE ##
    if(!geo.sim){
        ## initialization
        res$seq <- as.matrix(seq.gen())
        rownames(res$seq) <- "1"
        res$dates[1] <- 0
        res$ances[1] <- NA
        toExpand <- TRUE

        ## simulations: isn't simplicity beautiful?
        while(any(toExpand)){
            idx <- min(which(toExpand))
            expand.one.strain(res$seq[idx,], res$dates[idx], idx)
            resize.result()
        }


        ## SHAPE AND RETURN OUTPUT ##
        res$ances <- as.character(res$ances)
        names(res$dates) <- rownames(res$seq)
        class(res) <- "haploSim"
        return(res)

    } # END NON-SPATIAL SIMULATIONS




    ## PERFORM SIMULATIONS - SPATIAL CASE ##
    if(geo.sim){
        ## some checks
        if(!is.null(lambdaMat.xy)) {
            if(nrow(lambdaMax.xy) != ncol(lambdaMat.xy)) stop("lambdaMat.xy is not a square matrix")
            if(nrow(lambdaMax.xy) != grid.size) stop("dimension of lambdaMat.xy does not match grid size")
        }

        ## initialization
        res$seq <- as.matrix(seq.gen())
        rownames(res$seq) <- "1"
        res$dates[1] <- 0
        res$ances[1] <- NA
        res$xy <- matrix(xy.gen(), nrow=1)
        colnames(res$xy) <- c("x","y")
        toExpand <- TRUE

        ## simulations: isn't simplicity beautiful?
        while(any(toExpand)){
            idx <- min(which(toExpand))
            expand.one.strain.xy(res$seq[idx,], res$dates[idx], idx, res$xy[idx,])
            resize.result.xy()
        }


        ## SHAPE AND RETURN OUTPUT ##
        res$ances <- as.character(res$ances)
        names(res$dates) <- rownames(res$seq)

        class(res) <- "haploSim"
        res$call <- match.call()
        return(res)

    } # end SPATIAL SIMULATIONS


} # end haploSim








##################
## print.haploSim
##################
print.haploSim <- function(x, ...){

    cat("\t\n========================")
    cat("\t\n= simulated haplotypes =")
    cat("\t\n=  (haploSim object)   =")
    cat("\t\n========================\n")

    cat("\nSize :", length(x$ances),"haplotypes")
    cat("\nHaplotype length :", ncol(x$seq),"nucleotids")
    cat("\nProportion of NA ancestors :", signif(mean(is.na(x$ances)),5))
    cat("\nNumber of known ancestors :", sum(!is.na(x$ances)))

    cat("\n\n= Content =")
    for(i in 1:length(x)){
        cat("\n")

        cat(paste("$", names(x)[i], sep=""),"\n")
        if(names(x)[i] %in% c("seq","call")) {
            print(x[[i]])
        } else if(names(x)[i]=="xy"){
            print(head(x[[i]]))
            if(nrow(x[[i]]>6)) cat("    ...\n")
        } else cat(head(x[[i]],6), ifelse(length(x[[i]])>6,"...",""),"\n")
    }


    return(NULL)
} # end print.haploSim






##############
## [.haploSim
##############
"[.haploSim" <- function(x,i,j,drop=FALSE){
    res <- x
    res$seq <- res$seq[i,]
    res$ances <- res$ances[i]
    res$dates <- res$dates[i]
    if(!is.null(res$xy)) res$xy <- res$xy[i,]

    return(res)
}





####################
## na.omit.haploSim
####################
##
## ACTUALLY THIS FUNCTION MAKES NO SENSE FOR NOW
## AS STRAINS WITH NO ANCESTOR MAY BE ANCESTORS OF
## OTHER STRAINS.
##
## na.omit.haploSim <- function(object, ...){
##     res <- object
##     isNA <- is.na(res$ances)
##     res$seq <- res$seq[!isNA,]
##     res$ances <- res$ances[!isNA]
##     res$dates <- res$dates[!isNA]
##     if(!is.null(res$xy)) res$xy <- res$xy[!isNA,]

##     return(res)
## }




##################
## labels.haploSim
##################
labels.haploSim <- function(object, ...){
    return(rownames(object$seq))
}



#######################
## as.POSIXct.haploSim
#######################
as.POSIXct.haploSim <- function(x, tz="", origin=as.POSIXct("2000/01/01"), ...){
    res <- as.POSIXct(x$dates*24*3600, origin=origin)
    return(res)
}




#####################
## seqTrack.haploSim
#####################
seqTrack.haploSim <- function(x, optim=c("min","max"), proxMat=NULL, ...){
    myX <- dist.dna(x$seq, model="raw")
    seq.names <- labels(x)
    seq.dates <- as.POSIXct(x)
    res <- seqTrack(myX, seq.names=seq.names, seq.dates=seq.dates, optim=optim, proxMat=proxMat,...)
    return(res)
}





##############################
## optimize.seqTrack.haploSim
##############################
optimize.seqTrack.haploSim <- function(x, thres=0.2, optim=c("min","max"),
                              prox.mat=NULL, nstep=10, step.size=1e3, rMissDate=.rUnifTimeSeq, ...){

    myX <- dist.dna(x$seq, model="raw")
    seq.names <- labels(x)
    seq.dates <- as.POSIXct(x)
    seq.length <- ncol(x$seq)
    prevCall <- as.list(x$call)
    if(is.null(prevCall$mu)){
        mu0 <- 0.0001
    } else {
        mu0 <- eval(prevCall$mu)
    }

    res <- optimize.seqTrack.default(x=myX, seq.names=seq.names, seq.dates=seq.dates,
                                     thres=thres, optim=optim, prox.mat=prox.mat,
                                     nstep=nstep, step.size=step.size, mu0=mu0,
                                     seq.length=seq.length, rMissDate=rMissDate, ...)
} # end optimize.seqTrack.haploSim





########################
## as.seqTrack.haploSim
########################
as.seqTrack.haploSim <- function(x){
    ## x.ori <- x
    ## x <- na.omit(x)
    toSetToNA <- x$dates==min(x$dates)
    res <- list()
    res$id <- labels(x)
    res <- as.data.frame(res)
    res$ances <- x$ances
    res$ances[toSetToNA] <- NA
    res$weight <- 1 # ??? have to recompute that...
    res$weight[toSetToNA] <- NA
    res$date <- as.POSIXct(x)[labels(x)]
    res$ances.date <- as.POSIXct(x)[x$ances]

    ## set results as indices rather than labels
    res$ances <- match(res$ances, res$id)
    res$id <- 1:length(res$id)

    return(res)
}




################
## plotHaploSim
################
plotHaploSim <- function(x, annot=FALSE, dateRange=NULL, col=NULL, bg="grey", add=FALSE, ...){

    ## SOME CHECKS ##
    if(class(x)!="haploSim") stop("x is not a haploSim object")
    if(is.null(x$xy)) stop("x does not contain xy coordinates; try to simulate date")


    ## ## CONVERSION TO A SEQTRACK-LIKE OBJECT ##
    xy <- x$xy
    res <- as.seqTrack.haploSim(x)

    ##     res <- list()
    ##     res$id <- labels(x)
    ##     res <- as.data.frame(res)
    ##     res$ances <- x$ances
    ##     res$ances[toSetToNA] <- NA
    ##     res$weight <- 1 # ??? have to recompute that...
    ##     res$weight[toSetToNA] <- NA
    ##     res$date <- as.POSIXct(x.ori)[labels(x)]
    ##     res$ances.date <- as.POSIXct(x.ori)[x$ances]
    ##     ## set results as indices rather than labels
    ##     res$ances <- match(res$ances, res$id)
    ##     res$id <- 1:length(res$id)


    ## CALL TO PLOTSEQTRACK ##
    plotSeqTrack(res, xy=xy, annot=annot, dateRange=dateRange,
                        col=col, bg=bg, add=add, showAmbiguous=FALSE, ...)

    return(invisible(res))

} # end plotHaploSim
