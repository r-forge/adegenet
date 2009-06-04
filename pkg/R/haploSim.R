############
## haploSim
############
##
## N: number of sequences to simulate
## mu: mutation rate per nucleotid per year
## Tmax: periode of time to simulatet
## mean.gen.time, sd.gen.time: average time for transmission and its standard deviation (normal dist)
## mean.repro, sd.repro: average number of transmissions and its standard deviation (normal dist)
##
haploSim <- function(seq.length=1000, mu=0.0001,
                    Tmax=30, mean.gen.time=5, sd.gen.time=1,
                    mean.repro=2, sd.repro=1,
                    max.nb.haplo=1e4){

    ## GENERAL VARIABLES ##
    NUCL <- c("a","t","c","g")
    res <- list(seq=list(), dates=integer(), ances=integer())
    toExpand <- logical()
    mu <- mu/365 # mutation rate by day

    ## AUXILIARY FUNCTIONS ##
    ## generate sequence from scratch
    gen.seq <- function(){
        return(sample(NUCL, size=seq.length, replace=TRUE))
    }

    ## create substitutions for defined SNPs
    substi <- function(snp){
        res <- sapply(snp, function(e) sample(setdiff(NUCL,e),1))
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

    ## check result size and resize it if needed
    resize.result <- function(){
        curSize <- length(res$date)
        if(curSize < max.nb.haplo) return(NULL)
        toKeep <- sample(1:curSize, size=max.nb.haplo, replace=FALSE)
        removed.strains <- res$seq[!toKeep]
        res$seq <<- res$res[toKeep]
        res$date <<- res$date[toKeep]
        res$ances <<- res$ances[toKeep]
        res$ances[res$ances %in% removed.strains] <- NA

        return(NULL)
    }


    ## MAIN SUB-FUNCTION: EXPANDING FROM ONE SEQUENCE ##
    expand.one.strain <- function(seq, date, idx){
        toExpand[idx] <<- FALSE # this one is no longer to expand
        nbDes <- nb.desc()
        if(nbDes==0) return(NULL) # stop if no descendant
        newDates <- sapply(1:nbDes, function(i) date.dupli(date)) # find dates for descendants
        newDates <- newDates[newDates <= Tmax] # don't store future sequences
        nbDes <- length(newDates)
        if(nbDes==0) return(NULL) # stop if no suitable date
        newSeq <- lapply(1:nbDes, function(i) seq.dupli(seq)) # generate new sequences
        res$seq <<- c(res$seq, newSeq) # append to general output
        res$dates <<- c(res$dates, newDates) # append to general output
        res$ances <<- c(res$ances, rep(idx, nbDes)) # append to general output
        toExpand <<- c(toExpand, rep(TRUE, nbDes))
        return(NULL)
    }


    ## PERFORM SIMULATIONS ##

    ## initialization
    res$seq[[1]] <- gen.seq()
    res$dates[1] <- 0
    res$ances[1] <- NA
    toExpand <- TRUE

    ## simulations: isn't simplicity beautiful?
    while(any(toExpand)){
        idx <- min(which(toExpand))
        expand.one.strain(res$seq[[idx]], res$dates[idx], idx)
        resize.result()
    }


    ## RETURN OUTPUT ##
    return(res)

} # end haploSim
