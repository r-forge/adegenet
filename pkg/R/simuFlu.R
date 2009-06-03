###########
## simuFlu
###########
##
## N: number of sequences to simulate
## Tmax: periode of time to simulate
## mean.t.dupli, sd.t.dupli: average time for transmission and its standard deviation (normal dist)
## mean.nb.dupli, sd.nb.dupli: average number of transmissions and its standard deviation (normal dist)
##
simuFlu <- function(N=100, seq.length=100, mu=0.0035,
                    Tmax=20, mean.t.dupli=8, sd.t.dupli=2,
                    mean.nb.dupli=3, sd.nb.dupli=1){

    ## GENERAL VARIABLES ##
    NUCL <- c("a","t","c","g")
    res <- list(seq=list(),date=integer())
    toExpand <- logical()


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
        res <- round(rnorm(1, mean=mean.t.dupli, sd=sd.t.dupli))
        res[res<0] <- 0
        return(res)
    }

    ## when duplication occurs?
    date.dupli <- function(curDate){
        res <- curDate + time.dupli()
        return(res)
    }

    ## how many duplication/transmission occur?
    nb.dupli <- function(){
        res <- round(rnorm(1, mean=mean.nb.dupli, sd=sd.nb.dupli))
        res[res<0] <- 0
        return(res)
    }


    ## MAIN SUB-FUNCTION: EXPANDING FROM ONE SEQUENCE ##
    expand.one.strain <- function(seq, date, idx){
        toExpand[idx] <<- FALSE # this one is no longer to expand
        nbDes <- nb.dupli()
        if(nbDes==0) return(NULL) # stop if no descendant
        newDates <- date.dupli(date) # find dates for descendants
        newDates <- newDates[newDates <= Tmax] # don't store future sequences
        nbDes <- length(newDates)
        if(nbDes==0) return(NULL) # stop if no suitable date
        newSeq <- lapply(1:nbDes, function(i) seq.dupli(seq)) # generate new sequences
        res$seq <<- c(res$seq, newSeq) # append to general output
        res$dates <<- c(res$dates, newDates) # append to general output
        toExpand <<- c(toExpand, rep(TRUE, nbDes))
    }


    ## PERFORM SIMULATIONS ##

    ## initialization
    res$seq[[1]] <- gen.seq()
    res$dates[1] <- 0
    toExpand <- TRUE

    while(any(toExpand)){
        idx <- which.min(toExpand)
        expand.one.strain(seq[[idx]], date[idx], idx)
    }

    return(res)

} # end simuFlu
