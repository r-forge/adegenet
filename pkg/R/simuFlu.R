simuFlu <- function(N=100, seq.length=100, mu=0.0035){
    NUCL <- c("a","t","c","g")
    res <- list()

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

}
