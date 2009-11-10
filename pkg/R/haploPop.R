############
## haploPop
############
##
## Simulate only SNPs, allow reverse mutations.
##
## - haplo.length: length of simulated haplotypes
## - mu: substitution rate / nucleotide / year
## - n.steps: number of generations to simulate
##
haploPop <- function(n.steps=20, haplo.length=1e6, mu=1e-4, gen.time=1,
                     n.snp.ini=10,
                     Rfunc=function(Nt){max(0, Nt * rnorm(1, mean=1.2, sd=.2))},
                     pop.ini.size=function(){1e1}, pop.max.size=function(){1e4}, p.new.pop=function(){1e-4},
                     max.nb.pop=100) {


    ## SOME CHECKS
    if(is.numeric(pop.ini.size)){
        pop.ini.size.val <- pop.ini.size
        pop.ini.size <- function(){pop.ini.size.val}
    }

    if(is.numeric(pop.max.size)){
        pop.max.size.val <- pop.max.size
        pop.max.size <- function(){pop.max.size.val}
    }

     if(is.numeric(p.new.pop)){
        p.new.pop.val <- p.new.pop
        p.new.pop <- function(){p.new.pop.val}
    }



    ## GLOBAL VARIABLES ##
    mu <- gen.time * (mu/365)
    SNP.POOL <- 1:haplo.length
    vecS <- 1 # will be redefined later, but needed for evolveOnePop definition

    ## AUXILIARY FUNCTIONS ##
    createMutations <- function(N){ # L:genome length; N: pop size
        nb.mutations <- sum(rbinom(N, size=haplo.length, prob=mu))
        return( sample(SNP.POOL, size=nb.mutations, replace=TRUE) )
    }

    assignMutations <- function(myPop, mutations){ # mypop: list of `haplotypes'; mutations: vector of SNPs
        if(length(mutations)==0) return(myPop)
        id <- sample(1:length(myPop), size=length(mutations), replace=TRUE)
        mutations <- split(mutations, id)
        myPop[as.integer(names(mutations))] <- mapply(c, myPop[as.integer(names(mutations))], mutations, SIMPLIFY=FALSE)
        return(myPop)
    }

    evolveOnePop <- function(myPop, myS){ # myPop: pop to evolve; myS: nb of susceptible in the pop
        ## sample and mutate
        sampSize <- round(min( Rfunc(Nt=length(myPop)), myS)) # number of strains for next step
        res <- myPop[sample(1:length(myPop), sampSize, replace=TRUE)] # sample strains
        res <- assignMutations(res, createMutations(sampSize)) # mutate strains

        ## possibly create a new pop
        if(sum(vecS>0) < max.nb.pop) { # total number of pop. limitation
            nbNewPop <- rbinom(1, sampSize, prob=p.new.pop())
        } else {
            nbNewPop <- 0
        }
        if(nbNewPop>0){
            newPop <- sample(listPop, size=nbNewPop, replace=TRUE)
            listPop <<- c(listPop, newPop)
            vecS <<- c(vecS, replicate(nbNewPop,pop.max.size()) )
        }
        return(list(pop=res, S=myS-sampSize ))
    }


    ## INITIATE SIMULATIONS ##
    vecS <- pop.max.size() # susceptibles
    haplo.ini <- sample(SNP.POOL, n.snp.ini, replace=TRUE)
    listPop <- list()
    listPop[[1]] <- lapply(1:pop.ini.size(), function(i) haplo.ini) # contains only one population of identical clones to start with


    ## MAKE SIMULATIONS ##

    ## evolve all populations
    i <- 1L
    while((sum(vecS)>0) & (i<(n.steps+1))){ # evolve all generations
        i <- i + 1L # update iterator

        ## purge non-susceptible pop
        listPop <- listPop[vecS>0]
        vecS <- vecS[vecS>0]

        ## evolve populations of one generation
        temp <- which(vecS>0)
        for(j in temp){
            temp<- evolveOnePop(listPop[[j]], vecS[j])
            listPop[[j]] <- temp$pop
            vecS[j] <- temp$S
        }
    } # end while


    ## CLEAN RESULTS ##
    ## handle reverse mutations
    cleanRes <- function(vec){
        temp <- table(vec)
        return(sort(as.integer(names(temp)[temp %% 2 != 0])))
    }

    for(i in 1:length(listPop)){
        listPop[[i]] <- lapply(listPop[[i]], cleanRes)
    }


    ## RETURN RESULTS ##
    res <- listPop
    class(res) <- "haploPop"
    res$call <- match.call()
    return(res)

} # end haploPop








##################
## summary.haploPop
##################
summary.haploPop <- function(object, ...){
    x <- object
    myCall <- x$call
    x$call <- NULL
    res <- list()

    ## cat("\t\n=======================================")
    ## cat("\t\n= simulated populations of haplotypes =")
    ## cat("\t\n=          (haploPop object)          =")
    ## cat("\t\n=======================================\n")

    cat("\nNumber of populations :", length(x))

    cat("\nPopulation sizes :\n")
    temp <- sapply(x,length)
    names(temp) <- 1:length(temp)
    print(temp)
    res$pop.size <- temp

    cat("\nNumber of SNPs per population :\n")
    temp <- sapply(x,function(e) length(unique(unlist(e))))
    names(temp) <- 1:length(temp)
    print(temp)
    res$n.snp <- temp

    return(invisible(res))
} # end print.haploPop






##################
## sample.haploPop
##################
sample.haploPop <- function(x, n){
    x <- unlist(x, recursive=FALSE)
    res <- list()
    res[[1]] <- sample(x, n)
    class(res) <- "haploPop"
    return(res)
} # end sample.haploPop






###############
## dist.haploPop
###############
dist.haploPop <- function(x){
    if(!inherits(x, "haploPop")) stop("x is not a haploPop object")

    x <- unlist(x, recursive=FALSE)
    f1 <- function(a,b){
        return(sum(!union(a,b) %in% intersect(a,b)))
    }

    res <- outer(x, x, FUN=f1)
    return(as.dist(res))
} # end dist.haploPop
