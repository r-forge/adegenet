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
haploPop <- function(n.steps=20, haplo.length=1e6, mu=1e-4, n.snp.ini=10,
                     r.func=function(Nt){max(0, Nt * rnorm(1, mean=1.2, sd=.2))}, gen.time=1,
                     pop.ini.size=function(){1e1}, pop.max.size=function(){1e4}, max.nb.pop=100,
                     p.new.pop=function(){1e-4}, kill.func=function(age){age>1}, quiet=FALSE) {


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

     if(is.numeric(kill.func)){
        kill.func.val <- kill.func[1]
        kill.func <- function(age){age>p.new.pop.val}
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

    evolveOnePop <- function(myPop, myS, myAge){ # myPop: pop to evolve; myS: nb of susceptible in the pop; myAge: vector of ages
        ## kill 'em bastards (= old strains)
        myAge <- myAge + 1
        myPop[kill.func(myAge)] <- NULL

        ## sample and mutate
        sampSize <- round(min( r.func(Nt=length(myPop)), myS)) # number of strains for next step
        if(sampSize<1){ # if no sample
            return(list(pop=myPop, S=myS, age=myAge))
        }
        res <- myPop[sample(1:length(myPop), sampSize, replace=TRUE)] # sample strains
        res <- assignMutations(res, createMutations(sampSize)) # mutate strains
        myAge <- rep(0, pop.ini.size()) # new ages for newborns

        ## possibly create one or more new pop
        if(sum(vecS>0) < max.nb.pop) { # total number of pop. limitation
            nbNewPop <- rbinom(1, sampSize, prob=p.new.pop())
        } else {
            nbNewPop <- 0
        }
        if(nbNewPop>0){
            newPop <- sample(listPop, size=nbNewPop, replace=TRUE)
            listPop <<- c(listPop, newPop)
            vecS <<- c(vecS, replicate(nbNewPop,pop.max.size()) )
            listAges <<- c(listAges, lapply(1:nbNewPop, function(i) rep(0, pop.ini.size())) )
        } # end new pop
        return(list(pop=res, S=myS-sampSize, age=myAge ))
    }


    ## INITIATE SIMULATIONS ##
    vecS <- pop.max.size() # susceptibles
    haplo.ini <- sample(SNP.POOL, n.snp.ini, replace=TRUE)
    listPop <- list()
    listPop[[1]] <- lapply(1:pop.ini.size(), function(i) haplo.ini) # contains only one population of identical clones to start with
    listAges <- list() # will contain vectors of ages of haplotypes (a time of appearance, age=0)
    listAges[[1]] <- rep(0, pop.ini.size())


    ## MAKE SIMULATIONS ##

    ## evolve all populations
    i <- 1L
    if(!quiet){
        cat("\nSimulating populations of haplotypes through time: \n")
    }
    ##while((sum(vecS)>0) & (i<(n.steps+1))){ # evolve all generations
    while(i<(n.steps+1)){ # evolve all generations
        i <- i + 1L # update iterator
        if(!quiet){
            cat(ifelse((i%%10)==0, i, "."))
        }

        ## ## purge non-susceptible pop
        ## listPop <- listPop[vecS>0]
        ## vecS <- vecS[vecS>0]

        ## purge empty populations
        toKeep <- sapply(listPop, length)>0
        listPop <- listPop[toKeep]
        vecS <- vecS[toKeep]
        listAges <- listAges[toKeep]

        ## stop if all pop go extinct
        if(length(listPop)==0L){
            cat("\n All populations went extinct at time",i,"\n")
            return(invisible(NULL))
        }

        ## make populations evolve of one generation
        temp <- 1:length(listPop) # make sure that new pop won't evolve this time
        for(j in temp){
            temp<- evolveOnePop(listPop[[j]], vecS[j], listAges[[j]])
            listPop[[j]] <- temp$pop
            vecS[j] <- temp$S
            listAges[[j]] <- temp$age
        }
    } # end while

    if(!quiet){
        cat("\n... done! \n")
    }

    ## END OF SIMULATIONS ##


    ## CLEAN RESULTS ##
    if(!quiet){
        cat("\n... Cleaning haplotypes (handling reverse mutations)\n")
    }

    ## handle reverse mutations
    cleanRes <- function(vec){
        temp <- table(vec)
        return(sort(as.integer(names(temp)[temp %% 2 != 0])))
    }

    for(i in 1:length(listPop)){
        listPop[[i]] <- lapply(listPop[[i]], cleanRes)
    }

    if(!quiet){
        cat("\n... done! \n")
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
    ## x <- as.array(x)
    n <- length(x)

    f1 <- function(a,b){
        return(sum(!union(a,b) %in% intersect(a,b)))
    }

    ## res <- outer(x, x, FUN=f1)
    res <- matrix(0, ncol=n, nrow=n)
    for(i in 1:(n-1)){
        for(j in (i+1):n){
            res[i,j] <- f1(x[[i]], x[[j]])
        }
    }

    res <- res+t(res)

    return(as.dist(res))
} # end dist.haploPop
