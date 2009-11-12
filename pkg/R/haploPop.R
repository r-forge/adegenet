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
haploPop <- function(n.steps=20, haplo.length=1e6, mu=1e-5, n.snp.ini=1,
                     birth.func=function(){ sample(0:3, 1)},
                     ini.pop.size=function(){1e1}, max.pop.size=function(){1e4}, max.nb.pop=100,
                     p.new.pop=function(){1e-4}, kill.func=function(age){age>1},
                     quiet=FALSE, clean.haplo=FALSE) {


    ## SOME CHECKS
    if(is.numeric(ini.pop.size)){
        ini.pop.size.val <- ini.pop.size
        ini.pop.size <- function(){ini.pop.size.val}
    }

    if(is.numeric(max.pop.size)){
        max.pop.size.val <- max.pop.size
        max.pop.size <- function(){max.pop.size.val}
    }

     if(is.numeric(p.new.pop)){
        p.new.pop.val <- p.new.pop
        p.new.pop <- function(){p.new.pop.val}
    }

    if(is.numeric(birth.func)){
        birth.func.val <- birth.func[1]
        birth.func <- function(){birth.func.val}
    }

    if(is.numeric(kill.func)){
        kill.func.val <- kill.func[1]
        kill.func <- function(age){age>kill.func.val}
    }


    ## GLOBAL VARIABLES ##
    SNP.POOL <- 1:haplo.length
    vecS <- 1 # will be redefined later, but needed for evolveOnePop definition

    ## AUXILIARY FUNCTIONS ##
    createMutations <- function(N){ # L:genome length; N: pop size
        nb.mutations <- sum(rbinom(N, size=haplo.length, prob=mu))
        return( sample(SNP.POOL, size=nb.mutations, replace=TRUE) )
    }

    assignMutations <- function(myPop, mutations){ # mypop: list of `haplotypes'; mutations: vector of SNPs
        if(length(mutations)==0 | length(myPop)==0) return(myPop)
        id <- sample(1:length(myPop), size=length(mutations), replace=TRUE)
        mutations <- split(mutations, id)
        myPop[as.integer(names(mutations))] <- mapply(c, myPop[as.integer(names(mutations))], mutations, SIMPLIFY=FALSE)
        return(myPop)
    }

    evolveOnePop <- function(myPop, myS, myAge){ # myPop: pop to evolve; myS: nb of susceptible in the pop; myAge: vector of ages
        ## kill 'em bastards (= old strains)
        myAge <- myAge + 1
        toKill <- kill.func(myAge)
        myPop[toKill] <- NULL
        myAge <- myAge[!toKill]

        ## generate new strains for new generation
        sampSize <- round(min( length(myPop)*birth.func(), myS)) # number of strains for next step
        if(sampSize<1){ # if no sample
            return(list(pop=myPop, S=myS, age=myAge))
        }
        newGen <- myPop[sample(1:length(myPop), sampSize, replace=TRUE)] # sample strains for new generations
        newGen <- assignMutations(newGen, createMutations(sampSize)) # mutate strains
        newAge <- rep(0, sampSize) # new ages for newborns

        ## merge old and new generation
        myPop <- c(myPop,newGen)
        myAge <- c(myAge, newAge)

        ## possibly create one or more new pop
        if((length(listPop) < max.nb.pop) & (p.new.pop()>0)) { # total number of pop. limitation
            nbNewPop <- rbinom(1, length(myPop), prob=p.new.pop())
        } else {
            nbNewPop <- 0
        }
        if(nbNewPop>0){
            ## newPop <- sample(listPop, size=nbNewPop, replace=TRUE) # wrong
            newPop <- lapply(sample(myPop, size=nbNewPop, replace=TRUE), as.list)
            listPop <<- c(listPop, newPop)
            vecS <<- c(vecS, replicate(nbNewPop, max.pop.size()) )
            listAges <<- c(listAges, replicate(nbNewPop, 0, simplify=FALSE) )
        } # end new pop
        return(list(pop=myPop, S=myS-sampSize, age=myAge))
    }


    ## INITIATE SIMULATIONS ##
    vecS <- max.pop.size() -  n.snp.ini # susceptibles
    haplo.ini <- sample(SNP.POOL, n.snp.ini, replace=TRUE)
    listPop <- list()
    listPop[[1]] <- lapply(1:ini.pop.size(), function(i) haplo.ini) # contains only one population of identical clones to start with
    listAges <- list() # will contain vectors of ages of haplotypes (a time of appearance, age=0)
    listAges[[1]] <- rep(0, ini.pop.size())


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
            catStep <- max(round(n.steps/200), 10)
            cat(ifelse((i %% catStep)==0, paste(" ...", i), ""))
        }


        ## make populations evolve of one generation
        idx <- which(vecS>0) # make sure that new pop won't evolve this time
        if(length(idx)>0){
            for(j in idx){
                temp <- evolveOnePop(listPop[[j]], vecS[j], listAges[[j]])
                listPop[[j]] <- temp$pop
                vecS[j] <- temp$S
                listAges[[j]] <- temp$age
            }
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

        ## FOR DEBUGGING
        ## cat("\n=== ",i," ===")
        ## cat("\nlistPop")
        ## print(listPop)
        ## cat("\nvecS")
        ## print(vecS)
        ## cat("\nlistAges")
        ## print(listAges)
        ## END DEBUGGING
    } # end while

    if(!quiet){
        cat("\n... done! \n")
    }

    ## END OF SIMULATIONS ##


    ## CLEAN RESULTS ##
    ## handle reverse mutations
    if(clean.haplo){
        if(!quiet){
            cat("\n... Cleaning haplotypes (handling reverse mutations)\n")
        }

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
sample.haploPop <- function(x, n, n.pop=NULL){
    x$call <- NULL
    if(!is.null(n.pop)){ # pre-treatment: reduce to n.pop populations with same size
        toKeep <- sapply(x,length)>n
        x <- x[toKeep] # keep only pop large enough
        x <- sample(x, n.pop, replace=FALSE) # keep n.pop populations
        x <- lapply(x, sample, n, replace=FALSE) # make them the same size
    }

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
