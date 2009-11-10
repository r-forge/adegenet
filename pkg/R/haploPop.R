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
haploPop <- function(n.steps=10, haplo.length=1e6, mu=0.0001, gen.time=1,
                     n.snp.ini=10,
                     Rfunc=function(Nt){max(0, Nt * rnorm(1, mean=1.2, sd=.2))},
                     pop.ini.size=function(){1e1}, pop.max.size=function(){1e4}, p.new.pop=function(){1e-4} ) {


    ## GLOBAL VARIABLES ##
    mu <- gen.time * (mu/365)
    SNP.POOL <- 1:haplo.length


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
        nbNewPop <- rbinom(1, sampSize, prob=p.new.pop())
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
    while(sum(vecS)>0 & i<(n.steps+1)){ # evolve all generations
        i <- i + 1L # update iterator

        ## purge non-susceptible pop
        listPop <- listPop[vecS>0]
        vecS <- vecS[vecS>0]

        ## evolve populations of one generation
        for(j in which(vecS>0)){
            temp<- evolveOnePop(listPop[[j]], vecS[j])
            listPop[[j]] <- temp$pop
            vecS[j] <- temp$S
        }
    } # end while

    ## RETURN RESULTS ##
    res <- listPop
    class(res) <- "haploPop"
    res$call <- match.call()
    return(res)

} # end haploPop








## ##################
## ## print.haploPop
## ##################
## print.haploPop <- function(x, ...){

##     cat("\t\n========================")
##     cat("\t\n= simulated haplotypes =")
##     cat("\t\n=  (haploPop object)   =")
##     cat("\t\n========================\n")

##     cat("\nSize :", length(x$ances),"haplotypes")
##     cat("\nHaplotype length :", ncol(x$seq),"nucleotids")
##     cat("\nProportion of NA ancestors :", signif(mean(is.na(x$ances)),5))
##     cat("\nNumber of known ancestors :", sum(!is.na(x$ances)))
##     nbAncInSamp <- sum(x$ances %in% labels(x))
##     cat("\nNumber of ancestors within the sample :", nbAncInSamp)
##     cat("\nDate range :", min(x$dates,na.rm=TRUE),"-",max(x$dates,na.rm=TRUE))
##     ##nUniqSeq <- length(unique(apply(as.character(x$seq),1,paste,collapse="")))
##     ##cat("\nNumber of unique haplotypes :", nUniqSeq)

##     cat("\n\n= Content =")
##     for(i in 1:length(x)){
##         cat("\n")

##         cat(paste("$", names(x)[i], sep=""),"\n")
##         if(names(x)[i] %in% c("seq","call")) {
##             print(x[[i]])
##         } else if(names(x)[i]=="xy"){
##             print(head(x[[i]]))
##             if(nrow(x[[i]]>6)) cat("    ...\n")
##         } else cat(head(x[[i]],6), ifelse(length(x[[i]])>6,"...",""),"\n")
##     }


##     return(NULL)
## } # end print.haploPop
