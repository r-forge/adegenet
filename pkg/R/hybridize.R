##
## Function hybridize takes two genind in inputs
## and generates hybrids individuals having one parent
## in both objects.
##

hybridize <- function(x1, x2, n, res.type=c("genind","df","STRUCTURE"), file=NULL, quiet=FALSE, sep="/"){
    ## checks
    if(!is.genind(x1)) stop("x1 is not a valid genind object")
    if(!is.genind(x2)) stop("x2 is not a valid genind object")
    n <- as.integer(n)
    res.type <- match.arg(res.type)
    if(!all(x1@loc.names==x2@loc.names)) stop("names of markers in x1 and x2 do not correspond")

    ## used variables
    n1 <- nrow(x1$tab)
    n2 <- nrow(x2$tab)
    k <- length(x1$loc.names)
 
    #### get frequencies for each locus
    x1 <- genind2genpop(x1,pop=factor(rep(1,n1)),quiet=TRUE)
    freq1 <- makefreq(x1,quiet=TRUE)$tab
    freq1 <- split(freq1, x1@loc.fac)

    x2 <- genind2genpop(x2,pop=factor(rep(1,n2)),quiet=TRUE)
    freq2 <- makefreq(x2,quiet=TRUE)$tab
    freq2 <- split(freq2, x2@loc.fac)

    #### sampling of gametes
    ## kX1 / kX2 are lists of tables of sampled gametes
    kX1 <- lapply(freq1, function(v) t(rmultinom(n,1,v)))
    names(kX1) <- x1$loc.names
    for(i in 1:k) { colnames(kX1[[i]]) <- x1$all.names[[i]]}
    kX2 <- lapply(freq2, function(v) t(rmultinom(n,1,v)))
    names(kX2) <- x2$loc.names
    for(i in 1:k) { colnames(kX2[[i]]) <- x2$all.names[[i]]}
  
    ## tab1 / tab2 are cbinded tables 
    tab1 <- cbind.data.frame(kX1)
    ## gam 1/2 are genind containing gametes
    ## gam 1
    gam1 <- genind(tab1)
    gam1@loc.names <- x1@loc.names
    gam1@loc.fac <- x1@loc.fac
    gam1@all.names <- x1@all.names
    gam1@loc.nall <- x1@loc.nall
    gam1 <- genind2df(gam1,sep="/")
    gam1 <- as.matrix(gam1)
    
    ## gam 2
    tab2 <- cbind.data.frame(kX2)
    ## gam 1/2 are genind containing gametes
    gam2 <- genind(tab2)
    gam2@loc.names <- x1@loc.names
    gam2@loc.fac <- x1@loc.fac
    gam2@all.names <- x1@all.names
    gam2@loc.nall <- x1@loc.nall
    gam2 <- genind2df(gam2,sep="/")
    gam2 <- as.matrix(gam2)

    #### construction of zygotes
    gam1 <-  gsub("/.*$","",gam1)
    gam2 <-  gsub("/.*$","",gam2)
 
    ## res.type=="STRUCTURE"
    if(res.type=="STRUCTURE"){
        res <- paste(gam1,gam2,sep="")
        res <- as.data.frame(matrix(res,ncol=k))
        names(res) <- x1@loc.names
        row.names(res) <- .genlab("",n)
        if(is.null(file)) {
            file <- gsub("[[:space:]]","-",date())
            file <- paste("hybrid",file,sep="_")
            file <- paste(file,"str",sep=".")
        }
        write.table(res, file=file,row.names = TRUE, col.names = TRUE)
        if(!quiet) cat("\nWrote results to file", file, "\n")

        return(invisible())
    }

    ## res.type=="df"
    if(res.type=="df"){
        res <- paste(gam1,gam2,sep=sep)
        res <- as.data.frame(matrix(res,ncol=k))
        names(res) <- x1@loc.names
        row.names(res) <- .genlab("",n)

        return(res)
    }

    ## res.type=="genind"
    if(res.type=="genind"){
        res <- paste(gam1,gam2,sep="")
        res <- as.data.frame(matrix(res,ncol=k))
        names(res) <- x1@loc.names
        row.names(res) <- .genlab("",n)
        res <- genind(res,prevcall=match.call())
        
        return(res)
    }
    
} # end hybridize
