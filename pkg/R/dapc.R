#############
## find.clusters
#############
find.clusters <- function (x, ...) UseMethod("find.clusters")

######################
## find.clusters.data.frame
######################
find.clusters.data.frame <- function(x, n.pca=NULL, n.clust=NULL, stat=c("BIC", "AIC", "WSS"), choose.n.clust=TRUE, criterion=c("min","diff"),
                                     max.n.clust=round(nrow(x)/10), n.iter=1e6, n.start=100, center=TRUE, scale=TRUE){

    ## CHECKS ##
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")
    if(!require(MASS, quiet=TRUE)) stop("MASS library is required.")
    if(!require(stats)) stop("package stats is required")

    stat <- match.arg(stat)


    ## SOME GENERAL VARIABLES ##
    N <- nrow(x)
    min.n.clust <- 2
    max.n.clust <- max(max.n.clust, 2)

    ## PERFORM PCA ##
    maxRank <- min(dim(x))

    pcaX <- dudi.pca(x, center = center, scale = scale, scannf = FALSE, nf=maxRank)

    ## select the number of retained PC for PCA
    if(is.null(n.pca)){
        cumVar <- 100 * cumsum(pcaX$eig)/sum(pcaX$eig)
        plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA")
        cat("Choose the number PCs to retain (>=1): ")
        n.pca <- as.integer(readLines(n = 1))
    }

     ## keep relevant PCs - stored in XU
    X.rank <- length(pcaX$eig)
    n.pca <- min(X.rank, n.pca)
    if(n.pca >= N) warning("number of retained PCs of PCA is greater than N")
    if(n.pca > N/3) warning("number of retained PCs of PCA may be too large (> N /3)")

    XU <- pcaX$li[, 1:n.pca, drop=FALSE] # principal components

    ## PERFORM K-MEANS
    if(is.null(n.clust)){
        nbClust <- min.n.clust:max.n.clust
        WSS <- numeric(0)

        for(i in 1:length(nbClust)){
            temp <- kmeans(XU, centers=nbClust[i], iter.max=min(n.iter, 100), nstart=min(n.start, 1e3))
            WSS[i] <- sum(temp$withinss)
        }


        ## DETERMINE THE NUMBER OF GROUPS
        ##TSS <- sum(pcaX$eig) * N
        ##betweenVar <- (1 - ((stat/(N-nbClust-1))/(TSS/(N-1)) )) *100
        ##WSS.ori <- sum(apply(XU, 2, function(v) sum((v-mean(v))^2) ))
        ##reducWSS <- -diff(c(WSS.ori, stat))
        ##reducWSS <- reducWSS/max(reducWSS)

        if(stat=="AIC"){
            WSS.ori <- sum(apply(XU, 2, function(v) sum((v-mean(v))^2) ))
            k <- nbClust
            myStat <- N*log(c(WSS.ori,WSS)/N) + 2*c(1,nbClust)
            myLab <- "AIC"
            myTitle <- "Value of AIC \nversus number of clusters"

        }
        if(stat=="BIC"){
            WSS.ori <- sum(apply(XU, 2, function(v) sum((v-mean(v))^2) ))
            k <- nbClust
            myStat <- N*log(c(WSS.ori,WSS)/N) + log(N) *c(1,nbClust)
            myLab <- "BIC"
            myTitle <- "Value of BIC \nversus number of clusters"
        }
        if(stat=="WSS"){
            WSS.ori <- sum(apply(XU, 2, function(v) sum((v-mean(v))^2) ))
            myStat <- c(WSS.ori, WSS)
            ##            reducWSS <- -diff(c(WSS.ori, stat))
            ##            myStat <- reducWSS/max(reducWSS)
            myLab <- "Within sum of squares"
            myTitle <- "Value of within SS\nversus number of clusters"
        }

        if(choose.n.clust){
            plot(c(1,nbClust), myStat, xlab="Number of clusters", ylab=myLab, main=myTitle, type="b", col="blue")
            abline(h=0, lty=2, col="red")
            cat("Choose the number of clusters (>=2: ")
            n.clust <- as.integer(readLines(n = 1))
        } else {
            if(criterion=="min") {
                n.clust <- which.min(myStat)
            }
            if(criterion=="diff") {
                temp <- diff(myStat)
                n.clust <- which.max( which( (temp-min(temp))<max(temp)/1e4))
            }
        }
    }

    ## get final groups
    best <-  kmeans(XU, centers=n.clust, iter.max=n.iter, nstart=n.start)


    ## MAKE RESULT ##
    if(is.null(n.clust)){
        names(myStat) <- paste("K",c(1,nbClust), sep="=")
    } else {
        myStat <- sum(best$withinss)
    }

    res <- list(stat=myStat, grp=factor(best$cluster), size=best$size)

    return(res)
} # end find.clusters.data.frame






###################
## find.clusters.genind
###################
find.clusters.genind <- function(x, n.pca=NULL, n.clust=NULL, stat=c("BIC", "AIC", "WSS"), choose.n.clust=TRUE,
                          max.n.clust=round(nrow(x@tab)/10), n.iter=1e6, n.start=100,
                          scale=TRUE, scale.method=c("sigma", "binom"), truenames=TRUE){

    ## CHECKS ##
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")
    if(!require(MASS, quiet=TRUE)) stop("MASS library is required.")
    if(!require(stats)) stop("package stats is required")
    if(!is.genind(x)) stop("x must be a genind object.")
    stat <- match.arg(stat)


    ## SOME GENERAL VARIABLES ##
    N <- nrow(x@tab)
    min.n.clust <- 2

    ## PERFORM PCA ##
    maxRank <- min(dim(x@tab))

    X <- scaleGen(x, center = TRUE, scale = scale, method = scale.method,
                  missing = "mean", truenames = truenames)

    ## CALL DATA.FRAME METHOD
    res <- find.clusters(X, n.pca=n.pca, n.clust=n.clust, stat=stat, max.n.clust=max.n.clust, n.iter=n.iter, n.start=n.start,
                         choose.n.clust=choose.n.clust, center=FALSE, scale=FALSE)
    return(res)
} # end find.clusters.genind





###################
## find.clusters.matrix
###################
find.clusters.matrix <- function(x, ...){
    return(find.clusters(as.data.frame(x), ...))
}






########
## dapc
########
dapc <- function (x, ...) UseMethod("dapc")

#################
## dapc.data.frame
#################
dapc.data.frame <- function(x, grp, n.pca=NULL, n.da=NULL,
                            center=TRUE, scale=TRUE, var.contrib=FALSE){

    ## FIRST CHECKS
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")
    if(!require(MASS, quiet=TRUE)) stop("MASS library is required.")
    grp <- as.factor(grp)
    if(length(grp) != nrow(x)) stop("Inconsistent length for grp")


    ## SOME GENERAL VARIABLES
    N <- nrow(x)

    ## PERFORM PCA ##
    maxRank <- min(dim(x))

    pcaX <- dudi.pca(x, center = center, scale = scale, scannf = FALSE, nf=maxRank)

    ## select the number of retained PC for PCA
    if(is.null(n.pca)){
        cumVar <- 100 * cumsum(pcaX$eig)/sum(pcaX$eig)
        plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA")
        cat("Choose the number PCs to retain (>=1): ")
        n.pca <- as.integer(readLines(n = 1))
    }

    ## keep relevant PCs - stored in XU
    X.rank <- length(pcaX$eig)
    n.pca <- min(X.rank, n.pca)
    if(n.pca >= N) stop("number of retained PCs of PCA is greater than N")
    if(n.pca > N/3) warning("number of retained PCs of PCA may be too large (> N /3)")

    U <- pcaX$c1[, 1:n.pca, drop=FALSE] # principal axes
    XU <- pcaX$li[, 1:n.pca, drop=FALSE] # principal components
    XU.lambda <- sum(pcaX$eig[1:n.pca])/sum(pcaX$eig) # sum of retained eigenvalues
    names(U) <- paste("PCA-pa", 1:ncol(U), sep=".")
    names(XU) <- paste("PCA-pc", 1:ncol(XU), sep=".")


     ## PERFORM DA ##
    ldaX <- lda(XU, grp)
    if(is.null(n.da)){
        barplot(ldaX$svd^2, xlab="Linear Discriminants", ylab="F-statistic", main="Discriminant analysis eigenvalues", col=heat.colors(length(levels(grp))) )
        cat("Choose the number discriminant functions to retain (>=1): ")
        n.da <- as.integer(readLines(n = 1))
    }

    n.da <- min(n.da, length(levels(grp))-1) # can't be more than K-1 disc. func.
    predX <- predict(ldaX, dimen=n.da)


    ## BUILD RESULT
    res <- list()
    res$n.pca <- n.pca
    res$n.da <- n.da
    res$tab <- XU
    res$grp <- grp
    res$var <- XU.lambda
    res$eig <- ldaX$svd^2
    res$disc.func <- ldaX$scaling[, 1:n.da, drop=FALSE]
    res$ind.coord <-predX$x
    res$grp.coord <- apply(res$ind.coord, 2, tapply, grp, mean)
    res$prior <- ldaX$prior
    res$posterior <- predX$posterior
    res$assign <- predX$class
    res$call <- match.call()

    ## optional: get loadings of alleles
    if(var.contrib){
        res$var.contr <- as.matrix(U) %*% as.matrix(ldaX$scaling)
        res$var.contr <- t(apply(res$var.contr, 1, function(e) e*e / sum(e*e)))
    }

    class(res) <- "dapc"
    return(res)
} # end dapc.data.frame





#############
## dapc.matrix
#############
dapc.matrix <- function(x, ...){
    return(dapc(as.data.frame(x), ...))
}




#############
## dapc.genind
#############
dapc.genind <- function(x, pop=NULL, n.pca=NULL, n.da=NULL,
                 scale=TRUE, scale.method=c("sigma", "binom"), truenames=TRUE, all.contrib=FALSE){

    ## FIRST CHECKS
    if(!is.genind(x)) stop("x must be a genind object.")

    if(is.null(pop)) {
        pop.fac <- pop(x)
    } else {
        pop.fac <- pop
    }

    if(is.null(pop.fac)) stop("x does not include pre-defined populations, and `pop' is not provided")
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")
    if(!require(MASS, quiet=TRUE)) stop("MASS library is required.")


    ## SOME GENERAL VARIABLES
    N <- nrow(x@tab)

    ## PERFORM PCA ##
    maxRank <- min(dim(x@tab))

    X <- scaleGen(x, center = TRUE, scale = scale, method = scale.method,
                  missing = "mean", truenames = truenames)

    ## CALL DATA.FRAME METHOD ##
   res <- dapc(X, grp=pop.fac, n.pca=n.pca, n.da=n.da,
                            center=FALSE, scale=FALSE, var.contrib=all.contrib)

    res$call <- match.call()

    return(res)
} # end dapc.genind






######################
# Function print.dapc
######################
print.dapc <- function(x, ...){
    cat("\t#########################################\n")
    cat("\t# Discriminant Analysis of Principal Components #\n")
    cat("\t#########################################\n")
    cat("class: ")
    cat(class(x))
    cat("\n$call: ")
    print(x$call)
    cat("\n$n.pca:", x$n.pca, "first PCs of PCA used")
    cat("\n$n.da:", x$n.da, "discriminant functions saved")
    cat("\n$varn (proportion of conserved variance):", round(x$var,3))
    cat("\n\n$eig (eigenvalues): ")
    l0 <- sum(x$eig >= 0)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5)
        cat(" ...\n\n")

    ## vectors
    sumry <- array("", c(4, 3), list(1:4, c("vector", "length", "content")))
    sumry[1, ] <- c('$eig', length(x$eig),  'eigenvalues')
    sumry[2, ] <- c('$grp', length(x$grp), 'prior group assignment')
    sumry[3, ] <- c('$prior', length(x$prior), 'prior group probabilities')
    sumry[4, ] <- c('$assign', length(x$assign), 'posterior group assignment')
    class(sumry) <- "table"
    print(sumry)

    ## data.frames
    cat("\n")
    sumry <- array("", c(5, 4), list(1:5, c("data.frame", "nrow", "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "retained PCs of PCA")
    sumry[2, ] <- c("$disc.func", nrow(x$disc.func), ncol(x$disc.func), "discriminant functions")
    sumry[3, ] <- c("$ind.coord", nrow(x$ind.coord), ncol(x$ind.coord), "coordinates of individuals")
    sumry[4, ] <- c("$grp.coord", nrow(x$grp.coord), ncol(x$grp.coord), "coordinates of groups")
    sumry[5, ] <- c("$posterior", nrow(x$posterior), ncol(x$posterior), "posterior membership probabilities")
    class(sumry) <- "table"
    print(sumry)

    cat("\nother elements: ")
    if (length(names(x)) > 13)
        cat(names(x)[14:(length(names(x)))], "\n")
    else cat("NULL\n")
}






##############
## summary.dapc
##############
summary.dapc <- function(object, ...){
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")

    x <- object
    res <- list()

    ## number of dimensions
    res$n.dim <- ncol(x$disc.func)
    res$n.pop <- length(levels(x$grp))

    ## assignment success
    temp <- as.character(x$grp)==as.character(x$assign)
    res$assign.prop <- mean(temp)
    res$assign.per.pop <- tapply(temp, x$grp, mean)

    ## group sizes
    res$prior.grp.size <- table(x$grp)
    res$post.grp.size <- table(x$assign)

    return(res)
} # end summary.dapc






##############
## scatter.dapc
##############
scatter.dapc <- function(x, xax=1, yax=2, col=rainbow(length(levels(x$grp))), posi="bottomleft", bg="grey", ratio=0.3, csub=1.2, ...){
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")
    axes <- c(xax,yax)
    par(bg=bg)
    s.class(x$ind.coord[,axes], fac=x$grp, col=col, ...)
    if(ratio>0.001) {
        add.scatter.eig(x$eig, ncol(x$disc.func), axes[1], axes[2], posi=posi, ratio=ratio, csub=csub)
    }
    return(invisible())
} # end scatter.dapc






############
## assignplot
############
assignplot <- function(x, only.grp=NULL, subset=NULL, cex.lab=.75, pch=3){
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")
    if(!inherits(x, "dapc")) stop("x is not a dapc object")

    if(!is.null(only.grp)){
        only.grp <- as.character(only.grp)
        ori.grp <- as.character(x$grp)
        x$grp <- x$grp[only.grp==ori.grp]
        x$assign <- x$assign[only.grp==ori.grp]
        x$posterior <- x$posterior[only.grp==ori.grp, , drop=FALSE]
    } else if(!is.null(subset)){
        x$grp <- x$grp[subset]
        x$assign <- x$assign[subset]
        x$posterior <- x$posterior[subset, , drop=FALSE]
    }


    ##table.paint(x$posterior, col.lab=ori.grp, ...)
    ## symbols(x$posterior)


    ## FIND PLOT PARAMETERS
    n.grp <- ncol(x$posterior)
    n.ind <- nrow(x$posterior)
    Z <- t(x$posterior)
    Z <- Z[,ncol(Z):1,drop=FALSE ]

    image(x=1:n.grp, y=seq(.5, by=1, le=n.ind), Z, col=rev(heat.colors(100)), yaxt="n", ylab="", xaxt="n", xlab="Clusters")
    axis(side=1, at=1:n.grp,tick=FALSE, label=colnames(x$posterior))
    axis(side=2, at=seq(.5, by=1, le=n.ind), label=rev(rownames(x$posterior)), las=1, cex.axis=cex.lab)
    abline(h=1:n.ind, col="lightgrey")
    abline(v=seq(0.5, by=1, le=n.grp))
    box()

    newGrp <- colnames(x$posterior)
    x.real.coord <- rev(match(x$grp, newGrp))
    y.real.coord <- seq(.5, by=1, le=n.ind)

    points(x.real.coord, y.real.coord, col="deepskyblue2", pch=pch)

    return(invisible())
} # end assignplot




###############
## randtest.dapc
###############
##randtest.dapc <- function(x, nperm = 999, ...){

##} # end randtest.dapc
