########
## dapc
########
dapc <- function(x, pop=NULL, n.pca=NULL, n.da=NULL, scale=TRUE,
                 scale.method=c("sigma", "binom"), truenames=TRUE, all.contrib=FALSE){

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

    pcaX <- dudi.pca(X, center = FALSE, scale = FALSE, scannf = FALSE, nf=maxRank)
    ## select the number of retained PC for PCA
    if(is.null(n.pca)){
        cumVar <- 100 * cumsum(pcaX$eig)/sum(pcaX$eig)
        plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative genetic variance (%)", main="Genetic variance explained by PCA")
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
    ldaX <- lda(XU, pop.fac)
    if(is.null(n.da)){
        barplot(ldaX$svd^2, xlab="Linear Discriminants", ylab="F-statistic", main="Discriminant analysis eigenvalues", col=heat.colors(length(levels(pop.fac))) )
        cat("Choose the number discriminant functions to retain (>=1): ")
        n.da <- as.integer(readLines(n = 1))
    }

    predX <- predict(ldaX, dimen=n.da)


    ## BUILD RESULT
    res <- list()
    res$tab <- XU
    res$fac <- pop.fac
    res$var.gen <- XU.lambda
    res$eig <- ldaX$svd^2
    res$disc.func <- ldaX$scaling[, 1:n.da, drop=FALSE]
    res$ind.coord <-predX$x
    res$pop.coord <- apply(res$ind.coord, 2, tapply, pop.fac, mean)
    res$prior <- ldaX$prior
    res$posterior <- predX$posterior
    res$assign <- predX$class

    ## optional: get loadings of alleles
    if(all.contrib){
        res$all.contr <- as.matrix(U) %*% as.matrix(ldaX$scaling)
        res$all.contr <- t(apply(all.contr, 1, function(e) e*e / sum(e*e)))
    }

    res$call <- match.call()
    class(res) <- "dapc"
    return(res)
} # end dapc






######################
# Function print.dapc
######################
print.dapc <- function(x, ...){
  cat("\t#########################################\n")
  cat("\t# Discriminant analysis of Principal Components #\n")
  cat("\t#########################################\n")
  cat("class: ")
  cat(class(x))
  cat("\n$call: ")
  print(x$call)
  cat("\n$n.pca:", x$n.pca, "first PCs of PCA used")
  cat("\n$n.da:", x$n.da, "discriminant functions saved")
  cat("\nProportion of conserved genetic variance: ", x$var.gen)

  cat("\nEigenvalues: ")
  l0 <- sum(x$eig >= 0)
  cat(signif(x$eig, 4)[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...\n")

  ## vectors
  sumry <- array("", c(4, 3), list(1:4, c("vector", "length", "content")))
  sumry[1, ] <- c('$eig', length(x$eig),  'eigenvalues')
  sumry[2, ] <- c('$fac', length(x$fac), 'prior group assignment')
  sumry[3, ] <- c('$prior', length(x$prior), 'prior group probabilities')
  sumry[4, ] <- c('$assign', length(x$assign), 'posterior group assignment')
  class(sumry) <- "table"
  print(sumry)

  ## data.frames
  cat("\n")
  sumry <- array("", c(5, 4), list(1:5, c("data.frame", "nrow", "ncol", "content")))
  sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "retained principal components of PCA")
  sumry[2, ] <- c("$disc.func", nrow(x$disc.func), ncol(x$disc.func), "discriminant functions")
  sumry[3, ] <- c("$ind.coord", nrow(x$ind.coord), ncol(x$ind.coord), "coordinates of individuals")
  sumry[4, ] <- c("$pop.coord", nrow(x$pop.coord), ncol(x$pop.coord), "coordinates of populations")
  sumry[5, ] <- c("$posterior", nrow(x$posterior), ncol(x$posterior), "posterior membership probabilities")
 class(sumry) <- "table"
  print(sumry)

  cat("\nother elements: ")
  if (length(names(x)) > 11)
    cat(names(x)[12:(length(names(x)))], "\n")
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
    res$n.pop <- length(levels(x$fac))

    ## assignment success
    temp <- as.character(x$fac)==as.character(x$assign)
    res$assign.prop <- mean(temp)
    res$assign.per.pop <- tapply(temp, x$fac, mean)

    return(res)
} # end summary.dapc






##############
## scatter.dapc
##############
scatter.dapc <- function(x, axes=1:2, col=rainbow(length(levels(x$fac))), posi="bottomleft", bg="lightgrey", ...){
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")
    par(bg=bg)
    s.class(x$ind.coord[,axes], fac=x$fac, col=col, ...)
    add.scatter.eig(x$eig, ncol(x$disc.func), axes[1], axes[2], posi=posi)
    return(invisible())
} # end scatter.dapc














########################
# Function summary.dapc
########################
summary.dapc <- function (object, ..., printres=TRUE) {
  if (!inherits(object, "dapc"))stop("to be used with 'dapc' object")
  if(!require(ade4,quietly=TRUE)) stop("the library ade4 is required; please install this package")
  if(!require(spdep,quietly=TRUE)) stop("the library spdep is required; please install this package")

  #util <- function(n) { ## no longer used
  #  x <- "1"
  #  for (i in 2:n) x[i] <- paste(x[i - 1], i, sep = "+")
  #  return(x)
  #}
  norm.w <- function(X, w) {
    f2 <- function(v) sum(v * v * w)/sum(w)
    norm <- apply(X, 2, f2)
    return(norm)
  }

  resfin <- list()

  if(printres) {
    cat("\nSpatial principal component analysis\n")
    cat("\nCall: ")
    print(object$call)
  }

  appel <- as.list(object$call)
  ## compute original pca
  # prepare data
  obj <- eval(appel$obj)
  if(is.null(appel$truenames)) truenames <- FALSE

  f1 <- function(vec){
    m <- mean(vec,na.rm=TRUE)
    vec[is.na(vec)] <- m
    return(vec)
  }

  if(is.genind(obj)) { X <- obj@tab }
  if(is.genpop(obj)) { X <- makefreq(obj, quiet=TRUE)$tab }

  X <- apply(X,2,f1)

  if(truenames){
    rownames(X) <- rownames(truenames(obj))
    colnames(X) <- colnames(truenames(obj))
  }

  nfposi <- object$nfposi
  nfnega <- object$nfnega

  dudi <- dudi.pca(X, center=TRUE, scale=FALSE, scannf=FALSE, nf=nfposi+nfnega)
  ## end of pca

  lw <- object$lw

  # I0, Imin, Imax
  n <- nrow(X)
  I0 <- -1/(n-1)
  L <- listw2mat(lw)
  eigL <- eigen(0.5*(L+t(L)))$values
  Imin <- min(eigL)
  Imax <- max(eigL)
  Ival <- data.frame(I0=I0,Imin=Imin,Imax=Imax)
  row.names(Ival) <- ""
  if(printres) {
    cat("\nConnection network statistics:\n")
    print(Ival)
  }

  Istat <- c(I0,Imin,Imax)
  names(Istat) <- c("I0","Imin","Imax")
  resfin$Istat <- Istat


  # les scores de l'analyse de base
  nf <- dudi$nf
  eig <- dudi$eig[1:nf]
  cum <- cumsum(dudi$eig)[1:nf]
  ratio <- cum/sum(dudi$eig)
  w <- apply(dudi$l1,2,lag.listw,x=lw)
  moran <- apply(w*as.matrix(dudi$l1)*dudi$lw,2,sum)
  res <- data.frame(var=eig,cum=cum,ratio=ratio, moran=moran)
  row.names(res) <- paste("Axis",1:nf)
  if(printres) {
    cat("\nScores from the centred PCA\n")
    print(res)
  }

  resfin$pca <- res


  # les scores de l'analyse spatiale
  # on recalcule l'objet en gardant tous les axes
  eig <- object$eig
  nfposimax <- sum(eig > 0)
  nfnegamax <- sum(eig < 0)

  ms <- multispati(dudi=dudi, listw=lw, scannf=FALSE,
                     nfposi=nfposimax, nfnega=nfnegamax)

  ndim <- dudi$rank
  nf <- nfposi + nfnega
  agarder <- c(1:nfposi,if (nfnega>0) (ndim-nfnega+1):ndim)
  varspa <- norm.w(ms$li,dudi$lw)
  moran <- apply(as.matrix(ms$li)*as.matrix(ms$ls)*dudi$lw,2,sum)
  res <- data.frame(eig=eig,var=varspa,moran=moran/varspa)
  row.names(res) <- paste("Axis",1:length(eig))

  if(printres) {
    cat("\ndapc eigenvalues decomposition:\n")
    print(res[agarder,])
  }

  resfin$dapc <- res

  return(invisible(resfin))
}



#####################
# Function plot.dapc
#####################
plot.dapc <- function (x, axis = 1, useLag=FALSE, ...){
    if (!inherits(x, "dapc")) stop("Use only with 'dapc' objects.")

    if(!require(ade4)) stop("ade4 package is required.")
    if(!require(spdep)) stop("spdep package is required.")
    if(axis>ncol(x$li)) stop("wrong axis required.")

    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mar = rep(.1,4), mfrow=c(3,2))

    n <- nrow(x$li)
    xy <- x$xy

    ## handle useLag argument
    if(useLag){
        z <- x$ls[,axis]
    } else {
        z <- x$li[,axis]
    } # end if useLag
    nfposi <- x$nfposi
    nfnega <- x$nfnega
    ## handle neig parameter - hide cn if nore than 100 links
    nLinks <- sum(card(x$lw$neighbours))
    if(nLinks < 500) {
        neig <- nb2neig(x$lw$neighbours)
    } else {
        neig <- NULL
    }

    sub <- paste("Score",axis)
    csub <- 2

    # 1
    if(n<30) clab <- 1 else clab <- 0
    s.label(xy, clab=clab, include.ori=FALSE, addaxes=FALSE, neig=neig,
            cneig=1, sub="Connection network", csub=2)

    # 2
    s.image(xy,z, include.ori=FALSE, grid=TRUE, kgrid=10, cgrid=1,
            sub=sub, csub=csub, possub="bottomleft")
    box()

    # 3
    if(n<30) {neig <- nb2neig(x$lw$neighbours)} else {neig <- NULL}
    s.value(xy,z, include.ori=FALSE, addaxes=FALSE, clegend=0, csize=.6,
            neig=neig, sub=sub, csub=csub, possub="bottomleft")

    # 4
    s.value(xy,z, include.ori=FALSE, addaxes=FALSE, clegend=0, csize=.6,
            method="greylevel", neig=neig, sub=sub, csub=csub, possub="bottomleft")

    # 5
    omar <- par("mar")
    par(mar = c(0.8, 2.8, 0.8, 0.8))
    m <- length(x$eig)
    col.w <- rep("white", m) # elles sont toutes blanches
    col.w[1:nfposi] <- "grey"
    if (nfnega>0) {col.w[m:(m-nfnega+1)] <- "grey"}
    j <- axis
    if (j>nfposi) {j <- j-nfposi +m -nfnega}
    col.w[j] <- "black"
    barplot(x$eig, col = col.w)
    scatterutil.sub(cha ="Eigenvalues", csub = 2.5, possub = "topright")
    par(mar=rep(.1,4))
    box()
    par(mar=omar)

    # 6
    par(mar=c(4,4,2,1))
    screeplot(x,main="Eigenvalues decomposition")
    par(mar=rep(.1,4))
    box()
    return(invisible(match.call()))
}

