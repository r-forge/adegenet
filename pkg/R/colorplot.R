colorplot <- function(xy, X, axes=1:ncol(X), add.plot=FALSE, defaultLevel=0, ...){

    ## some checks
    if(any(is.na(xy))) stop("NAs exist in xy")
    if(!is.numeric(xy)) stop("xy is not numeric")
    if(nrow(xy) != nrow(X)) stop("xy and X have different row numbers")
    X <- as.matrix(X[,axes,drop=FALSE])
    if(any(is.na(X))) stop("NAs exist in X")
    if(!is.numeric(X)) stop("X is not numeric")
    if(defaultLevel < 0 | defaultLevel>1) stop("defaultLevel must be between 0 and 1")

    ## function mapping x to [0,+inf[
    f1 <- function(x){
        x <- x + abs(min(x))
        return(x)
    }

    ## apply f1 to X
    X <- apply(X, 2, f1)

    v1 <- X[,1]
    if(ncol(X)==2) {v2 <- X[,2]} else {v2 <- defaultLevel}
    if(ncol(X)==3) {v2 <- X[,3]} else {v3 <- defaultLevel}

    ## find the colors
    col <- rgb(v1, v2, v3, maxColorValue=max(X))

    ## plot data
    if(!add.plot) {
        plot(xy, pch=20, col=col, ...)
    } else {
        points(xy, pch=20, col=col, ...)
    }
}
