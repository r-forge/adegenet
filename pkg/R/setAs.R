#############
# S4 methods
#############
setAs("genind", "data.frame", function(from, to) {
    return(from@tab)
})



setAs("genpop", "data.frame", function(from, to) {
    return(from@tab)
})



setAs("genind", "matrix", function(from, to) {
    return(from@tab)
})



setAs("genpop", "matrix", function(from, to) {
    return(from@tab)
})



setAs("genind", "genpop", function(from, to) {
    if(!is.genind(from)) stop("object is not a valid genind")

    x <- genind2genpop(from, quiet=TRUE)
    warning("You had better use genind2genpop to specify treatment of NAs")

    return(x@tab)
})




##############
# S3 versions
##############
as.data.frame.genind <- function(x,...){
    return(as(x,"data.frame"))
}



as.data.frame.genpop <- function(x,...){
    return(as(x,"data.frame"))
}



as.matrix.genind <- function(x,...){
    return(as(x,"matrix"))
}



as.matrix.genpop <- function(x,...){
    return(as(x,"matrix"))
}



as.genpop.genind <- function(x,...){
    return(as(x,"genpop"))
}
