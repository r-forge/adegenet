
##
## HACK TO INSTALL R PACKAGES ON WINDOWS SYSTEMS
## WITHOUT ADMIN RIGHTS
##
## Thibaut Jombart, March 2013
## tjombart@imperial.ac.uk
##
## usage:
## hackLib()
## install.packages("foo")
##
hackLib <- function(myPath="C:/Users/Public/Rlibs"){
    if(file.exists(myPath)) myPath <- paste(myPath,"Rlibs", sep="/")
    myPath <- gsub("/+", "/", myPath)
    if(!file.exists(myPath)){
        if(!dir.create(myPath)) stop(paste(myPath, "could not be created."))
    }

    .libPaths(myPath)
    cat(paste("\nR library set to:",myPath,"\n"))
    if(as.numeric(file.info(myPath)$mode)<500) warning(paste(myPath, "may not be 'writable'."))
    return(invisible())
}
