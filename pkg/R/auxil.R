###########################
#
# Auxiliary functions for
# adegenet objects
#
# T. Jombart
###########################


#######################
# Function rmspaces
#######################
# removes spaces and tab at the begining and the end of each element of charvec
.rmspaces <- function(charvec){
    charvec <- gsub("^([[:blank:]]*)([[:space:]]*)","",charvec)
    charvec <- gsub("([[:blank:]]*)([[:space:]]*)$","",charvec)
    return(charvec)
}





###################
# Function readExt
###################
.readExt <- function(char){
    temp <- as.character(char)
    temp <- unlist(strsplit(char,"[.]"))
    res <- temp[length(temp)]
    return(res)
}





###################
# Function .genlab
###################
# recursive function to have labels of constant length
# base = a character string
# n = number of labels
.genlab <- function(base, n) {
  f1 <- function(cha,n){
    if(nchar(cha)<n){
      cha <- paste("0",cha,sep="")
      return(f1(cha,n))
    } else {return(cha)}
  }
  w <- as.character(1:n)
  max0 <- max(nchar(w))
  w <- sapply(w, function(cha) f1(cha,max0))
  return(paste(base,w,sep=""))
}





#######################
# Function adegenetWeb
#######################
adegenetWeb <- function(){
    cat("Opening url \"http://adegenet.r-forge.r-project.org/\" ...\n")
    browseURL("http://adegenet.r-forge.r-project.org/")
}





############################
# Function adegenetTutorial
############################
adegenetTutorial <- function(which=c("general","spca")){
    which <- match.arg(which)
    if(which=="general"){
        url <- "http://adegenet.r-forge.r-project.org/files/adegenet.pdf"
        cat("\n")
        cat("  >> Seeking the general tutorial for adegenet.\n")
        cat("  >> Opening url \"",url,"\".\n ", sep="")
        cat("\n")
        browseURL(url)
    }
    if(which=="spca"){
        url <- "http://adegenet.r-forge.r-project.org/files/tutorial-spca.pdf"
        cat("\n")
        cat("  >> Seeking the sPCA tutorial for adegenet.\n")
        cat("  >> Opening url \"",url,"\". \n", sep="")
        cat("\n")
        browseURL(url)
    }
}





############
# checkType
############
checkType <- function(markType=x@type){

    if(markType=="codom") return() # always ok for codominant markers

    currCall <- match.call()
    currFunction <- sub("[[:space:]]*[(].*","",currCall)

    ## names of functions which are ok for dominant markers
    PAOk <- c("genind","genpop","genind2genpop","summary",
                 "truenames","seppop","na.replace","nLoc")

    PAWarn <- c("df2genind")

    ## function exists but is experimental
    if(currFunction %in% PAWarn){
        msg <- paste(currFunction,"is implemented but experimental presence/absence markers")
        warning(msg)
        return()
    }

    ## function not implemented
    if(! currFunction %in% PAOk){
        msgError <- paste(currFunction,"is not implemented for presence/absence markers")
        stop(msgError)
    } else return() # else, ok.
} # end checkType

