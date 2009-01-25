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
checkType <- function(markType){
    if(markType=="codom") return()

    currCall <- match.call()
    currFunction <- sub("[[:space:]]*[(].*","",currCall)

    ## names of functions which are ok for dominant markers
    dominOk <- c("genind","genpop","genind2genpop","na.replace","nLoc")

    if(! currFunction %in% dominOk){
        msgError <- paste(currFunction,"is not implemented for dominant markers")
        stop(msgError)
    }
} # end checkType

