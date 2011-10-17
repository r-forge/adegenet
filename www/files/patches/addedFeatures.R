###########################
#
# Added features for
# adegenet - not part of 
# the official release.
#
# T. Jombart
###########################


##########################
## as("seqTrack", "graphNEL")
##########################
if(require(graph)){
setOldClass("seqTrack")
setAs("seqTrack", "graphNEL", def=function(from){
    if(!require(ape)) stop("package ape is required")
     if(!require(graph)) stop("package graph is required")

     ori.labels <- rownames(from)
     from <- from[!is.na(from$ances),,drop=FALSE]


      ## CONVERT TO GRAPH
     res <- ftM2graphNEL(ft=cbind(ori.labels[from$ances], ori.labels[from$id]), W=from$weight, edgemode = "directed", V=ori.labels)
     return(res)
 })
}

