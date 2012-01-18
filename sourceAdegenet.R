sourceAdegenet <- function(path=getwd()){
    ## EXPECTED FILES ##
    expFiles <- c("classes.R", "basicMethods.R", "handling.R", "auxil.R", "setAs.R", "SNPbin.R", "glHandle.R", "glFunctions.R", "glSim.R", "find.clust.R", "hybridize.R", "scale.R", "fstat.R", "import.R", "seqTrack.R", "chooseCN.R", "genind2genpop.R", "loadingplot.R", "sequences.R", "gstat.randtest.R", "makefreq.R", "colorplot.R", "monmonier.R", "spca.R", "coords.monmonier.R", "haploGen.R", "old2new.R", "spca.rtests.R", "dapc.R", "haploPop.R", "PCtest.R", "dist.genpop.R", "Hs.R", "propShared.R", "export.R", "HWE.R", "propTyped.R", "inbreeding.R", "glPlot.R")

    ## READ FILES IN DIR ##
    foundFiles <- dir(path, pattern=".R", full=FALSE)
    foundPaths <- dir(path, pattern=".R", full=TRUE)


    ## CHECK FILES IN DIR ##
    temp <- expFiles %in% foundFiles
    if(all(!temp)) stop(paste("none of the expected files are found in", path))
    if(!all(temp)) {
        warning(paste("some expected source files are missing from" , path))
        cat("\nMissing files:")
        print(expFiles[!temp])
    }


    ## SOURCE FILES ##
    expFiles <- expFiles[temp]
    cat("\nLoading sources of adegenet 1.3-3...\n")
    ord <- match(expFiles,foundFiles)
    invisible(sapply(foundPaths[ord], source))
    cat("\n... done.")

    cat("\n\nNOTE! This kludge allowed you to load only source files of adegenet, without data or documentation.")
    cat("\nPlease consider installing the package for regular uses.\n")

    return(invisible())
}
