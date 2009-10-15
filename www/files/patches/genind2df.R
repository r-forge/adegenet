#####################
# Function genind2df
#####################
genind2df <- function(x, pop=NULL, sep="", usepop=TRUE, oneColPerAll=FALSE){

  if(!is.genind(x)) stop("x is not a valid genind object")
  ## checkType(x)

  if(is.null(pop)) {
      pop <- x@pop
      levels(pop) <- x@pop.names
  }

  if(oneColPerAll){
      sep <- "/"
  }

  ## PA case ##
  if(x@type=="PA"){
      temp <- truenames(x)
      if(is.list(temp) & usepop){
          res <- cbind.data.frame(pop=temp[[2]],temp[[1]])
      } else{
          if(is.list(temp)) {
              res <- temp[[1]]
          } else{
              res <- temp
          }
      }

      return(res) # exit here
  }

  ## codom case ##
  # make one table by locus from x@tab
  kX <- seploc(x,res.type="matrix")
  kX <- lapply(kX, function(X) round(X*x@ploidy)) # take data as numbers of alleles
  ## (kX is a list of nloc tables)

  ## function to recode a genotype in form "A1[sep]...[sep]Ak" from frequencies
  recod <- function(vec,lab){
      if(any(is.na(vec))) return(NA)
      res <- paste( rep(lab,vec), collapse=sep)
      return(res)
  }


  # kGen is a list of nloc vectors of genotypes
  kGen <- lapply(1:length(kX), function(i) apply(kX[[i]],1,recod,x@all.names[[i]]))
  names(kGen) <- x@loc.names

  ## if use one column per allele
  if(oneColPerAll){
      f1 <- function(vec){ # to repeat NA with seperators
          vec[is.na(vec)] <- paste(rep("NA",x@ploidy), collapse=sep)
          return(vec)
      }
      temp <- lapply(kGen, f1)
      temp <- lapply(temp, strsplit,sep)

      res <- lapply(temp, function(e) matrix(unlist(e), ncol=x@ploidy, byrow=TRUE))
      res <- data.frame(res,stringsAsFactors=FALSE)
      names(res) <- paste(rep(locNames(x),each=x@ploidy), 1:x@ploidy, sep=".")

      ## handle pop here
      if(!is.null(pop) & usepop) res <- cbind.data.frame(pop,res,stringsAsFactors=FALSE)

      return(res) # exit here
  } # end if oneColPerAll

  ## build the final data.frame
  res <- cbind.data.frame(kGen,stringsAsFactors=FALSE)

  ## handle pop here
  if(!is.null(pop) & usepop) res <- cbind.data.frame(pop,res,stringsAsFactors=FALSE)

  return(res)
}
