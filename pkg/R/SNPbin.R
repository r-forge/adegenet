## each byte takes a value on [0,255]

## function to code multiple SNPs on a byte
## 7 combinations of SNPs can be coded on a single byte
## we use bytes values from [1,128]
## 200 is a missing value
SNPCOMB <- expand.grid(rep(list(c(0,1)), 7))


f1 <- function(vecSnp){
    nbBytes <- 1+ length(vecSnp) %/% 7
    out.length <- 7*nbBytes
    temp <- c(vecSnp, rep(0, out.length-length(vecSnp))) # fill the end with 0 of necessary
    sapply(seq(1, by=7, length=nbBytes), function(i) which(apply(SNPCOMB,1, function(e) all(vecSnp[i:(i+7)]==e))) )


    as.raw(length(vecSnp))
}
