testExamples <- function(dir="~/dev/adegenet/pkg"){
    setwd(dir)
    setwd("man")
    toRead <- dir()
    for(e in toRead){
        txt <- readLines(e)
        temp <- txt[grep("alias",txt)]
        
    }
}
