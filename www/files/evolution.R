pdf(file="~/Desktop/adegenetEvol.pdf")

dat <- read.table(file="evolution.txt", head=T)

par(lend=2)

## plots
plot(dat$Functions, type="b", col="blue", lwd=2, xaxt="n",xlab="Versions",ylab="Number of items", main="Evolution of adegenet", ylim=c(0,110))

points(dat$DocSize, type="b", pch=3, col="red", lwd=2)

points(dat$Datasets, type="b", pch=5, col="orange", lwd=2)

## annot/legend
legend("topleft", pch=c(1,3,5), lwd=2, leg=c("Functions","Manpages","Datasets"), col=c("blue", "red", "orange"))

text(2.5,55, "Formal (S4) \nclasses", font=2)

text(5.5,85, "Different ploidy \nlevels", font=2)

## axis
axis(side=1, at=1:nrow(dat), lab=dat$Version)



dev.off()
