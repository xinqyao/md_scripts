#source("~/bio3d/new_funs/myfuns.R")

#"slideMean" <- function(x, windowsize=3, slide=1) {
#   idx1 <- seq(1, length(x)-windowsize+1, by=slide)
#   idx2 <- idx1 + windowsize
#   cx <- c(0,cumsum(x))
#   return((cx[idx2]-cx[idx1])/windowsize)
#}

## more balanced
"slideMean" <- function(x, windowsize=3, slide=1) {
   ws <- seq(1, length(x), 2)
   windowsize <- ws[ which.min(abs(windowsize - ws)) ]
#   idx1 <- seq(1, length(x)-windowsize+1, by=slide)
   idx1 <- seq(1, length(x), by=slide) - (windowsize-1)/2
   idx2 <- idx1 + windowsize 
   idx1[idx1 < 1] <- 1
   idx2[idx2 > (length(x)+1)] <- length(x)+1
   len <- idx2 - idx1 
   cx <- c(0,cumsum(x))
   return((cx[idx2]-cx[idx1])/len)
}

plot_all <- function(fin="./summary.EPTOT",
      fout="eptot.png", xlab="Time (ns)",
      ylab="Ep (kcal/mol)", zoom=FALSE) {
  
   mydata <- read.table(fin)
   max.frame <- 10000 # max. data points to plot
   if(nrow(mydata) > max.frame) {
      fac <- ceiling(nrow(mydata) / max.frame)
      mydata <- mydata[seq(1, nrow(mydata), fac), ]
   }
   
   ## pretty x-axis ticks
   ## Assume time unit ps: True for AMBER output file 
   xx <- mydata$V1/1000
   xlab <- "Time (ns)"
   if(max(xx) > 1000) {
      xx <- xx / 1000
      xlab <- "Time (us)"
   } else if(max(xx) < 1) {
      xx <- xx * 1000
      xlab <- "Time (ps)"
   }
   xticks <- pretty(xx, n=10)

   png(fout, width=8.2, height=8.2, unit="cm", 
           res=300, pointsize=8)
   if(zoom) {
      mm <- mean(mydata$V2[(nrow(mydata)-9):nrow(mydata)])
      ss <- sd(mydata$V2[(nrow(mydata)-9):nrow(mydata)])
      plot(xx, mydata$V2, type="l", xaxt="n", col="black", xlab=xlab, ylab=ylab, 
           ylim=c(min(mydata$V2), mm+5*ss)) 
   } else {
      plot(xx, mydata$V2, type="l", xaxt="n", col="black", xlab=xlab, ylab=ylab)
   }
   axis(1, at=xticks)
   ws <- max(floor(length(mydata$V1)*0.1), 3)
#   lines(mydata$V1[1:(length(mydata$V1)-ws+1)]/1000, 
   lines(xx, slideMean(mydata$V2, windowsize=ws), type="l", col="red")
   dev.off()
}

plot_all(fin="./summary.ETOT", fout="etot.png", ylab="Etot (kcal/mol)")
plot_all(fin="./summary.EPTOT", fout="eptot.png", ylab="Ep (kcal/mol)")
plot_all(fin="./summary.DIHEDRAL", fout="edihe.png", ylab="Edihe (kcal/mol)")
plot_all(fin="./summary.VOLUME", fout="volume.png", ylab="Volume (A^3)")
# zoom in VOLUME
plot_all(fin="./summary.VOLUME", fout="volume2.png", ylab="Volume (A^3)", zoom=TRUE)
#nline <- length(readLines("./summary.VOLUME"))
#system(paste("sed -n '", floor(nline/10), ",$'p summary.VOLUME > summary2.VOLUME", sep=""))
#plot_all(fin="./summary2.VOLUME", fout="volume2.png", ylab="Volume (A^3)")
plot_all(fin="./summary.DENSITY", fout="density.png", ylab="Density (g/mL)")
plot_all(fin="./summary.PRES", fout="pres.png", ylab="Pressure (bar)")
plot_all(fin="./summary.TEMP", fout="temp.png", ylab="Temperature (K)")
