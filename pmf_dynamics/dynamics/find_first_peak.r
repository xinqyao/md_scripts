root <- "../us/resampling/sep1_-30_210"

tol <- 0
# tol <- c(0, 0, 0, 0)

## Which peak to use? (give the peak ID)
ipeaks <- 1
#ipeaks <- c(1, 1, 1, 1)

#library(bio3d)
dirs <- list.files(root, "^[0-9\\.]*$", full.names=TRUE)
dirs <- dirs[order(as.numeric(basename(dirs)))]

source("funs.R")
if(length(tol) == 1) {
   tol <- rep(tol, length(dirs))
}
if(length(ipeaks) == 1) {
   ipeaks <- rep(ipeaks, length(dirs))
}

cat("Root path: ", root, "\n")
for(j in 1:length(dirs)) {
   i <- dirs[j]
   cat("Processing ", basename(i), "...\n", sep="")
   files <- list.files(i, "^fep\\.[0-9]*\\.out$", full.names=TRUE)
   inds <- as.numeric(sub("fep\\.([0-9]*)\\.out", "\\1", basename(files)))
   file <- files[order(inds, decreasing=TRUE)[1]]
   dat <- read.table(file)
   myout <- get_wells_peaks(dat$V2, tol=tol[j])
  
#   ind <- which.min(dat$V2)
#   cat("Well (Gmin): ", dat$V1[ind], "\n")
#   ddat <- diff(dat$V2)
#   bs1 <- bounds(which(ddat[1:(ind-1)] < 0))
#   peak1.ind <- bs1[nrow(bs1), 1]
#   bs2 <- bounds(which(ddat[ind:length(ddat)] > 0))
#   peak2.ind <- bs2[1, 2] + ind
   cat("Peak TS: ", dat$V1[myout$peaks[ipeaks[j]]], "\n")
#   cat("Peak 2: ", dat$V1[peak2.ind], "\n")
   cat("\n\n") 
}
cat("Done.\n")
