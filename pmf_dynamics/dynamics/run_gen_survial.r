library(bio3d)
threslife <- 20
outdir <- paste0("ts_time_life", threslife)
if(!dir.exists(outdir)) {
   dir.create(outdir)
}

for(para in c("3.5", "4.2", "4.9", "5.6")) {
   print(para)
   files <- list.files(para, 'ang.*.out', full.names=TRUE)
   ofile <- paste(outdir, "/ts_time_", para, ".txt", sep="")
   cat("", file=ofile)
   #tmpfile <- tempfile()
   for (i in files) {
   #  print(i)
     dat <- read.table(i)

     # start from cis (0)
#     n <- which(dat$V2 > 90)
#     m <- which(dat$V2 < -90)

     # start from trans (180)
     n <- which(dat$V2 < 90)
     m <- which(dat$V2 > 270)

     if(length(n) > 0) {
        segn <- bounds(n)
        segn <- segn[segn[, "length"] >= threslife, ,drop=FALSE]
     }
     else {
        segn <- matrix(integer(0), ncol=3)
     }
     if(length(m) > 0) {
        segm <- bounds(m)
        segm <- segm[segm[, "length"] >= threslife, ,drop=FALSE]
     }
     else {
        segm <- matrix(integer(0), ncol=3)
     }
 
     if(nrow(segn)>0) {
        if(nrow(segm)==0 || (nrow(segm)>0 && segm[1, "start"]>segn[1, "start"])) {
           cat(segn[1, "start"], " 1\n", file=ofile, append=TRUE)
   #        x11()
   #        par(mfrow = c(2, 1))
   #        plot(dat$V2, typ='p', pch=20, cex=0.3)
   #        abline(v=c(n[1], m[1]), lty=1, col='red')
   #        plot(dat$V2, typ='p', pch=20, xlim=c(n[1]-100, n[1]+100), cex=0.3)
   ##        title(main=paste('life=',life, ", type=", type, sep=""))
   #        abline(v=c(n[1]), lty=1, col='red')
        }
        else {
           cat(segm[1, "start"], " 0\n", file=ofile, append=TRUE)
        }
     }
     else if(nrow(segm)>0) {
        cat(segm[1, "start"], " 0\n", file=ofile, append=TRUE)
     }
   }
#if(file.exists(ofile)) {
#   system(paste("paste", ofile, tmpfile, "> tmp"))
#   system(paste("mv tmp", ofile))
#}
#else {
#   file.copy(tmpfile, ofile)
#}
}
