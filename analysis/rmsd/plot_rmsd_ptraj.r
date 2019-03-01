fin <- "rms.out" # rms file generated by cpptraj
fout <- "rmsd.pdf"
max.frame <- 10000   # max. # frames to plot

# time per frame (Unit: ns)
tpf <- 0.001


dat <- read.table(fin)

## reduce plot if data points are too many
if(nrow(dat) > max.frame) {
   fac <- ceiling(nrow(dat) / max.frame)
   dat <- dat[seq(1, nrow(dat), fac), ]
   tpf <- tpf * fac
}

## pretty x-axis ticks
xx <- seq_along(dat[, 1]) * tpf
xunit <- "ns"
if(max(xx) > 1000) {
   tpf <- tpf / 1000
   xx <- xx / 1000
   xunit <- "us"
} else if(max(xx) < 1) {
   tpf <- tpf * 1000
   xx <- xx * 1000
   xunit <- "ps"
}
xticks <- pretty(xx, n=10)

pdf(fout, width=3.23, height=3, colormodel="srgb", pointsize=10)
plot(x=xx, y=dat[, 2], typ='l', col="black", xaxt="n", 
     xlab=paste("Time (", xunit, ")", sep=""),
     ylab="RMSD (A)")
axis(1, at=xticks)
axis(2, at=round(dat[1, 2], 1), col="red")
dev.off()
