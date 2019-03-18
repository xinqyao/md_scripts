## For plotting PCA results from cpptraj
## Assume output files of cpptraj: covar.dat, project.dat

############################################################################
################ Adjustable Parameters ##################################### 
## # bins along each dimension for density calculations
breaks <- 120

## # cores to use; default=NULL means use all available cores
ncore <- NULL

## pdb files of initial conformations of simulations; NULL means ignore
## Example: 
##   pdbfiles <- c('XXX', 'YYY', ... )
pdbfiles <- NULL

## Define boundaries of simulations; NULL means single trajectory
## Example:
##                     start,  end
##   bounds <- matrix(c(1,     XXX, 
##                      XXX+1, YYY,
##                          ...   ), ncol=2, byrow=TRUE)
bounds <- matrix(c(1,     10000,
                   10001, 20000,
                   20001, 30000), ncol=2, byrow=TRUE)

## pdb file of averged conformation from pca; NULL, ignored
avgpdb <- NULL

## color for trajectories; Recommend lighter colors
col.traj <- c("lightblue", "pink", "lightgreen")

## color for contour lines and outlines; Recommend darker colors
col.line <- c("blue", "red", "darkgreen")

## plot contour line or not?
contour <- TRUE

## Should the last trajectory be outlined?
outline.last <- TRUE

## customized xlim and ylim; NULL, ignored
## Example:
##   myxlim <- c(-30, 40)
##   myylim <- c(-40, 60)
myxlim <- NULL
myylim <- NULL

############## END Adjustable Parameters ################################# 
##########################################################################


## !! DO NOT CHANGE FOLLOWING LINES UNLESS YOU KNOW WHAT IT MEANS !! ##

library(bio3d)
library(ggplot2)
library(parallel)

## A function to convert variance-covariance matrix into a bio3d 'pca'-like object
mypca <- function(S, z, avg=NULL) {
  prj  <- eigen(S, symmetric = TRUE)
  L <- prj$values
  U <- prj$vectors

  L[L<0]<-0
  sdev <- sqrt(L)
  if(ncol(U) %% 3 == 0) {
     au <- apply(U, 2, function(x) {
       sqrt(colSums(matrix(x^2, nrow=3))) })
  } else {
     au <- NULL
  }

  class(U)="pca.loadings"
  out <- list(L=L, U=U, z=z, au=au, sdev=sdev, mean=avg)

  class(out)="pca"
  out
}

## a function copied from the 'entropy' package
discretize2d <-
   function (x1, x2, numBins1, numBins2, r1 = range(x1), r2 = range(x2))
   {
       b1 = seq(from = r1[1], to = r1[2], length.out = numBins1 +
           1)
       b2 = seq(from = r2[1], to = r2[2], length.out = numBins2 +
           1)
       y2d = table(cut(x1, breaks = b1, include.lowest = TRUE),
           cut(x2, breaks = b2, include.lowest = TRUE))
       return(y2d)
   }


######## START   CALCULATIONS #################

ncore <- setup.ncore(ncore)
vcov <- read.table('covar.dat')
z <- read.table('project.dat')

if(!is.null(pdbfiles)) {
   pdb.list <- mclapply(pdbfiles, read.pdb, mc.cores=ncore)
} else {
   pdb.list <- NULL
}

if(is.null(bounds)) {
   bounds <- matrix(c(1, nrow(z)), ncol=2, byrow=TRUE)
}

if(!is.null(avgpdb)) {
   avg <- read.pdb(avgpdb)
   avg <- trim(avg, "protein", elety=c("N", "CA", "C", "O")) ## assume PCA is performed on mainchain
} else {
   avg <- NULL
}

if(length(col.traj) != nrow(bounds)) {
   stop("Colors do not match trajectories")
}

pc <- mypca(S=vcov, z=as.matrix(z[, 2:4]), avg=as.vector(avg$xyz))
rm(z, vcov)

# check signs of eigenvectors
lines <- readLines('evecs.dat')
inds <- which(lines == " ****")
myU <- matrix(NA, nrow=nrow(pc$U), ncol=3)
for(i in 1:3) {
   myinds <- c((inds[i]+2) : (inds[i+1]-1))
   myU[, i] <- as.numeric(unlist(strsplit(trimws(lines[myinds]), split="\\s+")))
}
chk <- round(t(pc$U[, 1:3]) %*% myU, 3)
if(! (all(chk[upper.tri(chk)]==0) || all(abs(diag(chk))==1)) ) {
   stop("Calculated eigenvectors do not match those obtained from cpptraj")
}
pc$U[, 1:3] <- t( t(pc$U[, 1:3]) * diag(chk) )
#pc$z[, 1:3] <- t( t(pc$z[, 1:3]) * diag(chk) )
######### End of check #############

vcx1 <- round(pc$L[1]/sum(pc$L)*100, 1)
vcx2 <- round(pc$L[2]/sum(pc$L)*100, 1)

z <- pc$z[, 1:2]
xr = range(z[, 1])
yr = range(z[, 2])

xlim = c(xr[1] - diff(xr)*0.05, xr[2] + diff(xr)*0.05)
ylim = c(yr[1] - diff(yr)*0.05, yr[2] + diff(yr)*0.05)
xlab = paste('PC1 (', format(vcx1, nsmall=1),'%)', sep='')
ylab = paste('PC2 (', format(vcx2, nsmall=1),'%)', sep='')

## build data for each simulation
rets <- mclapply(1:nrow(bounds), function(i) {
   data <- discretize2d(z[bounds[i, 1]:bounds[i, 2], 1], z[bounds[i, 1]:bounds[i, 2], 2], 
      breaks, breaks, r1=xlim, r2=ylim)
   xx = seq(xlim[1], xlim[2], length.out = breaks + 1)
   yy = seq(ylim[1], ylim[2], length.out = breaks + 1)
   xx = xx[-length(xx)]
   yy = yy[-length(yy)]
   data <- data / sum(data) / ((diff(xlim)/breaks) * (diff(ylim)/breaks))
   class(data) <- "matrix"
   datab <- data
   data[data > 0] <- 1
   
   ## arbitary plot() for xspline()
   plot(1)
   
   inds <- matrix(which(data>0, arr.ind=TRUE), ncol=2)
   ch <- chull(inds)
   ch <- xspline(xx[inds[ch, 1]], yy[inds[ch, 2]], shape=0.3, open=FALSE, draw=FALSE)
   ch <- do.call(cbind, ch)
   dd <- data.frame(xx=xx[inds[, 1]], 
                    yy=yy[inds[, 2]], 
                    density=as.numeric(datab[datab>0]))

   ## for annotation of initial conformations
   if(!is.null(pdb.list)) {
      pdb <- pdb.list[[i]]
      pdb <- trim(pdb, "protein", elety=c("N", "CA", "C", "O"))  ## assume PCA is performed on mainchain
      if(!is.null(avg)) {
         myz <- project.pca(pdb$xyz, pc, fit=TRUE, 
                            fixed.inds=1:ncol(pdb$xyz), mobile.inds=1:ncol(pdb$xyz))
      } else {
         warning("Average pdb file not provided. Annotation of initial conformations ignored")
         myz <- NULL
      }
   } else {
      myz <- NULL
   }
   list(ch=ch, dd=dd, myz=myz)
}, mc.cores=1)  ## has to be 1; otherwise, errors because of xspline() call!!
dev.off()
   
############ shaded area, contour lines, outlines, and annotated initial conformations (optional) ###########
pdf(onefile=TRUE, file='pca.pdf', width=3, height=3)

p <- ggplot(rets[[1]]$dd, aes(xx, yy)) + geom_raster(fill=col.traj[1])
if(length(rets) > 1) {
  for(i in 2:length(rets)) {
     p <- p + geom_raster(data=rets[[i]]$dd, aes(x=xx, y=yy), fill=col.traj[i])
  }
  to.outline <- 1:( ifelse(outline.last, length(rets), length(rets)-1) )
  for(i in to.outline) {
     p <- p + geom_polygon(data=data.frame(x=rets[[i]]$ch[, 1], y=rets[[i]]$ch[, 2]), 
        aes(x, y), fill=NA, col=col.line[i], size=0.4)
  }
}
if(contour) {
  p <- p + geom_contour(aes(z=density), col=col.line[1], size=0.4)
  if(length(rets) > 1) {
     for(i in 2:length(rets)) {
       p <- p + geom_contour(data=rets[[i]]$dd, aes(x=xx, y=yy, z=density), col=col.line[i], size=0.4)
     }
  }
}
inds <- which(sapply(rets, function(x) !is.null(x$myz)))
if(length(inds)>0) {
   xx <- sapply(rets[inds], function(x) x$myz[1])
   yy <- sapply(rets[inds], function(x) x$myz[2])
   p <- p + annotate("point", x=xx, y=yy, color=col.line[inds], size=0.5)
}
if(!is.null(myxlim)) {
   xlim <- myxlim
}
if(!is.null(myylim)) {
   ylim <- myylim
}
p <- p + xlab(xlab) + ylab(ylab) + 
#         xlim(-50, 30) + ylim(-40, 40)+
         xlim(xlim[1], xlim[2]) + ylim(ylim[1], ylim[2])+
         theme_bw()

print(p)

dev.off()


## scree plot
pdf(onefile=TRUE, file='scree.pdf', width=4, height=4)
plot.pca.scree(pc)
dev.off()

## PC trajectories
if(!is.null(avg)) {
   mktrj(pc, pc=1, pdb=avg, file="pc1.pdb")
   mktrj(pc, pc=2, pdb=avg, file="pc2.pdb")
}

## Time-series of PCs
max.frame <- 10000   # max. # frames to plot
# time per frame (Unit: ns)
tpf <- 0.001

dat <- mclapply(1:nrow(bounds), function(i) {
   pc$z[bounds[i, 1]:bounds[i, 2], ]
}, mc.cores=ncore)

ll <- sapply(dat, nrow)

## reduce plot if data points are too many
if(max(ll) > max.frame) {
   fac <- ceiling(max(ll) / max.frame)
   dat <- mclapply(dat, function(x) x[seq(1, nrow(x), fac), ], mc.cores=ncore)
   tpf <- tpf * fac
}

## pretty x-axis ticks
xx <- seq_along(dat[[which.max(ll)]][, 1]) * tpf
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

pdf(onefile=TRUE, file='ts-pcs.pdf', width=4, height=6)
layout(matrix(1:3, nrow=3, ncol=1))

par(mar=c(2,4,1,1))
plot(x=xx[1:nrow(dat[[1]])], y=dat[[1]][, 1], typ='l', col=col.line[1], xaxt="n",
     xlim=range(xx), ylim=range(unlist(sapply(dat, "[", ,1))),
     xlab=paste("Time (", xunit, ")", sep=""),
     ylab="PC1 (A)")
axis(1, at=xticks)
if(length(dat) > 1) {
   for(i in 2:length(dat)) {
      lines(x=xx[1:nrow(dat[[i]])], y=dat[[i]][, 1], col=col.line[i])
   }
}

par(mar=c(2,4,1,1))
plot(x=xx[1:nrow(dat[[1]])], y=dat[[1]][, 2], typ='l', col=col.line[1], xaxt="n",
     xlim=range(xx), ylim=range(unlist(sapply(dat, "[", ,2))),
     xlab=paste("Time (", xunit, ")", sep=""),
     ylab="PC2 (A)")
axis(1, at=xticks)
if(length(dat) > 1) {
   for(i in 2:length(dat)) {
      lines(x=xx[1:nrow(dat[[i]])], y=dat[[i]][, 2], col=col.line[i])
   }
}

par(mar=c(4,4,1,1))
plot(x=xx[1:nrow(dat[[1]])], y=dat[[1]][, 3], typ='l', col=col.line[1], xaxt="n",
     xlim=range(xx), ylim=range(unlist(sapply(dat, "[", ,3))),
     xlab=paste("Time (", xunit, ")", sep=""),
     ylab="PC3 (A)")
axis(1, at=xticks)
if(length(dat) > 1) {
   for(i in 2:length(dat)) {
      lines(x=xx[1:nrow(dat[[i]])], y=dat[[i]][, 3], col=col.line[i])
   }
}

dev.off()

