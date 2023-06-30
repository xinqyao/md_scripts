library(ggplot2)
library(bio3d)

###### !!! MAY NEED TO UPDATE FOR "MULTI-PEAK" CASES !!! ###########

#rootpath <- "resampling/sep1_10k_err"
rootpath <- "resampling/sep1_-30_210"
b_extrapolate <- FALSE

## Set TRUE to include parameters under 'trash' for extrapolation 
## If TRUE, need to set trash path
b_extrapolate_include_trash <- FALSE

#rootpath_trash <- "trash/resampling/sep1_10k_err"
#rootpath_trash <- "trash/resampling/sep1_-30_210"
rootpath_trash <- file.path("trash", rootpath)


####################
tol <- 0
# tol <- c(0, 0, 0, 0)

## Which peak to use? (give the peak ID)
ipeaks <- 1
#ipeaks <- c(1, 1, 1, 1)

## Which well to use for cis? (give the well ID)
iwells <- 1
#iwells <- c(1, 1, 1, 1)

## Which well to use for trans? (give the well ID)
jwells <- 2
#jwells <- c(2, 2, 2, 2)
####################

datfile <- list.files(rootpath, pattern='^fep.*\\.out$', recursive=TRUE, full.names=TRUE)

if( (!b_extrapolate) & b_extrapolate_include_trash ) {
   b_extrapolate_include_trash <- FALSE
}
if(b_extrapolate_include_trash) {
   if(!is.null(rootpath_trash)) {
      datfile_trash <- list.files(rootpath_trash, pattern='^fep.*\\.out$', recursive=TRUE, full.names=TRUE)
   }
   else {
      stop("Set trash path")
   }
   datfile <- c(datfile, datfile_trash)
}
if(length(datfile)==0) {
   stop("No data file found.")
}

## Sort files
alphas <- as.numeric(basename(dirname(datfile)))
numbers <- as.numeric(sub("fep.([0-9]*).out$", "\\1", basename(datfile)))
inds <- alphas * 10^6 + numbers
datfile <- datfile[order(inds)]
alphas <- as.numeric(basename(dirname(datfile)))
####

# take the last file for each alpha to plot
rl <- rle2(alphas)
datfile <- datfile[rl$inds]

dat <- NULL
for(i in datfile) {
   tdat <- read.table(i)
#   fref <- sum(tdat[tdat[, 1]%in% c(178.75, 181.25), 2])/2 # reference to trans
   fref <- tdat[tdat[, 1]%in% 180, 2] # reference to trans
   fref2 <- tdat[tdat[, 1]%in% 0, 2] # reference to cis 
   tdat <- cbind(tdat, as.numeric(sub(".*\\/([\\.0-9]+)\\/fep.*", "\\1", i))/14.0)
   tdat[, 2] <- tdat[, 2] - (fref + fref2)/2
   dat <- rbind(dat, tdat)
}
colnames(dat) <- c("Angle", "Energy", "EnergySD", "Alpha")
dat$Alpha <- factor(format(round(dat$Alpha, 2), nsmall=2), levels=format(seq(0, 1, 0.01), nsmall=2))
dat$Alpha2 <- 1-as.numeric(as.character(dat$Alpha))

if(b_extrapolate) {
   coefs <- tapply(1:nrow(dat), dat$Angle, function(i) {
       inds <- which(dat$Alpha[i]!="1.00")
       if(length(inds)>0) {
          if(any(!is.finite(dat$EnergySD[i][inds]))) {
             summary(lm(Energy~Alpha2, data=dat[i, ][inds, ]))$coefficients
          }
          else {
             summary(lm(Energy~Alpha2, data=dat[i, ][inds, ], weights=1/(dat$EnergySD[i][inds]^2)))$coefficients
          }
       }
       else {
          NA
       }
   })

   dat.pred <- data.frame(Angle=as.numeric(names(coefs)),
      Energy=unname(sapply(coefs, "[", "(Intercept)", 1)),
      EnergySD=unname(sapply(coefs, "[", "(Intercept)", 2)),
      Alpha="1.00", Alpha2=0.0, predict=TRUE)

   dat2 <- cbind(dat, predict=FALSE)
   dat2 <- rbind(dat2, dat.pred)

   if(b_extrapolate_include_trash) {
      ## remove trash in output
      alphas_trash <- format(round(as.numeric(basename(dirname(datfile_trash)))/14.0, 2), nsmall=2)
      dat2 <- dat2[!dat2$Alpha %in% alphas_trash, ]
   }

} else {

   ## comment out if extrapolation is used
   dat2 <- dat
   dat.pred <- dat[dat$Alpha=="1.00", ]

}


source("funs.R")
if(length(tol) == 1) {
   tol <- rep(tol, length(unique(as.character(dat2$Alpha))))
}
names(tol) <- unique(as.character(dat2$Alpha))
if(length(iwells) == 1) {
   iwells <- rep(iwells, length(unique(as.character(dat2$Alpha))))
}
if(length(jwells) == 1) {
   jwells <- rep(jwells, length(unique(as.character(dat2$Alpha))))
}
if(length(ipeaks) == 1) {
   ipeaks <- rep(ipeaks, length(unique(as.character(dat2$Alpha))))
}
names(iwells) <- unique(as.character(dat2$Alpha))
names(jwells) <- unique(as.character(dat2$Alpha))
names(ipeaks) <- unique(as.character(dat2$Alpha))

dg <- tapply(1:nrow(dat2), as.character(dat2$Alpha), function(i) {
   dat <- dat2[i, ]
   alpha <- as.character(dat2$Alpha[i[1]])
   tol <- tol[alpha]
   iwell <- iwells[alpha]
   jwell <- jwells[alpha]
   ipeak <- ipeaks[alpha]

   # find the minima/maxima
   myout <- get_wells_peaks(dat$Energy, tol=tol)

#   rl <- rle2(sign(diff(dat$Energy))) 
#   cis <- dat$Angle[rl$inds[1] + 1]
#   ts  <- dat$Angle[rl$inds[2] + 1]
#   trans <- dat$Angle[rl$inds[3] + 1]
   cis <- dat$Angle[myout$wells[iwell]]
   ts <- dat$Angle[myout$peaks[ipeak]]
   trans <- dat$Angle[myout$wells[jwell]]
   cat("Cis,   Ts,   Trans\n")
   cat(cis, ts, trans, "\n", sep=", ")
#   dat <- dat[dat$Run==88, ]
#   dg1 <- dat$Energy[dat$Angle==89.4] - dat$Energy[dat$Angle==179.4]
#   dg2 <- dat$Energy[dat$Angle==89.4] - dat$Energy[dat$Angle==-0.6]
#   dg3 <- dat$Energy[dat$Angle==-0.6] - dat$Energy[dat$Angle==179.4]
   dg1 <- dat$Energy[dat$Angle==ts] - dat$Energy[dat$Angle==trans]
   dg1err <- sqrt(dat$EnergySD[dat$Angle==ts]^2 + dat$EnergySD[dat$Angle==trans]^2)
   dg2 <- dat$Energy[dat$Angle==ts] - dat$Energy[dat$Angle==cis]
   dg2err <- sqrt(dat$EnergySD[dat$Angle==ts]^2 + dat$EnergySD[dat$Angle==cis]^2)
   dg3 <- dat$Energy[dat$Angle==cis] - dat$Energy[dat$Angle==trans]
   dg3err <- sqrt(dat$EnergySD[dat$Angle==cis]^2 + dat$EnergySD[dat$Angle==trans]^2)
   c(trans=dg1, cis=dg2, dg=dg3, errtrans=dg1err, errcis=dg2err, errdg=dg3err)
})
dg <- as.data.frame(do.call(rbind, dg))
#mdg <- reshape2::melt(dg)
#mdg <- cbind(Run=rep(1:86, 3), mdg)

#ggplot(mdg, aes(x=Run, y=value, col=variable)) + geom_line() + facet_wrap(~variable, ncol=1, scales="free_y")
write.table(dg, file="dg.txt")
