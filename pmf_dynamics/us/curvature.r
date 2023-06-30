library(bio3d)

#rootpath <- "resampling/sep1_10k_err"
rootpath <- "resampling/sep1_-30_210"
b_extrapolate <- FALSE

## Set TRUE to include parameters under 'trash' for extrapolation 
## If TRUE, need to set trash path
b_extrapolate_include_trash <- FALSE
#rootpath_trash <- "trash/resampling/sep1_10k_err"
#rootpath_trash <- "trash/resampling/sep1_-30_210"
rootpath_trash <- file.path("trash", rootpath)

##
tol <- 0
#tol=c(4, 4, 4, 0)

## angle range around the min(ax)imum to draw the curve
offset <- 5  # trans
offset2 <- NULL  # TS
#offset2 <- c(3, 5, 5, 5)

## pick up where to draw the curve (No. wells/peaks identified)
iwells <- 2
ipeaks <- 1
#ipeaks <- c(2, 1, 1, 1)
###

if(is.null(offset2)) {
   offset2 <- offset
}

datfile <- list.files(rootpath, pattern='^fep\\.[0-9]*\\.out$', recursive=TRUE, full.names=TRUE)

if( (!b_extrapolate) & b_extrapolate_include_trash ) {
   b_extrapolate_include_trash <- FALSE
}
if(b_extrapolate_include_trash) {
   if(!is.null(rootpath_trash)) {
      datfile_trash <- list.files(rootpath_trash, pattern='^fep\\.[0-9]*\\.out$', recursive=TRUE, full.names=TRUE)
   } else {
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

#ggplot(dat, aes(x=Angle, y=Energy, col=Alpha)) + geom_line()

#ggplot(dat, aes(x=1-as.numeric(as.character(Alpha)), y=Energy)) + 
#   geom_point() + 
#   geom_smooth(method="lm", se=FALSE, size=0.5) + 
#   facet_wrap(~Angle)

if(b_extrapolate) {
   dat <- dat[which(dat$Alpha!="1.00"), ]
   if(length(unique(dat$Alpha))<2) {
      stop("Not enough Alpha's to expolate (must be >1)")
   }
   coefs <- tapply(1:nrow(dat), dat$Angle, function(i) {
   #   summary(lm(Energy~Alpha2, data=dat[i, ], weights=1.0/(dat$EnergySD[i]^2)))$coefficients
   #   summary(lm(Energy~Alpha2, data=dat[i, ]))$coefficients
       inds <- which(dat$Alpha[i]!="1.00")
       if(length(inds)>0) {
          if(any(!is.finite(dat$EnergySD[i][inds]))) {
             summary(lm(Energy~Alpha2, data=dat[i, ][inds, ]))$coefficients
          } else {
             summary(lm(Energy~Alpha2, data=dat[i, ][inds, ], weights=1/(dat$EnergySD[i][inds]^2)))$coefficients
          }
       } else {
          NA
       }
   })

   dat.pred <- data.frame(Angle=as.numeric(names(coefs)),
      Energy=unname(sapply(coefs, "[", "(Intercept)", 1)),
      EnergySD=unname(sapply(coefs, "[", "(Intercept)", 2)),
      Alpha="1.00", Alpha2=0.0, predict=TRUE)

   dat2 <- cbind(dat, predict=FALSE)
   dat2 <- rbind(dat2, dat.pred)
} else {
   dat2 <- dat
   if("1.00" %in% dat$Alpha) {
      dat.pred <- dat[dat$Alpha=="1.00", ]
   } else {
      dat.pred <- NULL
   }
}

if(b_extrapolate_include_trash) {
   ## remove trash in plotting
   alphas_trash <- format(round(as.numeric(basename(dirname(datfile_trash)))/14.0, 2), nsmall=2)
   dat2 <- dat2[!dat2$Alpha %in% alphas_trash, ]
}


#inds <- dat2$Alpha==0.25
#plot(dat2$Angle[inds]+40, dat2$Energy[inds]+1, typ="l")
#mydat <- dat2[inds, ]
##dat <- mydat[mydat$Angle >=60 & mydat$Angle <=100, ]
#dat <- mydat[mydat$Angle >=-25 & mydat$Angle <=25, ]
#lines(dat$Angle, dat$Energy, col="blue", lw=2)
##fit <- lm(log(-Energy+10) ~ 1 + offset(2*log(Angle+60)), data=dat)
#fit <- lm(log(Energy+1) ~ 1 + offset(2*log(Angle+40)), data=dat)
#exp(coef(fit))
#sqrt(exp(coef(fit)))

if(FALSE) {

inds <- dat2$Alpha=="0.25"
plot(dat2$Angle[inds], dat2$Energy[inds], typ="l")
mydat <- dat2[inds, ]
#dat <- mydat[mydat$Angle >=45 & mydat$Angle <=115, ]
dat <- mydat[mydat$Angle >=120 & mydat$Angle <=210, ]
lines(dat$Angle, dat$Energy, col="blue", lw=2)
#fit <- lm(log(-Energy+10) ~ 1 + offset(2*log(Angle+60)), data=dat)
#fit <- nls(Energy ~ A*Angle^2+B, data=dat)
fit <- nls(Energy ~ AA*Angle^2 + BB*Angle + CC, data=dat)
lines(dat$Angle, predict(fit), col='red', lwd=3)

omega <- sqrt(abs(coef(fit)[1])) * sqrt(2)
omega

}

if(!dir.exists("curvature")) {
   dir.create("curvature")
}

source("funs.R")
if(length(tol) == 1) {
   tol <- rep(tol, length(unique(as.character(dat2$Alpha))))
}
names(tol) <- unique(as.character(dat2$Alpha))
if(length(offset) == 1) {
   offset <- rep(offset, length(unique(as.character(dat2$Alpha))))
}
if(length(offset2) == 1) {
   offset2 <- rep(offset2, length(unique(as.character(dat2$Alpha))))
}
names(offset) <- unique(as.character(dat2$Alpha))
names(offset2) <- unique(as.character(dat2$Alpha))
if(length(iwells) == 1) {
   iwells <- rep(iwells, length(unique(as.character(dat2$Alpha))))
}
if(length(ipeaks) == 1) {
   ipeaks <- rep(ipeaks, length(unique(as.character(dat2$Alpha))))
}
names(iwells) <- unique(as.character(dat2$Alpha))
names(ipeaks) <- unique(as.character(dat2$Alpha))
out <- tapply(1:nrow(dat2), as.character(dat2$Alpha), function(i) {

#   cat(as.character(dat2$Alpha[i[1]]), "...\n", sep="")  
   inds <- i 
   alpha <- as.character(dat2$Alpha[inds[1]])
   mydat <- dat2[inds, ]
   tol <- tol[alpha]
   offset <- offset[alpha]
   offset2 <- offset2[alpha]
   iwell <- iwells[alpha]
   ipeak <- ipeaks[alpha]

   myout <- get_wells_peaks(mydat$Energy, tol=tol)

   abound1 <- mydat$Angle[myout$wells[iwell] + offset]
   val <- mydat$Energy[mydat$Angle==abound1]
   tind <- which.min(abs(mydat$Energy - val)[myout$peaks[ipeak]:myout$wells[iwell]])
#   tinds <- tinds & mydat$Angle < amin
   abound2 <- mydat$Angle[myout$peaks[ipeak]:myout$wells[iwell]][tind]
 
   dat <- mydat[mydat$Angle >=abound2 & mydat$Angle <=abound1, ]
   pdf(paste("curvature/w0_", alpha, ".pdf", sep=""), height=4, width=4) 
   plot(dat2$Angle[inds], dat2$Energy[inds], typ="l")
   lines(dat$Angle, dat$Energy, col="blue", lw=2)
   suppressWarnings( fit <- nls(Energy ~ AA*Angle^2 + BB*Angle + CC, data=dat) )
   lines(dat$Angle, predict(fit), col='red', lwd=3)
   dev.off()
 
   omega <- unname(sqrt(abs(coef(fit)[1])) * sqrt(2))

   
   abound2 <- mydat$Angle[myout$peaks[ipeak] - offset2]
   val <- mydat$Energy[mydat$Angle==abound2]
   tind <- which.min(abs(mydat$Energy - val)[myout$peaks[ipeak]:myout$wells[iwell]])
#   tinds <- tinds & mydat$Angle < amin
   abound1 <- mydat$Angle[myout$peaks[ipeak]:myout$wells[iwell]][tind]
#   tinds <- abs(mydat$Energy - val)<=0.001
#   tinds <- tinds & mydat$Angle < amax
#   abound2 <- which(tinds)[sum(tinds)]

   dat <- mydat[mydat$Angle >=abound2 & mydat$Angle <=abound1, ]
   pdf(paste("curvature/wb_", alpha, ".pdf", sep=""), height=4, width=4) 
   plot(dat2$Angle[inds], dat2$Energy[inds], typ="l")
   lines(dat$Angle, dat$Energy, col="blue", lw=2)
   suppressWarnings( fit <- nls(Energy ~ AA*Angle^2 + BB*Angle + CC, data=dat) )
   lines(dat$Angle, predict(fit), col='red', lwd=3)
   dev.off()
 
   omega <- c(w0=omega, wb=unname(sqrt(abs(coef(fit)[1])) * sqrt(2)))
   omega
})
out <- do.call("cbind", out)
out <- round(out, 3)

write.table(out, file="curvature.txt")
