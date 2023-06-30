library(ggplot2)
library(bio3d)

#rootpath <- "resampling/sep1_10k_err"
rootpath <- "resampling/sep1_-30_210"
b_extrapolate <- FALSE
b_trash <- FALSE

pred_dash_line <- FALSE

# ffspr, ffpspr_pin1, ffpspr_ffpspr, ffpspr_tpp
tol=0

# ffpspr
#tol=1

# ffpspr_ppiase
#tol=3

ylim <- c(-5.5, 22)

## only valid when b_extrapolate=TRUE (to compare prediction and result)
b_remove_original <- FALSE
color_scales <- factor(format(c(0.20, 0.25, 0.30, 0.35, 0.40, 0.43, 0.45, 0.50, 1.00), nsmall=2),
                       levels=format(seq(0, 1, 0.01), nsmall=2))

datfile <- list.files(rootpath, pattern='^fep.*\\.out$', recursive=TRUE, full.names=TRUE)

#if(!b_trash) {
#   rm.inds <- grep("trash", datfile)
#   if(length(rm.inds)>0) {
#      datfile <- datfile[-rm.inds]
#   }
#}
if(b_trash) {
   rootpath2 <- file.path("trash", rootpath)
   datfile2 <- list.files(rootpath2, pattern='^fep.*\\.out$', recursive=TRUE, full.names=TRUE)
   if(length(datfile2)>0) {
      datfile <- c(datfile, datfile2)
   }
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

#datfile <- datfile[-1]

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

if(b_extrapolate && b_remove_original) {
   # remove Alpha == 1
   dat <- dat[which(dat$Alpha != "1.00"), ]
   if(length(unique(dat$Alpha))<2) {
      stop("Not enough Alpha's to expolate (must be >1)")
   }
}
if(b_extrapolate) {
   coefs <- tapply(1:nrow(dat), dat$Angle, function(i) {
   #   summary(lm(Energy~Alpha2, data=dat[i, ], weights=1.0/(dat$EnergySD[i]^2)))$coefficients
   #   summary(lm(Energy~Alpha2, data=dat[i, ]))$coefficients
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

   dat2$Energy <- dat2$Energy - dat2$Energy[dat2$Alpha=="1.00" & dat2$Angle==0.00]

   ## plot extrpolated
   p <- ggplot(dat2, aes(x=Angle, y=Energy, col=Alpha)) +
      geom_errorbar(aes(ymin=Energy-EnergySD, ymax=Energy+EnergySD), width=3, size=0.3)
   if(any(dat2$Alpha=="1.00")) {
      if(pred_dash_line) {
         p <- p + geom_line(aes(linetype=predict))
      } else {
         p <- p + geom_line()
      }
   } else {
      p <- p + geom_line()
   }
   p <- p +
      xlab("Omega Angle (degree)") +
      ylab("Free Energy (kcal/mol)") +
      scale_colour_discrete(drop=TRUE, limits=color_scales, breaks=unique(dat2$Alpha)) +
      ylim(ylim[1], ylim[2])
  
   if(b_trash) {  
      pdf(height=4, width=5, file="pmf_extrapolated_include_trash.pdf")
   } else {
      pdf(height=4, width=5, file="pmf_extrapolated.pdf")
   }
   print(p)
   dev.off()

} else {

   ## comment out if extrapolation is used
   dat2 <- dat
   dat.pred <- dat[dat$Alpha=="1.00", ]

   dat2$Energy <- dat2$Energy - dat2$Energy[dat2$Alpha=="1.00" & dat2$Angle==0.00]
   
   p <- ggplot(dat2, aes(x=Angle, y=Energy, col=Alpha)) + 
      geom_errorbar(aes(ymin=Energy-EnergySD, ymax=Energy+EnergySD), width=3, size=0.3) +
   #   geom_line(aes(linetype=predict)) +
      geom_line() +
      xlab("Omega Angle (degree)") +
      ylab("Free Energy (kcal/mol)") +
      scale_colour_discrete(drop=TRUE, limits=color_scales, breaks=unique(dat2$Alpha)) +
      ylim(ylim[1], ylim[2])
   
   if(b_trash) {  
      pdf(height=4, width=5, file="pmf_include_trash.pdf") 
   } else {
      pdf(height=4, width=5, file="pmf.pdf") 
   }
   print(p)
   dev.off()
}

save(dat2, file="dat2.RData")

source("funs.R")
out <- get_wells_peaks(dat.pred$Energy, tol=tol)
cat("Cis: ", dat.pred$Angle[out$wells[1]], "\n")
cat("TS: ", dat.pred$Angle[out$peaks[1]], "\n")
cat("Trans: ", dat.pred$Angle[out$wells[2]], "\n")

cat("DG_TS_cis: ", dat.pred$Energy[out$peaks[1]] - dat.pred$Energy[out$wells[1]], "\n")
cat("DG_TS_trans: ", dat.pred$Energy[out$peaks[1]] - dat.pred$Energy[out$wells[2]], "\n")
cat("DG_cis_trans: ", dat.pred$Energy[out$wells[1]] - dat.pred$Energy[out$wells[2]], "\n")
cat("DG_err_TS_cis: ", sqrt(dat.pred$EnergySD[out$peaks[1]]^2 + dat.pred$EnergySD[out$wells[1]]^2), "\n")
cat("DG_err_TS_trans: ", sqrt(dat.pred$EnergySD[out$peaks[1]]^2 + dat.pred$EnergySD[out$wells[2]]^2), "\n")
cat("DG_err_cis_trans: ", sqrt(dat.pred$EnergySD[out$wells[1]]^2 + dat.pred$EnergySD[out$wells[2]]^2), "\n")


#cat("DG_TS_cis: ", dat.pred$Energy[dat.pred$Angle==90] - dat.pred$Energy[dat.pred$Angle==0], "\n")
#cat("DG_TS_trans: ", dat.pred$Energy[dat.pred$Angle==90] - dat.pred$Energy[dat.pred$Angle==180], "\n")
#cat("DG_cis_trans: ", dat.pred$Energy[dat.pred$Angle==0] - dat.pred$Energy[dat.pred$Angle==180], "\n")
#cat("DG_err_TS_cis: ", sqrt(dat.pred$EnergySD[dat.pred$Angle==90]^2 + dat.pred$EnergySD[dat.pred$Angle==0]^2), "\n")
#cat("DG_err_TS_trans: ", sqrt(dat.pred$EnergySD[dat.pred$Angle==90]^2 + dat.pred$EnergySD[dat.pred$Angle==180]^2), "\n")
#cat("DG_err_cis_trans: ", sqrt(dat.pred$EnergySD[dat.pred$Angle==0]^2 + dat.pred$EnergySD[dat.pred$Angle==180]^2), "\n")

