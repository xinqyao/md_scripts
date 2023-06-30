library(ggplot2)

#rootpath <- "resampling/sep1_10k_err"
rootpath <- "resampling/sep1_-30_210"
b_extrapolate <- FALSE

## Set TRUE to include parameters under 'trash' for extrapolation 
## If TRUE, need to set trash path
b_extrapolate_include_trash <- FALSE
#rootpath_trash <- "trash/resampling/sep1_10k_err"
#rootpath_trash <- "trash/resampling/sep1_-30_210"
rootpath_trash <- file.path("trash", rootpath)

## obsolete; run the script under 'trash' directly to check trash
b_trash <- FALSE

datfile <- list.files(rootpath, pattern='^fep.*\\.out$', recursive=TRUE, full.names=TRUE)

if(!b_trash) {
   rm.inds <- grep("trash", datfile)
   if(length(rm.inds)>0) {
      datfile <- datfile[-rm.inds]
   }
}

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

#datfile <- datfile[grep("wham", datfile)]
#datfile <- datfile[order(as.numeric(sub("fep.(.*).out$", "\\1", basename(datfile))), na.last=FALSE)]
#datfile <- datfile[-grep('bak', datfile)]
#datfile <- datfile[-grep('88', datfile)]
#datfile <- datfile[-grep('fep.45.out', datfile)]
#datfile <- datfile[1:length(datfile) %in% c(grep("4.2", datfile), grep("4.9", datfile), grep("5.6", datfile), grep("6.3", datfile))]
#datfile <- datfile[-c(61:63)]
#datfile <- datfile[-grep("fep.21.out", datfile)]

dat <- NULL
for(i in datfile) {
   tdat <- read.table(i)
##   fref <- sum(tdat[tdat[, 1]%in% c(178.75, 181.25), 2])/2 # reference to trans
   fref <- tdat[tdat[, 1]%in% 180, 2] # reference to trans
   fref2 <- tdat[tdat[, 1]%in% 0, 2] # reference to cis 
##   tdat <- cbind(tdat, as.numeric(sub("\\.\\/(.*)\\/wham.*", "\\1", i))/14.0)
##   tdat <- cbind(tdat, as.numeric(basename(dirname(i)))/14.0)
   tdat <- cbind(tdat, as.numeric(sub(".*\\/([\\.0-9]+)\\/fep.*", "\\1", i))/14.0)
   runs <- as.numeric(sub("fep\\.(.*)\\.out$", "\\1", basename(i)))
   runs[is.na(runs)] <- 1
   tdat <- cbind(tdat, runs)
   tdat[, 2] <- tdat[, 2] - (fref + fref2)/2
   dat <- rbind(dat, tdat)
}
colnames(dat) <- c("Angle", "Energy", "EnergySD", "Alpha", "Run")
dat$Alpha <- factor(format(round(dat$Alpha, 2), nsmall=2), levels=format(seq(0, 1, 0.01), nsmall=2))
dat$Alpha2 <- 1-as.numeric(as.character(dat$Alpha))

if(b_extrapolate) {
   # remove Alpha == 1
   dat <- dat[which(dat$Alpha != "1.00"), ]
   if(length(unique(dat$Alpha))<2) {
      stop("Not enough Alpha's to expolate (must be >1)")
   }

   coefs <- tapply(1:nrow(dat), dat$Run, function(j) {
      tapply(1:nrow(dat[j, ]), dat$Angle[j], function(i) {
   #      cat(dat$Run[j[1]], "; ", dat$Angle[j][i[1]], "\n")
   #   coef(lm(Energy~Alpha2, data=dat[j[i], ]))
   #      summary(lm(Energy~Alpha2, data=dat[j[i], ]))$coefficients
         inds <- which(dat$Alpha[j[i]] != "1.00")
         if(length(inds)>0) {
            if(any(!is.finite(dat$EnergySD[j[i]][inds]))) {
               summary(lm(Energy~Alpha2, data=dat[j[i], ][inds, ]))$coefficients
            }
            else {
               summary(lm(Energy~Alpha2, data=dat[j[i], ][inds, ], weights=1/(dat$EnergySD[j[i]][inds]^2)))$coefficients
            }
         }
         else {
            NA
         }
   }) })
   maxruns <- tapply(1:nrow(dat), dat$Alpha, function(i) max(dat$Run[i]))
   rmin <- min(maxruns[-length(maxruns)], na.rm=TRUE)
   rmin.inds <- which(maxruns==rmin)
   coefs <- coefs[as.numeric(names(coefs))<=rmin]
   
   dat.pred <- data.frame(Angle=as.numeric(sapply(coefs, names)),
      Energy=unname(as.numeric(sapply(coefs, sapply, "[", "(Intercept)", 1))),
      EnergySD=unname(as.numeric(sapply(coefs, sapply, "[", "(Intercept)", 2))),
      Alpha="1.00", Run=subset(dat, Alpha==names(rmin.inds)[1])$Run, Alpha2=0.0, predict=TRUE)
#      Alpha="1.00", Run=subset(dat, Alpha==dat$Alpha[1])$Run, Alpha2=0.0, predict=TRUE)
   
   dat2 <- cbind(dat, predict=FALSE)
   dat2 <- rbind(dat2, dat.pred)

} else {

   dat2 <- dat
   if("1.00" %in% dat$Alpha) {
      dat.pred <- dat[dat$Alpha=="1.00", ]
   }
   else {
      dat.pred <- NULL
   }
}

if(b_extrapolate_include_trash) {
   ## remove trash in plotting
   alphas_trash <- format(round(as.numeric(basename(dirname(datfile_trash)))/14.0, 2), nsmall=2)
   dat2 <- dat2[!dat2$Alpha %in% alphas_trash, ]
}
 
p <- ggplot(dat2, aes(x=Angle, y=Energy, group=Run, col=Run)) + geom_line() + facet_wrap(~Alpha, scales="free_y")

pdf(file="convergence.pdf", height=4, width=6)
print(p)
dev.off()

rmsd <- tapply(1:nrow(dat2), as.character(dat2$Alpha), function(i) {
   rmsd <- NULL
   for(j in 2:max(dat2$Run[i])) {
      rmsd <- c(rmsd, sqrt(mean((dat2$Energy[i][dat2$Run[i]==j] - dat2$Energy[i][dat2$Run[i]==(j-1)])^2)))
   }
   rmsd
})
#as.data.frame(do.call(cbind, rmsd))
rmsd <- data.frame(variable=rep(names(rmsd), sapply(rmsd, length)), 
                   value = unlist(rmsd))
#rmsd <- reshape2::melt(rmsd)
rmsd <- cbind(rmsd, Run=unlist(lapply(rle(as.character(rmsd$variable))$lengths, function(x)2:(x+1))))
p <- ggplot(rmsd, aes(x=Run, y=value)) + geom_line() + ylim(0, NA) + facet_wrap(~variable, scales="free_y")
pdf(file="convergence_rmsd.pdf", height=4, width=6)
print(p)
dev.off()

mdgs <- tapply(1:nrow(dat2), as.character(dat2$Alpha), function(j) {
   mydat <- dat2[j, ]
   dg <- tapply(1:nrow(mydat), mydat$Run, function(i) {
      dat <- mydat[i, ]
      dg1 <- dat$Energy[dat$Angle==90] - dat$Energy[dat$Angle==180]
      dg1sd <- sqrt(dat$EnergySD[dat$Angle==90]^2 + dat$EnergySD[dat$Angle==180]^2)
      dg2 <- dat$Energy[dat$Angle==90] - dat$Energy[dat$Angle==0]
      dg2sd <- sqrt(dat$EnergySD[dat$Angle==90]^2 + dat$EnergySD[dat$Angle==0]^2)
      dg3 <- dat$Energy[dat$Angle==0] - dat$Energy[dat$Angle==180]
      dg3sd <- sqrt(dat$EnergySD[dat$Angle==0]^2 + dat$EnergySD[dat$Angle==180]^2)
      c(trans=dg1, trans_sd=dg1sd, cis=dg2, cis_sd=dg2sd, dg=dg3, dg_sd=dg3sd)
   })
   dg <- as.data.frame(do.call(rbind, dg))
   mdg <- reshape2::melt(dg[, c("trans", "cis", "dg")])
   mdg <- cbind(Run=rep(1:max(mydat$Run), 3), mdg)
   mdgsd <- reshape2::melt(dg[, c("trans_sd", "cis_sd", "dg_sd")])
   #mdgsd <- cbind(Run=rep(1:max(dat$Run), 3), mdg_sd)
   mdg <- cbind(mdg, sd=mdgsd$value, Alpha=mydat$Alpha[1])
})
mdgs <- do.call(rbind, mdgs)
p <- ggplot(mdgs, aes(x=Run, y=value, col=variable)) + 
   geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=1, size=0.3) +
   geom_line() + 
   facet_wrap(~Alpha, ncol=3, scales="free_y")

pdf(file="convergence_dg_all.pdf", height=4, width=6)
print(p)
dev.off()

if(!is.null(dat.pred)) {
   dg <- tapply(1:nrow(dat.pred), dat.pred$Run, function(i) {
      dat <- dat.pred[i, ]
      dg1 <- dat$Energy[dat$Angle==90] - dat$Energy[dat$Angle==180]
      dg1sd <- sqrt(dat$EnergySD[dat$Angle==90]^2 + dat$EnergySD[dat$Angle==180]^2)
      dg2 <- dat$Energy[dat$Angle==90] - dat$Energy[dat$Angle==0]
      dg2sd <- sqrt(dat$EnergySD[dat$Angle==90]^2 + dat$EnergySD[dat$Angle==0]^2)
      dg3 <- dat$Energy[dat$Angle==0] - dat$Energy[dat$Angle==180]
      dg3sd <- sqrt(dat$EnergySD[dat$Angle==0]^2 + dat$EnergySD[dat$Angle==180]^2)
      c(trans=dg1, trans_sd=dg1sd, cis=dg2, cis_sd=dg2sd, dg=dg3, dg_sd=dg3sd)
   })
   dg <- as.data.frame(do.call(rbind, dg))
   mdg <- reshape2::melt(dg[, c("trans", "cis", "dg")])
   mdg <- cbind(Run=rep(1:max(dat.pred$Run), 3), mdg)
   mdgsd <- reshape2::melt(dg[, c("trans_sd", "cis_sd", "dg_sd")])
   #mdgsd <- cbind(Run=rep(1:max(dat$Run), 3), mdg_sd)
   mdg <- cbind(mdg, sd=mdgsd$value)

   p <- ggplot(mdg, aes(x=Run, y=value, col=variable)) + 
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=1, size=0.3) +
      geom_line() + 
      facet_wrap(~variable, ncol=1, scales="free_y")
   
   pdf(file="convergence_dg.pdf", height=4, width=6)
   print(p)
   dev.off()
   
   write.table(mdg, file="mdg.txt")

}

# write dg time-series separately; dg/
if(!dir.exists("dg")) {
   dir.create("dg")
}

tapply(1:nrow(dat2), as.character(dat2$Alpha), function(inds) {
   dat <- dat2[inds, ]

   dg <- tapply(1:nrow(dat), dat$Run, function(i) {
      dat <- dat[i, ]
      dg1 <- dat$Energy[dat$Angle==90] - dat$Energy[dat$Angle==180]
      dg1sd <- sqrt(dat$EnergySD[dat$Angle==90]^2 + dat$EnergySD[dat$Angle==180]^2)
      dg2 <- dat$Energy[dat$Angle==90] - dat$Energy[dat$Angle==0]
      dg2sd <- sqrt(dat$EnergySD[dat$Angle==90]^2 + dat$EnergySD[dat$Angle==0]^2)
      dg3 <- dat$Energy[dat$Angle==0] - dat$Energy[dat$Angle==180]
      dg3sd <- sqrt(dat$EnergySD[dat$Angle==0]^2 + dat$EnergySD[dat$Angle==180]^2)
      c(trans=dg1, trans_sd=dg1sd, cis=dg2, cis_sd=dg2sd, dg=dg3, dg_sd=dg3sd)
   })
   dg <- as.data.frame(do.call(rbind, dg))
   mdg <- reshape2::melt(dg[, c("trans", "cis", "dg")])
   mdg <- cbind(Run=rep(1:max(dat$Run), 3), mdg)
   mdgsd <- reshape2::melt(dg[, c("trans_sd", "cis_sd", "dg_sd")])
   #mdgsd <- cbind(Run=rep(1:max(dat$Run), 3), mdg_sd)
   mdg <- cbind(mdg, sd=mdgsd$value)
   
   p <- ggplot(mdg, aes(x=Run, y=value, col=variable)) + 
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=1, size=0.3) +
      geom_line() + 
      guides(color="none") +
      facet_wrap(~variable, ncol=1, scales="free_y")
   
   pdf(file=file.path("dg", paste("convergence_dg_", dat$Alpha[1], ".pdf", sep="")), height=4, width=4)
   print(p)
   dev.off()
})
