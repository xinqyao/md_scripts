library(ggplot2)
library(pracma)
source("funs.R")

#mode <- "bin"
#mode <- "nonbin"
mode <- "km"

dir <- "ts_time_life1"
#paras <- c('3.5', '4.2', '4.9', '5.6', '6.3', '7.0')
paras <- c('4.9', '5.6', '6.3')
root <- "."

# Define the type of fitting:
# "single" for single exponential (linear fitting in the logrithm space)
# "double" for double exponentials
#fit <- "single"   
#fit <- "power"   
fit <- c("multiple", "multiple", "multiple")
#bb <- c(-1.0, -0.1, -0.01) ## exponentials
bbs <- list("4.9"=c(-0.1, -0.01),
           "5.6"=c(-0.1, -0.01),
           "6.3"=c(-0.1, -0.01)) ## exponentials


# power order
aa <- 1.0001

if(length(fit) == 1) {
   fit <- rep(fit, length(paras))
}
names(fit) <- paras

if(length(aa) == 1) {
   aa <- rep(aa, length(paras))
}
names(aa) <- paras

color_scales <- factor(format(c(0.20, 0.25, 0.30, 0.35, 0.40, 0.43, 0.45, 0.50, 1.00), nsmall=2),
                       levels=format(seq(0, 1, 0.01), nsmall=2))

alphas <- format(round(as.numeric(paras) / 14.0, 2), nsmall=2)
names(alphas) <- paras
kk <- NULL
kks <- NULL
err <- NULL
mydat.all <- NULL
ssq <- NULL  # residual sum of squares
coefs <- NULL  # coefficients
if(mode == "bin") {
   if(file.exists("random.RData")) {
      rinds_loaded <- TRUE
      load("random.RData")
   } else {
      rinds <- NULL
      rinds_loaded <- FALSE
   }
}

for(para in paras) {

#   nsim <- length(list.files(file.path(root, para), "ang.[0-9]*.out"))
#   nsim <- 100

   ## simpler
   dat <- read.table(paste(dir, "/ts_time_", para, ".txt", sep=""))
   dat0 <- dat
   nsim0 <- nrow(dat)
   dat <- dat[dat[, 3] != -1, ]
   dat <- dat[, -1]
   nsim <- nrow(dat)
   ntrash <- sum(dat[, 2]==0)
   dat <- dat[, 1][dat[, 2]==1]
   dat <- sort(dat)
   nn <- nsim - ntrash
   nn0 <- nsim0 - ntrash
   cat("Para: ", para, ", #Sim: ", nsim, "#EffSim: ", nn, "\n")
   if(mode == "bin") {
      mytimes <- (dat-1)*10*0.002*0.001  # sampled times: ns
#      nmax <- floor(nn0 * 0.9)
#      mytimes <- mytimes[1:ifelse(nmax>length(mytimes), length(mytimes), nmax)] # only take up to 10% survival (ignore small sample size cases)
      xx <- seq(min(mytimes), max(mytimes), length.out=25)
     
      if(!rinds_loaded) {
         inds <- lapply(1:1000, function(i) {
            round(runif(length(mytimes), 1, length(mytimes)))
         })
      } else {
         inds <- rinds[[para]]
      }
      yy <- sapply(inds, function(i) {
         mytimes <- mytimes[i]
         cats <- cut(mytimes, breaks=c(xx[1]-.001, xx))
         y_1 <- cumsum( table(cats) )
         y <- 1 - y_1/nn0
      })
#      yy0 <- (nn - sapply(xx, function(t) sum(mytimes <= t)))/nn
      cats <- cut(mytimes, breaks=c(xx[1]-.001, xx))
      y_1 <- cumsum( table(cats) )
      yy0 <- 1 - y_1/nn0
      mydat <- data.frame(time=xx, prob=yy0, SE=apply(log(yy), 1, sd)) # ns; SE in log scale
   } else if(mode == "nonbin") {
      mydat <- data.frame(time=(dat-1)*10*0.002*0.001, prob=(nn0-c(1:length(dat)))/nn0) # ns
   } else if(mode == "km") {
      dat0[, 2] <- dat0[, 2] - 1
      mytimes <- dat0[dat0[, 3]==1, 2]
      unit <- 10*0.002*0.001 # ns
#      mydat <- get_survival_km(dat0, unit=unit, t=seq(min(mytimes), max(mytimes), length.out=25))
      mydat <- get_survival_km(dat0, unit=unit, t=seq(0, max(mytimes), 10/unit))
#      mydat <- get_survival_km(dat0, unit=unit)
      mydat$SE <- mydat$lse_max - mydat$lse_min
   }
   mydat <- mydat[mydat$prob>0, ]
#   mydat <- mydat[-nrow(mydat), ]  ## always remove the last one (otherwise log(0))

   #ggplot(mydat, aes(x=time, y=log(prob))) + geom_point() + 
   ##   xlim(c(0, 10)) + ylim(c(-7, 0)) +
   #   geom_smooth(method="lm", formula=y~x+0, se=FALSE)

   myfit <- fit[para]
   a <- aa[para]
   if(!is.null(mydat$SE)) {
      ww <- 1/(mydat$SE)^2
      ww[is.na(ww) | !is.finite(ww)] <- 0
   } else {
      ww <- rep(1, nrow(mydat))
   }
   if(myfit == "single") {
#         out <- lm(log(prob)~time+0, data=mydat, weights=(c(1/(mydat$SE[-nrow(mydat)])^2, 0))) ## remove the last point for fitting (0 SE)
      out <- lm(log(prob)~time+0, data=mydat, weights=ww)
#      out <- lm(log(prob)~time+0, data=mydat)
      #summary(out)
      kk <- c(kk, abs(coef(out)))
      err <- c(err, summary(out)$coefficients[1, "Std. Error"])
      ssq <- c(ssq, summary(out)$r.squared)
      coefs <- c(coefs, list(1.0))
      kks <- c(kks, list(abs(coef(out))))
      pred <- exp(predict(out))
   } else if(myfit == "multiple") {
      initval <- bbs[[para]]
      out <- mexpfit(mydat$time, mydat$prob, p0=initval, w=ww, const=FALSE)
     
      Ks <- abs(out$b)
      As <- abs(out$a)
      kk <- c(kk, sum(Ks*As)/sum(As))
      err <- c(err, 0)
      ssq <- c(ssq, out$ssq)
      coefs <- c(coefs, list(As))
      kks <- c(kks, list(Ks))
      pred <- 0
      for(i in 1:length(Ks)) {
         pred <- pred + As[i]*exp(-Ks[i]*mydat$time)
      }
   } else if(myfit == "power") {
      initval <- list(k1=0.001)
      out <- nls(prob~((a-1)*k1*time+1)^(-1/(a-1)), data=mydat, start=initval, weights=ww)
#      out <- nls(prob~(k1*time+1)^(-1), data=mydat)
      kk <- c(kk, abs(coef(out)["k1"]))
      ssq <- c(ssq, sum(summary(out)$residuals^2))
      pred <- predict(out)
   }
   mydat <- cbind(mydat, pred=pred)
   mydat.all <- rbind(mydat.all, 
        cbind(mydat, alpha=factor(alphas[para], levels=format(seq(0, 1, 0.01), nsmall=2))))
   if(mode == "bin" && !rinds_loaded) {
      rinds <- c(rinds, list(inds))
      names(rinds)[length(rinds)] <- para
   }
}

write.table(kk, file="rate.txt", row.names=FALSE, col.names=FALSE)
write.table(err, file="rate_err.txt", row.names=FALSE, col.names=FALSE)
write.table(ssq, file="rss.txt", row.names=FALSE, col.names=FALSE)
lmax <- max(sapply(kks, length))
tmat <- matrix(0, nrow=length(kks), ncol=lmax)
for(i in 1:length(kks)) tmat[i, 1:length(kks[[i]])] <- kks[[i]]
write.table(tmat, file="kks.txt", row.names=FALSE, col.names=FALSE)
for(i in 1:length(coefs)) tmat[i, 1:length(coefs[[i]])] <- coefs[[i]]
write.table(tmat, file="coefs.txt", row.names=FALSE, col.names=FALSE)

if(mode == "bin" && !rinds_loaded) {
   save(rinds, file="random.RData")
}

### fake data just to make color consistent
#mydat.all <- rbind(data.frame(time=c(20, 40), prob=c(0.2, 0.1), alpha=c(2.8, 2.8)), mydat.all)
#mydat.all <- rbind(mydat.all, data.frame(time=c(20, 40), prob=c(0.2, 0.1), alpha=c('5.6', '5.6')))
#mydat.all <- rbind(mydat.all, data.frame(time=c(20, 40), prob=c(0.2, 0.1), alpha=c('14.0', '14.0')))

p <- ggplot(mydat.all, aes(x=time, y=log(prob), color=alpha, group=alpha)) +
#   xlim(0, 90) + ylim(-4, 0) +
#   geom_smooth(method="lm", formula=y~x+0, se=FALSE) +
#   geom_line(aes(x=time, y=log(pred), color=alpha, group=alpha)) +
   geom_point()

if(!is.null(mydat.all$lse_min)) {
   p <- p + geom_ribbon(aes(ymin=lse_min, ymax=lse_max, fill=alpha), color=NA, alpha=.3, show.legend=FALSE)
} else if(!is.null(mydat.all$SE)) {
   p <- p + geom_errorbar(aes(ymin=log(prob)-SE, ymax=log(prob)+SE), width=1, size=.3)
}
p <- p +
   geom_line(aes(x=time, y=log(pred)), color="gray30") +
   xlab("Time (ns)") + 
   ylab("ln(S(t))") +
   scale_colour_discrete(drop=TRUE, limits=color_scales, breaks=unique(mydat.all$alpha))+
   scale_fill_discrete(drop=TRUE, limits=color_scales, breaks=unique(mydat.all$alpha))
#   scale_y_continuous(trans="log")
#   scale_y_continuous(trans="log", breaks=scales::trans_breaks("log", function(x) exp(x)))


pdf(file="survival.pdf", width=5, height=4)
print(p)
dev.off()

