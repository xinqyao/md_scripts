library(ggplot2)
#library(deming)

## use the following formula for error estimate (see wiki)
## sf = sqrt((af/ax)^2*sx^2 + (af/ay)^2*sy^2 + ...)

ww <- read.table("../us/curvature.txt", header=TRUE)
dg <- read.table("../us/dg.txt", header=TRUE)
kk <- unlist(read.table("rate.txt"))
kk_err <- unlist(read.table("rate_err.txt"))
kb = 0.0019872041
temp = 300

xx <- dg$trans / (kb*temp)
xx_err <- dg$errtrans / (kb*temp)

dat <- data.frame(dg=xx[-length(xx)], logk=log(unname(kk)/unname((as.numeric(ww[1, -ncol(ww)]*ww[2, -ncol(ww)])))), 
                  dg_err=xx_err[-length(xx)], logk_err=kk_err/kk)

#fit <- lm(logk ~ 1 + offset(-dg), data=dat)
#fit <- deming(logk ~ 1 + offset(-dg), data=dat, xstd=dg_err, ystd=logk_err) 
# deming function does not support offset model fitting

#fit <- lm(logk ~ 1 + offset(-dg), data=dat, weights=1/(logk_err^2+dg_err^2))
#fit <- lm(dg ~ 1 + offset(-logk), data=dat, weights=1/dg_err^2)
fit <- lm(logk ~ 1 + offset(-dg), data=dat, weights=1/logk_err^2)

pred <- -xx + coef(fit)
#pred_err <- sqrt(xx_err^2+summary(fit)$coefficients[1, 'Std. Error']^2)
pred_err <- rep(summary(fit)$coefficients[1, 'Std. Error'], length(pred))

#exp(pred) * ww$X1.0[1] * ww$X1.0[2] * 10^9 * 10^3  #unmodified rate estimate; Unit: 10^3/s

## write out data for combined plot later
write.table(rbind(xx, pred, pred_err, cbind(t(dat$logk), NA)), col.names=FALSE, row.names=FALSE, file="plot_kinetics.txt")

## rate and diffusion coefficient 
out <- c(kk, exp(pred[length(pred)])*ww[1, ncol(ww)]*ww[2, ncol(ww)])   # rate unit: 1/ns
out <- rbind(out, c(kk_err, exp(pred[length(pred)])*ww[1, ncol(ww)]*ww[2, ncol(ww)]*pred_err[length(pred)]))   # rate err
out <- rbind(out, 1/out[1, ])   # time: ns
out <- rbind(out, out[2, ]/(out[1, ]^2))   # time err
out <- rbind(out, exp(coef(fit)) * 2 * pi * kb * temp * 10^9)  # diffusion unit: deg^2/s
out <- rbind(out, exp(coef(fit)) * 2 * pi * kb * temp * 10^9 * summary(fit)$coefficients[1, 'Std. Error']) # Std. Error of diffusion coef.
colnames(out) <- rownames(dg)
rownames(out) <- c("k_(1/ns)", "k_Std._Error", "tau_(ns)", "tau_Std._Error", "Deff_(deg^2/s)", "Deff_Std._Error")
write.table(format(out, scientific=TRUE), quote=FALSE, file="kinetics.txt")

preddat <- data.frame(dg=xx[length(xx)], logk=pred[length(pred)], dg_err=xx_err[length(xx)], logk_err=pred_err[length(pred)])
p <- ggplot(data=dat, aes(x=dg, y=logk)) + geom_point(shape=19) +
#   geom_smooth(method="lm", formula=y~1+offset(-x), se=FALSE) +
#   geom_errorbar(aes(ymin=logk-logk_err, ymax=logk+logk_err), width=.1, size=.3) +
#   geom_errorbar(aes(xmin=dg-dg_err, xmax=dg+dg_err), width=.1, size=.3) +
   geom_line(aes(x,y), data=data.frame(x=xx, y=pred), col="blue") +
   annotate("point", x=xx[length(xx)], y=pred[length(pred)], shape=1) +
#   geom_errorbar(data=preddat, aes(ymin=logk-logk_err, ymax=logk+logk_err), width=.1, size=.3, linetype=3) +
#   geom_errorbar(data=preddat, aes(xmin=dg-dg_err, xmax=dg+dg_err), width=.1, size=.3) +
   xlab("DG/kBT") + ylab("ln(k/w0wk)") #+

pdf(file="kinetics.pdf", width=5, height=4)
print(p)
dev.off()

outdat <- list(dat=dat, xx=xx, xx_err=xx_err, pred=pred, pred_err=pred_err)
save(outdat, file="dat.RData")

