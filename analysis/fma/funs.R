## partial least square ##
## Ref: Krivobokova, et al., Biophys J (2012) 103:786-796
pls <- function(X, f, k=1) {
   require(MASS)
   # center data
   Xmean <- colMeans(X)
   X <- sweep(X, 2, Xmean)
   f <- f -  mean(f)
   Wk <- matrix( t(X) %*% f )
   Tk <- matrix( X %*% Wk )
   if(k > 1) {
      for (i in 2:k) {
#         Tk_1 <- Tk[, 1:(i-1), drop=FALSE]
#         tWk <- t(X) %*% (f - Tk_1 %*% ginv(t(Tk_1) %*% Tk_1) %*% t(Tk_1) %*% f)
         tWk <- t(X) %*% (f - Tk %*% ginv(t(Tk) %*% Tk) %*% t(Tk) %*% f)
         Wk <- cbind(Wk, tWk)
         Tk <- cbind(Tk, X %*% tWk)
      }
   }
   return(list(Wk=Wk, Tk=Tk, Xmean=Xmean)) 
}

## for test purpose (normalized every step following Kramer et al., J Chemo Intel Lab Sys 2008, 94:60-69.
.pls <- function(X, f, k=1) {
   require(MASS)
   # center data
   Xmean <- colMeans(X)
   X <- sweep(X, 2, Xmean)
   X0 <- X
   f <- f -  mean(f)
   Wk <- NULL; Tk <- NULL 
   for (i in 1:k) {
      tWk <- t(X) %*% f
      tWk <- tWk / sqrt(sum(tWk^2))
      tTk <- X %*% tWk
      tTk <- tTk / sqrt(sum(tTk^2))
      Wk <- cbind(Wk, tWk)
      Tk <- cbind(Tk, tTk)
      Pt <- tTk %*% ginv(t(tTk)%*%tTk) %*% t(tTk)
      X <- X - Pt %*% X
   }
   Tk <- X0 %*% Wk # convert to the same def of above
   return(list(Wk=Wk, Tk=Tk, Xmean=Xmean)) 
}

## test code ##
#out1 <- pls(xyz[, inds$xyz], angs, k=10)
#out2 <- .pls(xyz[, inds$xyz], angs, k=10)
#m1 <- plsfit(out1, angs, xyz[, inds$xyz], k=2)
#m2 <- plsfit(out2, angs, xyz[, inds$xyz], k=2)
#testthat::expect_equal(m1$tmpz, m2$tmpz)
#testthat::expect_equal(m1$newv, m2$newv)
# All pass for k=1, 2, 10

## no centering?
## very different from pls()
.pls2 <- function(X, f, k=1) {
   require(MASS)
   Wk <- matrix( t(X) %*% f )
   Tk <- matrix( X %*% Wk )
   if(k > 1) {
      for (i in 2:k) {
#         Tk_1 <- Tk[, 1:(i-1), drop=FALSE]
#         tWk <- t(X) %*% (f - Tk_1 %*% ginv(t(Tk_1) %*% Tk_1) %*% t(Tk_1) %*% f)
         tWk <- t(X) %*% (f - Tk %*% ginv(t(Tk) %*% Tk) %*% t(Tk) %*% f)
         Wk <- cbind(Wk, tWk)
         Tk <- cbind(Tk, X %*% tWk)
      }
   }
   return(list(Wk=Wk, Tk=Tk)) 
}

# x is returned value from pls
plsfit <- function(x, f, xyz, k=1) {
   depx <- x$Tk[, 1:k, drop=FALSE]
   w <- x$Wk[, 1:k, drop=FALSE]
   colnames(depx) <- paste0("x", 1:ncol(depx))
   mydat <- as.data.frame(cbind(depx, omega=f))
   fit <- lm(omega ~., data=mydat)
   o <- summary(fit)
   beta <- o$coefficients[-1, 1]
   newv <- w %*% beta
   newv <- newv/sqrt(sum(newv^2))
#   tmpz <- sweep(xyz, 2, colMeans(xyz)) %*% newv
   tmpz <- sweep(xyz, 2, x$Xmean) %*% newv # sweeping not necessary; just different by a constant
   return(list(fit=fit, beta=beta, newv=newv, tmpz=tmpz, mean=x$Xmean))
}

## generate VMD arrows
gen_arrows <- function(x, pdb, refval=c(0.0, 0.18), mag=10, cut=1.0, file="mode_arrow.vmd") {
   source("vmd.modes.R")
   newv <- x$newv
   tmpz <- x$tmpz
   
   tpc <- list(U=matrix(newv, ncol=1), L=var(tmpz), z=matrix(tmpz, ncol=1), 
               mean = x$mean)
   class(tpc) <- "pca"

   ## take the maximal length within a residue
   au <- tapply(tpc$U[, 1], rep(pdb$atom$resno, each=3), function(x) {
#      sqrt(sum(x^2))
      max(sqrt(colSums(matrix(x, nrow=3)^2)))
   })
   #plot.bio3d(t(au), resno=pdb, typ="h")
   
#   mktrj(tpc, pc=1, pdb=pdb, mag=diff(range(tmpz))/2/sqrt(tpc$L[1]), b=vec2resno(au, pdb$atom$resno), file="pcr1.pdb")
   
   # set ref values to match a reference range
   # then set occup to trick VMD colormap
   # take the maximal length within a residue
   # should use occupancy to color backbone !!
   o <- tapply(tpc$U[, 1], rep(pdb$atom$resno, each=3), function(x) {
#      sqrt(mean(colSums(matrix(x, nrow=3)^2)))
       max(sqrt(colSums(matrix(x, nrow=3)^2)))
   })
   o <- vec2resno(o, pdb$atom$resno)
   o[c(1,3)] <- refval # do not change CA, which is #2
   vmd.modes(tpc, mode=1, pdb=pdb, mag=mag, cut=cut, color.map="linear", ref.0=refval[1], ref.1=refval[2], o=o, b=vec2resno(au, pdb$atom$resno), file=file)
}
 
## least square ###
lsfit <- function(mydat, pc, Ls, xyz, inds) {
   #fit <- lm(omega ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data=mydat)
   fit <- lm(omega ~., data=mydat)
   o <- summary(fit)
   #filter.i <- which(o$coefficients[, 4]<0.05)
   #filter.i <- which(o$coefficients[, 4]<1)
   #filter.i <- setdiff(filter.i, 1)

   #beta <- o$coefficients[filter.i, 1]
   beta <- o$coefficients[-1, 1]
   #newv <- pc$U[, filter.i-1] %*% beta
   newv <- pc$U[, Ls, drop=FALSE] %*% beta
#   newv <- pc$U[, 1:nL] %*% beta
   newv <- newv/sqrt(sum(newv^2))
   tmpz <- t(t(xyz[, inds$xyz]) - colMeans(xyz[, inds$xyz])) %*% newv
   return(list(fit=fit, beta=beta, newv=newv, tmpz=tmpz))
}

### LASSO ####
lasso <- function(mydat, pc, Ls, xyz, inds, cartesian=FALSE) {
   require(glmnet)
   grid <- 10^seq(10, -2, length=100)
   fit <- glmnet(as.matrix(mydat[, -ncol(mydat), drop=FALSE]), as.numeric(mydat[, ncol(mydat)]), alpha=1, lambda=grid)
   
   #set.seed(1)
   # use internal lambda
   #cv.out <- cv.glmnet(as.matrix(mydat[, 1:(ncol(mydat)-1)]), as.numeric(mydat[, ncol(mydat)]), alpha=1, lambda=grid)
   cv.out <- cv.glmnet(as.matrix(mydat[, -ncol(mydat), drop=FALSE]), as.numeric(mydat[, ncol(mydat)]), alpha=1)
#   plot(cv.out)
   #bestlam <- cv.out$lambda.1se
   bestlam <- cv.out$lambda.min
   #fit <- glmnet(as.matrix(mydat[, 1:(ncol(mydat)-1)]), as.numeric(mydat[, ncol(mydat)]), alpha=1, lambda=bestlam)
   coefs <- predict(fit, type="coefficients", s=bestlam)
#   coefs[coefs!=0]
   coefs <- coefs[-1]

   if(cartesian) {
      newv <- coefs/sqrt(sum(coefs^2))
   } else {
      newv <- pc$U[, Ls, drop=FALSE] %*% coefs
   }
#   newv <- pc$U[, 1:nL] %*% coefs
   newv <- newv/sqrt(sum(newv^2))
   tmpz <- t(t(xyz[, inds$xyz]) - colMeans(xyz[, inds$xyz])) %*% newv
   return(list(fit=fit, beta=coefs, cv.out=cv.out, bestlam=bestlam, newv=newv, tmpz=tmpz))
}

## Cross Validation
crossval <- function(mydat, nfold=10, method=c("lsq", "lasso"), bestlam=NULL, verbose=TRUE) {
   require(glmnet)
   method <- match.arg(method)
   mse <- NULL
   if(nfold > nrow(mydat)) {
      warning("nfold is larger than the number of data points")
      nfold <- nrow(mydat)
   }
   if(method == "lasso") {
      grid <- 10^seq(10, -2, length=100)
      if(is.null(bestlam)) {
         stop("Need a value for `bestlam`")
      }
   }
   foldid <- rep(1:nfold, length.out=nrow(mydat))
   for(j in 1:nfold) {
   #   print(which(foldid==10))
      if(verbose) cat(j, "...")
      i <- which(foldid==nfold)
      tmydat <- mydat[-i, ]
      mse <- c(mse, switch(method,
                lsq = {
                   fit <- lm(omega ~., data=tmydat)
                   pred <- predict(fit, newdata=mydat[i, ])
                   sum((pred - mydat$omega[i])^2)/length(i)
                },
                lasso = {
                   fit <- glmnet(as.matrix(tmydat[, -ncol(mydat), drop=FALSE]), as.numeric(tmydat[, ncol(mydat)]), alpha=1, lambda=grid)
                   coefs <- predict(fit, type="coefficients", s=bestlam)
                   pred <- as.matrix(cbind(1, mydat[i, -ncol(mydat), drop=FALSE])) %*% as.matrix(coefs)
                   sum((pred - mydat$omega[i])^2)/length(i)
                }
      ))
      foldid <- foldid + 1
      foldid[foldid>nfold] <- 1
   }
   if(verbose) cat("\n")
   return(list(rmse=sqrt(mean(mse)), mse=mse))
}

pca.array <- function(x) {
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
    mean <- colMeans(x)

    Q <- t(t(x) - mean)/sqrt(n - 1)
    prj <- svd(Q)
    L <- prj$d^2
    U <- prj$v
    
    L[L < 0] <- 0
    sdev <- sqrt(L)
    z <- sweep(x, 2, mean) %*% (U)
    class(U) = "pca.loadings"
    out <- list(L = L, U = U, z = z, sdev = sdev, mean = mean)
    class(out) = "pca"
    return(out)
}

## Cross Validation, version 2 (for PLS)
crossval2 <- function(xyz, f, k=1, nfold=10, verbose=TRUE) {
   if(nfold > nrow(xyz)) {
      warning("nfold is larger than the number of data points")
      nfold <- nrow(xyz)
   }
   foldid <- rep(1:nfold, length.out=nrow(xyz))
   preds <- matrix(NA, nrow=length(f), ncol=k)
   for(j in 1:nfold) {
   #   print(which(foldid==10))
      if(verbose) {
         cat(j, "...")
      }
      i <- which(foldid==nfold)
      plsdata <- pls(xyz[-i, ], f=f[-i], k=k)
      mypred <- sapply(1:k, function(kk) {
         out <- try( plsfit(plsdata, f[-i], xyz[i, ], k=kk), silent=TRUE)
         if(inherits(out, "try-error")) {
            return(rep(NA, length(i)))
         }
         newdata <- sweep(xyz[i, ], 2, plsdata$Xmean) %*% plsdata$Wk[, 1:kk, drop=FALSE]
         colnames(newdata) <- paste0("x", 1:ncol(newdata))
         newdata <- as.data.frame(newdata)
         pred <- predict(out$fit, newdata=newdata)
#         out$tmpz ## check if correlation bewteen predicted f and f on test set is 
                  ## equivalent to that between tmpz (mapped to PLS using training) and f on the test set. (YES) 
#         sum((pred - f[i])^2)
      })
#      preds <- rbind(preds, mypred)
      preds[i, ] <- mypred
#      sse <- cbind(sse, mysse)
      foldid <- foldid + 1
      foldid[foldid>nfold] <- 1
   }
   if(any(is.na(preds))) {
      ind <- min(which(is.na(preds), arr.ind=TRUE)[, 2]) - 1
      if(ind<=0) {
         stop("No fitting returned")
      }
      preds <- preds[, 1:ind, drop=FALSE]
   } 
   rmse <- sqrt(colMeans((preds - f)^2))
   corr <- apply(preds, 2, cor, f)
   if(verbose) cat("\n")
   return(list(rmse=rmse, corr=corr))
}

