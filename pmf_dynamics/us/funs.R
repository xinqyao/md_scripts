## Identify wells and peaks (with some tolerance of noises)
get_wells_peaks <- function(x, tol=0) {
   require(bio3d)
   if(!is.numeric(x) | length(x) == 0) {
      stop("x is non-numeric or has zero length")
   }
   rl <- rle2(sign(diff(x)))
   while(min(rl$lengths) <= min(tol, sum(rl$lengths)-1)) {
      ind <- which.min(rl$lengths)
      rl$values[ind] <- -1 * rl$values[ind]
      rl <- rle2(inverse.rle(rl))
   }
   # remove the last chunk
   rl$lengths <- rl$lengths[-length(rl$lengths)]
   rl$values <- rl$values[-length(rl$values)]
   rl$inds <- rl$inds[-length(rl$inds)]
   
   if(length(rl$lengths)==0) {
      warning("No wells or peaks detected")
      return(NULL)
   }
   peaks <- rl$inds[rl$values==1] + 1
   wells <- rl$inds[rl$values==-1] + 1
   
   return(list(peaks=peaks, wells=wells))
}

# calculate survival probability using the K-M method
#  and 95% confidence interval using the exponential Greenwood method
get_survival_km <- function(x, unit=1, t=NULL) {
   x <- x[order(x[, 2]), ]
   ti <- c(0, x[x[, 3]==1, 2])
   if(is.null(t)) {
      t <- ti
   }
   di <- rep(1, length(ti)); di[1] <- 0
   ci <- rep(0, length(ti)) # censored 
   cti <- x[x[, 3]!=1, 2] 
   if(length(cti)>0) {
      ind <- table(sapply(cti, function(i) {j <- which(ti<=i); j[length(j)]}))
      ci[as.numeric(names(ind))] <- as.numeric(ind)
   }

   # risk set
   ni <- nrow(x)
   for(i in 2:length(ti)) {
      ni_1 <- ni[length(ni)]
      ni <- c(ni, ni_1 - di[i-1] - ci[i-1])
   }

   s <- sapply(t, function(t) prod((1-di/ni)[ti<=t]))

   v <- sapply(1:length(t), function(i) {
      rm.inds <- (ni==di | ni==0 | s[i]==1)
      di <- di[ti<=t[i] & !rm.inds]
      ni <- ni[ti<=t[i] & !rm.inds]
      if(length(di)==0) {
         NA
      } else {
         1/((log(s[i]))^2) * sum(di/(ni*(ni-di)))
      }
   })
   se_min <- exp(-exp(log(-log(s)) + 1.96 * sqrt(v))) 
   se_max <- exp(-exp(log(-log(s)) - 1.96 * sqrt(v)))
  
   lse_min <- -exp(log(-log(s)) + 1.96 * sqrt(v)) 
   lse_max <- -exp(log(-log(s)) - 1.96 * sqrt(v)) 
   
   data.frame(time=t*unit, prob=s, se_min=se_min, se_max=se_max, lse_min=lse_min, lse_max=lse_max)
}
