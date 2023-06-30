## maximal gap shifting
mgshift <- function(x, n=73, bw=250) {
  require(circular)
  den <- density.circular(circular(x/180*pi), bw=bw, from=circular(-pi), to=circular(pi), n=n, kernel="vonmises")
  ind <- which.min(den$y)
  angle <- den$x[ind]/pi*180
  if(angle > 0) {
     angle.range <- c(angle - 360, angle)
  } else {
     angle.range <- c(angle, angle + 360)
  }
  x[x>max(angle.range)] <- x[x>max(angle.range)] - 360
  x[x<min(angle.range)] <- x[x<min(angle.range)] + 360 
  list(x=x, angle.range=angle.range)
}

dpca <- function(x, ...) {
   dat <- do.call("rbind", list(x, ...))
   nangles <- ncol(dat) - 1

   dat <- parallel::mclapply(dat[, -1], function(x) mgshift(x)$x, mc.cores=10)
   dat <- do.call("cbind", dat)

   mean <- colMeans(dat)
   S <- var(dat)
   prj <- eigen(S, symmetric = TRUE)
   L <- prj$values
   U <- prj$vectors

   L[L < 0] <- 0
   sdev <- sqrt(L)
   z <- sweep(dat, 2, mean) %*% U
   au <- U
   class(U) = "pca.loadings"
   
   out <- list(L = L, U = U, z = z, au = au, sdev = sdev, mean = mean)
   class(out) = "pca"

   ## write out data
   write(S, file="covar.dat", ncolumns = ncol(S))
   write.table(cbind(1:nrow(z), z[, 1:3]), file="project.dat", row.names=FALSE, col.names=c("#Frame", "Mode1", "Mode2", "Mode3"), quote=FALSE)
   cat("", file="evecs.dat")
   for(i in 1:4) {
      cat(" ****\n", file="evecs.dat", append=TRUE)
      cat("    ", i, "  ", L[i], "\n", sep="", file="evecs.dat", append=TRUE)
      write(U[, i], file="evecs.dat", ncolumns = 7, append=TRUE)
   }

   out
}


