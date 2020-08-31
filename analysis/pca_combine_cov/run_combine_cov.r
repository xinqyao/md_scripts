library(bio3d)

## Define number of frames for all simulations included
N <- c(10^5, 10^6) #...

## Define file names of average PDB
favg <- c("XXX/analysis/pca/avg.pdb", "YYY/analysis/pca/avg.pdb") #...

## Define file names of covariance matrix
fcov <- c("XXX/analysis/pca/covar.dat", "YYY/analysis/pca/covar.dat") #...


############ !! DO NOT CHANGE FOLLOWING LINES !! ####################

avgs <- lapply(favg, read.pdb)
covs <- lapply(fcov, read.table)

# new average structure
xyz <- do.call("rbind", lapply(avgs, "[[", "xyz"))
xyz <- as.xyz( colMeans(xyz) )

write.pdb(pdb=avgs[[1]], xyz=xyz, file="avg.pdb")

xys <- lapply(avgs, function(x) {
   t(x$xyz) %*% x$xyz
})

xy <- t(xyz) %*% xyz

out <- matrix(0, nrow = nrow(covs[[1]]), ncol = ncol(covs[[2]]))
for(i in 1:length(N)) {
   out <- out + (N[i]-1) / N[i] * covs[[i]] + xys[[i]] - xy
}
out <- as.matrix( out / (length(N) * (sum(N)-1) / sum(N)) )

write(format(round(out, 3), scientific=FALSE, justify="right", nsmall=3), file="covar.dat", ncolumns=ncol(out))

## Diagonalize matrix and write out results
avg <- xyz
prj <- eigen(out, symmetric=TRUE)
L <- prj$values
U <- prj$vectors

L[L<0] <- 0

f <- rep(1:ceiling(length(avg)/7), each=7)
f <- f[1:length(avg)]
cat(" Eigenvector file: COVAR nmodes 4 width 11\n", file="evecs.dat")
cat(" ", length(avg), " ", nrow(U), "\n",
   paste(sapply(split(sprintf("%11.5f", avg), f=f), paste, collapse=""), collapse="\n"), "\n",
   sep="", file="evecs.dat", append=TRUE)

for(i in 1:4) {
   cat(" ****\n",
      sprintf("%5d ", i), sprintf("%11.5f", L[i]), "\n",
      paste(sapply(split(sprintf("%11.5f", U[, i]), f=f), paste, collapse=""), collapse="\n"), "\n",
      sep="", file="evecs.dat", append=TRUE)
}
