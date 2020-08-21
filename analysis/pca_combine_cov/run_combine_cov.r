library(bio3d)

## Define number of frames for all simulations included
N <- c(10^6, 2*10^6) #...

## data from the first simulation PCA
avg1 <- read.pdb("XXX/analysis/pca/avg.pdb")
cov1 <- read.table("XXX/analysis/pca/covar.dat")

## data from the second simulation PCA
avg2 <- read.pdb("YYY/analysis/pca/avg.pdb")
cov2 <- read.table("YYY/analysis/pca/covar.dat")

## data from more simulation PCA if included...
# avg3 <- read.pdb("ZZZ/analysis/pca/avg.pdb")
# cov3 <- read.table("ZZZ/analysis/pca/covar.dat")
# ...

xyz <- rbind(avg1$xyz, avg2$xyz) # ...
xyz <- as.xyz( colMeans(xyz) )

write.pdb(pdb=avg1, xyz=xyz, file="avg.pdb")


xy1 <- t(avg1$xyz) %*% avg1$xyz
xy2 <- t(avg2$xyz) %*% avg2$xyz
# ...

xy <- t(xyz) %*% xyz

## adapt following for more simulations
out <- (N[1]-1) / N[1] * cov1 +
       (N[2]-1) / N[2] * cov2 +
       xy1 + xy2 - 2*xy

out <- as.matrix( out / (2 * (sum(N)-1) / sum(N)) )

write(format(round(out, 3), scientific=FALSE, justify="right", nsmall=3), file="covar.dat", ncolumns=ncol(cov1))

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
