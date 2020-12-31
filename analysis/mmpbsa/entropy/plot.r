#library(ggplot2)
RL <- read.table('sample/sys1/PBSA_RL.dat')
R <- read.table('sample/sys1/PBSA_R.dat')
L <- read.table('sample/sys1/PBSA_L.dat')

ene <- RL - R - L

##
kb = 0.0019872041
temp = 300

## constant chsi
dgc <- -kb*temp * log((6.022*10^-4)/(8 * pi**2))

## entropy change up to 2nd order
s2 <- 1/(2*kb*temp)*var(ene[, 1])

## entropy change up to 3rd order
s3 <- 1/(3*2*(kb*temp)**2)*colMeans((ene - colMeans(ene))**3)

## Binding free energy up to 2nd order entropy
dg2 <- ene + s2 + dgc

## Binding free energy up to 3rd order entropy
dg3 <- dg2 + s3

cat("Mean binding free energy:\n")
colMeans(ene)
cat("\n")

cat("Mean binding free energy + confinement:\n")
colMeans(ene + dgc)
cat("\n")

cat("Mean binding free energy + 2nd entropy:\n")
colMeans(ene+s2)
cat("\n")

cat("Mean binding free energy + 2nd entropy + confinement:\n")
colMeans(ene+s2+dgc)
cat("\n")

cat("Mean binding free energy + 3rd entropy:\n")
colMeans(ene+s2+s3)
cat("\n")

cat("Mean binding free energy + 3rd entropy + confinement:\n")
colMeans(ene+s2+s3+dgc)
cat("\n")
