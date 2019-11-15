library(ggplot2)
cis <- read.table('../MMPBSA_p151_cis.mmpbsa')
ts <- read.table('../MMPBSA_p151_ts.mmpbsa')
trans <- read.table('../MMPBSA_p151_trans.mmpbsa')
state <- rep(c("cis", "ts", "trans"), c(nrow(cis), nrow(ts), nrow(trans)))

##
kb = 0.0019872041
temp = 300

## constant chsi
dgc <- -kb*temp * log((6.022*10^-4)/(8 * pi**2))

## entropy change up to 2nd order
s2.cis  <- 1/(2*kb*temp)*var(cis[, 1])
s2.ts <- 1/(2*kb*temp)*var(ts[, 1])
s2.trans  <- 1/(2*kb*temp)*var(trans[, 1])

## entropy change up to 3rd order
s3.cis  <- 1/(3*2*(kb*temp)**2)*colMeans((cis - colMeans(cis))**3)
s3.ts <- 1/(3*2*(kb*temp)**2)*colMeans((ts - colMeans(ts))**3)
s3.trans <- 1/(3*2*(kb*temp)**2)*colMeans((trans - colMeans(trans))**3)

## Binding free energy up to 2nd order entropy
dg2.cis <- cis + s2.cis + dgc
dg2.ts <- ts + s2.ts + dgc
dg2.trans <- trans + s2.trans + dgc

## Binding free energy up to 3rd order entropy
dg3.cis <- dg2.cis + s3.cis
dg3.ts <- dg2.ts + s3.ts
dg3.trans <- dg2.trans + s3.trans

dat2 <- cbind(deltaG=rbind(dg2.cis, dg2.ts, dg2.trans), state=state)
colnames(dat2)[1] <- "energy"

pdf(onefile=TRUE, file="pbsa_2nd_entropy.pdf", width=5, height=4)
p <- ggplot(dat2, aes(x=energy, col=state)) + geom_line(stat="density")
print(p)
dev.off()


dat3 <- cbind(deltaG=rbind(dg3.cis, dg3.ts, dg3.trans), state=state)
colnames(dat3)[1] <- "energy"

pdf(onefile=TRUE, file="pbsa_3rd_entropy.pdf", width=5, height=4)
p <- ggplot(dat3, aes(x=energy, col=state)) + geom_line(stat="density")
print(p)
dev.off()
