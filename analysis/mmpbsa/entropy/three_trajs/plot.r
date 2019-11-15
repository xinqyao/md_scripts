library(ggplot2)
cis <- read.table('p151_cis_first.mmpbsa')
ts <- read.table('p151_ts_first.mmpbsa')
trans <- read.table('p151_trans_first.mmpbsa')
state <- rep(c("cis", "ts", "trans"), c(nrow(cis), nrow(ts), nrow(trans)))
dat <- cbind(deltaG=rbind(cis, ts, trans), state=state)
colnames(dat)[1] <- "energy"

pdf(onefile=TRUE, file="pbsa_first.pdf", width=5, height=4)
p <- ggplot(dat, aes(x=energy, col=state)) + geom_line(stat="density")
print(p)
dev.off()

cis <- read.table('p151_cis_second.mmpbsa')
ts <- read.table('p151_ts_second.mmpbsa')
trans <- read.table('p151_trans_second.mmpbsa')
state <- rep(c("cis", "ts", "trans"), c(nrow(cis), nrow(ts), nrow(trans)))
dat <- cbind(deltaG=rbind(cis, ts, trans), state=state)
colnames(dat)[1] <- "energy"

pdf(onefile=TRUE, file="pbsa_second.pdf", width=5, height=4)
p <- ggplot(dat, aes(x=energy, col=state)) + geom_line(stat="density")
print(p)
dev.off()

cis <- read.table('p151_cis_third.mmpbsa')
ts <- read.table('p151_ts_third.mmpbsa')
trans <- read.table('p151_trans_third.mmpbsa')
state <- rep(c("cis", "ts", "trans"), c(nrow(cis), nrow(ts), nrow(trans)))
dat <- cbind(deltaG=rbind(cis, ts, trans), state=state)
colnames(dat)[1] <- "energy"

pdf(onefile=TRUE, file="pbsa_third.pdf", width=5, height=4)
p <- ggplot(dat, aes(x=energy, col=state)) + geom_line(stat="density")
print(p)
dev.off()

