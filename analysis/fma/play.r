library(bio3d)
library(ggplot2)
source("funs.R")
source("vmd.modes.R")

## read structural data (PDB) and functional variable data (here, we use the prolyl omega angle as an example)
files <- list.files("pdbs", full.names=TRUE)
angs <- as.numeric(sub("pdbs\\/production\\.([-0-9]*)\\.pdb", "\\1", files))
files <- files[order(angs)]
angs <- sort(angs)
pdbs <- lapply(files, read.pdb, multi=TRUE)
xyz <- do.call("rbind", lapply(pdbs, function(x) x$xyz))
xyz <- as.xyz(xyz)

## remove periodicity; define color code of angular regions for presentation
angs <- rep(angs, sapply(pdbs, function(x) nrow(x$xyz)))
angs[angs<=-150] <- angs[angs<=-150] + 360
cols <- rep("blue", length(angs))
cols[angs>50 & angs<=130] <- "red"
cols[angs>130] <- "darkgreen"

## Quick check of correlation between XYZ coordinates and the functional variable (omega)
inds <- atom.select(pdbs[[1]], resno=c(51:163), elety=c("N", "CA", "C", "O", "CB"))
range(apply(xyz[, inds$xyz], 2, cor, angs))

## remove the first three -- specific to the system
inds <- atom.select(pdbs[[1]], resno=c(54:163), elety=c("N", "CA", "C", "O", "CB"))
inds.fit <- atom.select(pdbs[[1]], resno=c(54:163), elety=c("N", "CA", "C", "O", "CB"))

## Do PCA - for comparison purpose only
pdb <- trim(pdbs[[1]], inds=inds)
xyz <- fit.xyz(xyz[1, ], xyz, inds.fit$xyz, inds.fit$xyz)
pc <- pca(xyz[, inds$xyz], use.svd=TRUE)
mktrj(pc, pc=1, pdb=pdb, file="pc1.pdb")
mktrj(pc, pc=2, pdb=pdb, file="pc2.pdb")
save(pc, file="pc.RData")

## Quick check of correlation between PCs and omega
signifL <- which(pc$L>.Machine$double.eps) # non-zero eigenvalues
rhos <- sapply(1:length(pc$L), function(i) {
   cor(pc$z[, i], angs)
})
range(rhos)
mydat <- data.frame(x=signifL, y=rhos[signifL])
p <- ggplot(mydat, aes(x, y)) + geom_col() +
       geom_hline(yintercept=c(2, -2)*sd(mydat$y), col="blue", linetype="dashed") +
       xlab("PC") +
       ylab("Correlation") +
       ylim(-0.5, 0.5) +
       theme_bw() +
       theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
pdf(file="pc_correlation.pdf", height=4, width=4)
print(p)
dev.off()

## Pick up the functional PC - the PC with maximal correlation - and display the mode
npc <- 2
refval <- c(0.0, 0.20)
o <- tapply(pc$U[, npc], rep(pdb$atom$resno, each=3), function(x) {
    max(sqrt(colSums(matrix(x, nrow=3)^2)))
})
o <- vec2resno(o, pdb$atom$resno)
range(o)
o[c(1,3)] <- refval
vmd.modes(pc, mode=npc, pdb, mag=25/sqrt(pc$L[npc]), cut=1.0, color.map="linear",
   ref.0=refval[1], ref.1=refval[2], o=o, file=paste0("pc", npc, "_arrow.vmd"))



## PARTIAL LEAST SQUQRE, PLS ###
ntrain <- 1:607 # 9/10 training
ntest <- 608:674 # 1/10 testing
plsdata <- pls(xyz[ntrain, inds$xyz], f=angs[ntrain], k=100)
corrs <- sapply(1:100, function(k) {
   cat("k=", k, "\n")
   out <- plsfit(plsdata, angs[ntrain], xyz[ntest, inds$xyz], k=k)
   cor(out$tmpz, angs[ntest])
}) 
#corrs.m <- sapply(1:100, function(k) {
#   cat("k=", k, "\n")
#   out <- plsfit(plsdata, angs[ntrain], xyz[ntrain, inds$xyz], k=k)
#   cor(out$tmpz, angs[ntrain])
#}) 
#plot(abs(corrs), typ="o", ylim=c(0, 1))
#lines(corrs.m, typ="o", col="red")

cvout <- crossval2(xyz[, inds$xyz], f=angs, k=100, nfold=10)
save(cvout, file="cvout.RData")
#load("cvout.RData")

#plot(cvout$rmse, typ="o", xlab="# of PLS components", ylab="RMSE")
#par(new=TRUE)
#plot(cvout$corr, typ="o", col="red", axes=FALSE, bty="n", xlab="", ylab="")
#axis(side=4)

## generate a figure
mydat <- data.frame(x=seq_along(cvout$rmse), rmse=cvout$rmse, corr=cvout$corr)
scale <- diff(range(mydat$corr)) / diff(range(mydat$rmse))
shift <- min(mydat$corr) - min(mydat$rmse) * scale
#col <- "#F8766D"
col <- "red"
col2 <- "gray60"
#col2 <- "black"
p <- ggplot(mydat, aes(x=x, y=corr)) +
        geom_line(color=col) +
        geom_point(shape=1, color=col) +
        geom_line(aes(x, rmse*scale+shift), col=col2) +
        geom_point(aes(x, rmse*scale+shift), shape=1, col=col2) +
        scale_y_continuous(
           limits = c(0.4, 0.6),
           name = "Correlation coefficient",
           breaks = seq(0.4, 0.6, 0.05),
           sec.axis = sec_axis(~(.-shift)/scale, name = "RMSE")
        ) +
#        xlim(1, 100) +
        xlab("# of PLS components") +
        annotate("text", x=which.max(mydat$corr)+3, y=max(mydat$corr)+.005, color=col,
           label=paste0("(", which.max(mydat$corr), ", ", format(round(max(mydat$corr), 3), nsmall=3), ")")) +
        theme_bw() +
        theme(
           axis.title.y.left = element_text(color = col),
           axis.ticks.y.left = element_line(color = col),
           axis.text.y.left = element_text(color = col)
        )

pdf(file="cvout.pdf", width=4, height=4)
   print(p)
dev.off()

k <- which.max(mydat$corr)
plsdata1 <- pls(xyz[, inds$xyz], f=angs, k=k)
out <- plsfit(plsdata1, angs, xyz[, inds$xyz], k=k)

range(sqrt(colSums(matrix(out$newv, nrow=3)^2)))
sd(out$tmpz)

# atomic loading cutoff = 0.04 (1/25)
gen_arrows(out, pdb, mag=25/sd(out$tmpz), refval=c(0.0, 0.20))

## use occupancy to color backbone!

save(out, file="fma_out.RData")
## END PARTIAL LEAST SQUQRE, PLS ###


# compare FM and PC
ov <- overlap(pc, out$newv)
mydat <- data.frame(x=seq_along(ov$overlap), Overlap=ov$overlap, Cummulative=ov$overlap.cum)
p <- ggplot(mydat, aes(x, y=Overlap)) + geom_col(width=0.1) +
     geom_line(aes(x, y=Cummulative), color="red") +
     theme_bw() +
     xlab("PC") +
     ylab("Overlap between FM and PC") +
#     ylim(0.0, 0.2) +
     annotate("text", x=5, y=0.02, label="Cummulative overlap", color="red")+
     annotate("text", x=5, y=0.025, label="Individual overlap", color="black")

pdf(file="ov_fm_pc.pdf", height=4, width=4)
print(p)
dev.off()

ov_self <- ov
save(ov_self, file="ov_self.RData")

