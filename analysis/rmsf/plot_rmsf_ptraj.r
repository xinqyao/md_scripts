library(bio3d)

## pdb file for plotting SSE
pdb <- read.pdb("../../../prep/XXX_prot.pdb")
pdb <- trim(pdb, "protein")
sse <- pdb
if(is.null(sse$helix) && is.null(sse$sheet)) {
  if(check.utility("dssp")) {
     sse <- dssp(pdb)
  }
}

rr <- read.table('rms.out')
rr2 <- read.table('rms2.out')
#...

rrs <- rbind(rr[, 2], rr2[, 2])
pdf(onefile=TRUE, file='rmsf.pdf', width=3.8, height=3.2)
plot.fluct(rrs, col=c('black', 'red'), lab=c('WT', 'Mutant'), ylim2zero=FALSE,
   sse=sse, typ='l', sse.min.length=3, xlab="Residue Number", ylab="RMSF (A)", axes=FALSE)
axis(2)
box()
res <- pdb$atom$resno[pdb$calpha]
myres <- pretty(res, n=10)
if(! res[1] %in% myres) myres[1] <- res[1]
if(! res[length(res)] %in% myres) myres[length(myres)] <- res[length(res)]
xticks <- match(myres, res)
axis(1, at=xticks, labels=myres)
dev.off()

## write b factors
write.pdb(pdb, b=vec2resno(rrs[1, ], pdb$atom$resno), file="rmsf1.pdb")
write.pdb(pdb, b=vec2resno(rrs[2, ], pdb$atom$resno), file="rmsf2.pdb")


## differential rmsf
#drrs <- rbind(rr2[, 2] - rr[, 2], rr3[, 2] - rr[, 2], rr4[, 2] - rr[, 2])
drrs <- rbind(rr2[, 2] - rr[, 2])

pdf(onefile=TRUE, file='drmsf.pdf', width=6, height=4.5)
plot.fluct(drrs, col=c('red'), lab=c('Mutant - WT'), ylim2zero=FALSE,
   sse=sse, typ='l', sse.min.length=3, xlab="Residue Number", ylab="DRMSF (A)", axes=FALSE)
axis(2)
box()
res <- pdb$atom$resno[pdb$calpha]
myres <- pretty(res, n=10)
if(! res[1] %in% myres) myres[1] <- res[1]
if(! res[length(res)] %in% myres) myres[length(myres)] <- res[length(res)]
xticks <- match(myres, res)
axis(1, at=xticks, labels=myres)
dev.off()

## write b factors
write.pdb(pdb, b=abs(vec2resno(drrs[1, ], pdb$atom$resno)), o=vec2resno(drrs[1, ], pdb$atom$resno), file="drmsf.pdb")



#####################################################################################
## Following is the actual scaling procedure in Pymol for cartoon_putty_transform 0
range <- 2.0
power <- 1.5
ndr <- (range + (abs(rrs[1, ]) - mean(abs(rrs[1, ]))) / sd(abs(rrs[1, ])))/range
ndr[ndr<0] <- 0
ndr <- ndr**(power)

## Now, we use some universal scaling factors
ndr <- (range + (abs(rrs) - 0.7) / 0.3)/range
ndr[ndr<0] <- 0
ndr <- ndr**(power)

write.pdb(pdb, o=vec2resno(rrs[1, ], pdb$atom$resno), b=vec2resno(ndr[1, ], pdb$atom$resno), file="rmsf_scaled1.pdb")
write.pdb(pdb, o=vec2resno(rrs[2, ], pdb$atom$resno), b=vec2resno(ndr[2, ], pdb$atom$resno), file="rmsf_scaled2.pdb")

ndr <- (range + (abs(drrs[1, ]) - 0.07) / 0.1)/range
ndr[ndr<0] <- 0
ndr <- ndr**(power)

write.pdb(pdb, b=abs(vec2resno(ndr, pdb$atom$resno)), o=vec2resno(drrs[1, ], pdb$atom$resno), file="drmsf_scaled.pdb")
