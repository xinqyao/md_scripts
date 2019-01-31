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
plot.fluct(rrs, col=c('black', 'red'), lab=c('WT', 'Mutant'), 
   sse=sse, typ='l', sse.min.length=3, ylab="RMSF (A)", axes=FALSE)
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
