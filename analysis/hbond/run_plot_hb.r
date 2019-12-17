library(ggplot2)

keys <- read.table("key_hbs.txt", header=TRUE, stringsAsFactors=FALSE)
dat <- read.table("ch6_acceptor_detail.out")
head <- readLines("ch6_acceptor_detail.out", n=1)
head <- strsplit(head, split="\\s+")[[1]]
colnames(dat) <- head
#dat <- dat[, keys$ID[!keys$Bridged]]

res <- sapply(strsplit(colnames(dat)[-1], split="-"), function(x) {
   sub("@.*", "", x[2])
})

dat <- dat[, c(FALSE, res %in% keys$ID)]
res <- res[res %in% keys$ID]

collapsed <- matrix(0, nrow=nrow(dat), ncol=length(unique(res)))
colnames(collapsed) <- unique(res)
for(i in colnames(collapsed)) {
   for(j in which(res %in% i)) {
      collapsed[, i] <- collapsed[, i] | dat[, j]
   }
}
for(i in keys$ID[!keys$ID %in% res]) {
   collapsed <- cbind(collapsed, 0)
}
colnames(collapsed) <- c(unique(res), keys$ID[!keys$ID %in% res])

chunks <- rep(1:20, each=100000)
hb <- tapply(1:nrow(collapsed), chunks, function(i) colMeans(collapsed[i, ]))
mydat <- as.data.frame(do.call(rbind, hb))
mydat <- cbind(Time=rep(seq(0.1, 2, 0.1), ncol(mydat)), reshape::melt(mydat))
mydat <- cbind(mydat, Bridged=rep(FALSE, nrow(mydat)), System=rep("Holo", nrow(mydat)))

files <- list.files("water", ".*acceptor.*bridge.*", full.names=TRUE)
strings <- sapply(strsplit(keys$ID, split="_"), function(x) paste(x[2], x[1], sep=":"))

bridges <- lapply(strings, function(ss) {
   sapply(files, function(x) {
      lines <- readLines(x)
      inds <- grep("CH6", lines)
      inds2 <- grep(ss, lines)
      inds <- intersect(inds, inds2)
      str <- strsplit(lines[inds], split=",")
      sapply(str, function(y) as.numeric(sub("frames.", "", y[2]))
   ) })
})
bridges2 <- lapply(bridges, function(x) sapply(x, function(y) if(length(y)==0) 0 else sum(y)))
bridges2 <- bridges2[match(colnames(collapsed), keys$ID)]
tmpdat <- data.frame(Time=rep(seq(0.1, 2, 0.1), length(bridges2)), variable=rep(colnames(collapsed), each=20), 
   value=unname(unlist(bridges2))/1000, Bridged=rep(TRUE, 20*length(bridges2)), System=rep("Holo", 20*length(bridges2)))
mydat <- rbind(mydat, tmpdat)

# read apo
dat <- read.table("../../../apo/production/analysis/hbond/ch6_acceptor_detail.out")
head <- readLines("../../../apo/production/analysis/hbond/ch6_acceptor_detail.out", n=1)
head <- strsplit(head, split="\\s+")[[1]]
colnames(dat) <- head
#dat <- dat[, keys$ID[!keys$Bridged]]

res <- sapply(strsplit(colnames(dat)[-1], split="-"), function(x) {
   sub("@.*", "", x[2])
})

dat <- dat[, c(FALSE, res %in% keys$ID)]
res <- res[res %in% keys$ID]

collapsed <- matrix(0, nrow=nrow(dat), ncol=length(unique(res)))
colnames(collapsed) <- unique(res)
for(i in colnames(collapsed)) {
   for(j in which(res %in% i)) {
      collapsed[, i] <- collapsed[, i] | dat[, j]
   }
}
for(i in keys$ID[!keys$ID %in% res]) {
   collapsed <- cbind(collapsed, 0)
}
colnames(collapsed) <- c(unique(res), keys$ID[!keys$ID %in% res])

chunks <- rep(1:20, each=100000)
hb <- tapply(1:nrow(collapsed), chunks, function(i) colMeans(collapsed[i, ]))
mydat2 <- as.data.frame(do.call(rbind, hb))
mydat2 <- cbind(Time=rep(seq(0.1, 2, 0.1), ncol(mydat2)), reshape::melt(mydat2))
mydat2 <- cbind(mydat2, Bridged=rep(FALSE, nrow(mydat2)), System=rep("Apo", nrow(mydat2)))

files <- list.files("../../../apo/production/analysis/hbond/water", ".*acceptor.*bridge.*", full.names=TRUE)
strings <- sapply(strsplit(keys$ID, split="_"), function(x) paste(x[2], x[1], sep=":"))

bridges <- lapply(strings, function(ss) {
   sapply(files, function(x) {
      lines <- readLines(x)
      inds <- grep("CH6", lines)
      inds2 <- grep(ss, lines)
      inds <- intersect(inds, inds2)
      str <- strsplit(lines[inds], split=",")
      sapply(str, function(y) as.numeric(sub("frames.", "", y[2]))
   ) })
})
bridges2 <- lapply(bridges, function(x) sapply(x, function(y) if(length(y)==0) 0 else sum(y)))
bridges2 <- bridges2[match(colnames(collapsed), keys$ID)]
tmpdat <- data.frame(Time=rep(seq(0.1, 2, 0.1), length(bridges2)), variable=rep(colnames(collapsed), each=20), 
   value=unname(unlist(bridges2))/1000, Bridged=rep(TRUE, 20*length(bridges2)), System=rep("Apo", 20*length(bridges2)))
mydat2 <- rbind(mydat2, tmpdat)
mydat <- rbind(mydat, mydat2)

p <- ggplot(mydat, aes(x=Time, y=value, color=System, linetype=Bridged)) + 
     geom_line() + 
     ylab("Probability") +
     theme(legend.position="bottom") +
     facet_wrap(~variable, ncol=3)
pdf(file="Rcatcher_hb.pdf", height=5, width=5)
print(p)
dev.off()

