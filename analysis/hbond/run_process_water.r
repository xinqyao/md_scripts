files <- list.files("water", ".*acceptor.*bridge.*dat$", full.names=TRUE)
lines <- sapply(files, function(x) {y <- readLines(x); y[grep("CH6", y)]})
dat <- sapply(lines, function(x) { do.call(rbind, strsplit(x, split=",")) })
dat <- do.call(rbind, dat)
fr <- as.numeric(sub("^ ([0-9]*) frames.", "\\1", dat[, 2]))
dat <- data.frame(id=dat[, 1], fr=fr, stringsAsFactors=FALSE)
dat$id <- sub("Bridge Res (.*)\\s+$", "\\1", dat$id)
collapsed <- tapply(dat$fr, dat$id, sum)
collapsed <- sort(collapsed, decreasing=TRUE)
dat2 <- data.frame(id=names(collapsed), fr=unname(collapsed))
write.table(dat2, file="bridged.dat")


## manually process triple bridges
