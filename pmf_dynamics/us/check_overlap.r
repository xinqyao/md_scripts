library(ggplot2)
if(!dir.exists("overlap")) {
   dir.create("overlap")
}

skip=10

dirs <- dir(".", "^[0-9]+\\.[0-9]+$")
#dirs <- "5.6"
## MAY USE A LOT OF MEMORY IF USING MULTICORE !!
#for(j in c("3.5", "4.2", "4.9", "5.6")) {
#for(j in c("6.0", "6.3", "7.0", "14.0")) {
#for(j in c("6.0", "6.3", "7.0")) {
for(j in dirs) {
   cat("Processing ", j, "...", sep="")
   datfile <- list.files(j, pattern='^ang\\.[-0-9]*\\.out$', full.names=TRUE)
   i_prog <- 0
   dat <- parallel::mclapply(datfile, function(i) {
      gc()
      i_prog <<- i_prog + 1
      prog <- round(i_prog/length(datfile)*100, 0)
      if(prog %% 10 == 0) {
         cat(prog, "%...", sep="")
      }
      tdat <- read.table(i)
      ang <- sub('ang\\.(.*)\\.out', '\\1', basename(i))
      data.frame(y=tdat$V2[seq(1, nrow(tdat), skip)], ref=ang, stringsAsFactors=FALSE)
   }, mc.cores=1)
   dat <- do.call(rbind, dat)
   cat("Angles (", length(unique(dat$ref)), ")...", sep="")
   cat("Sample (", nrow(dat[dat$ref==dat$ref[1], ]), " rows)...", sep="")
   p <- ggplot(dat, aes(y, group=ref)) + geom_line(stat="density")
  
   pdf(file=paste0("overlap/overlap_", j, ".pdf"), width=6, height=6)
#   dev.copy2pdf(file="overlap_5.6.pdf")
   print(p)
   dev.off()

   cat("Done.\n")
}
