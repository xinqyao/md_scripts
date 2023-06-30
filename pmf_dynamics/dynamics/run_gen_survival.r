library(bio3d)

threslife <- 1
startover <- TRUE
paras <- c("4.9", "5.6", "6.3")
root <- "."
chkpoint <- c(90, 90, 90)
chkpoint2 <- c(270, 270, 270)
names(chkpoint) <- paras
names(chkpoint2) <- paras

#paras <- paras[1]

## if startover = FALSE
lchunk <- 4000000  # only check the last `lchunk` frames (10^6 = 20ns)



outdir <- paste0("ts_time_life", threslife)
if(!dir.exists(outdir)) {
   if(!startover) {
      cat('Output directory does not exist; Set startover = TRUE')
      startover <- TRUE
   }
   dir.create(outdir)
}

# init <- "trans"  # set manually
init <- NULL       # or detect automatically
if(is.null(init)) {
   dat <- read.table(file.path(file.path(root, paras[1]), "ang.1.out"), nrows=1)
   if(dat$V2[1] < 90) {
      init <- "cis"
   } 
   else {
      init <- "trans"
   } 
}

cat("\nThreshold: ", threslife, "\n")
cat("\nStartover: ", startover, "\n")
cat("Parameters: ", paste(paras, collapse=", "), "\n")
cat("Check point: ", paste(chkpoint[paras], collapse=", "), "\n")
cat("Check point2: ", paste(chkpoint2[paras], collapse=", "), "\n")
if(!startover) {
   cat("Chunk length: ", lchunk, "\n")
}
cat("Init: ", init, "\n\n")

cat("\n", date(), "\n\n")

for(para in paras) {
   cat("Processing ", para, " ...", sep="")
   files <- list.files(file.path(root, para), 'ang.*.out', full.names=TRUE)
   # sort files
   runs <- as.numeric( sub("ang.([0-9]+).out", "\\1", basename(files)) )
   files <- files[order(runs)]
   runs <- sort(runs)
   ofile <- paste(outdir, "/ts_time_", para, ".txt", sep="")
   if(!startover) {
      tbl <- suppressWarnings( tryCatch( read.table(ofile), error = function(e) {NULL} ) )
      if(!is.null(tbl)) {
         done_runs <- tbl$V1
         to_check <- tbl$V1[tbl$V3 == -1]
         if(length(to_check)>0) {
            done_runs <- setdiff(done_runs, to_check[ sapply(to_check, function(i) {
#                     tbl$V2[match(i, tbl$V1)] < length(readLines(file.path(para, paste0("ang.", i, ".out"))))
                     tbl$V2[match(i, tbl$V1)] < as.numeric(system(paste0("cat ", file.path(file.path(root, para), paste0("ang.", i, ".out")), " | wc -l"), intern=TRUE))
                  }) ])
         }
      } 
      else {
         done_runs <- NULL
      } 
      myruns <- setdiff(runs, done_runs)
      if(length(myruns) > 0) {
         files <- files[runs %in% myruns]
         runs <- myruns 
      } 
      else {
         cat("Done.\n")
         next
      }
   }
   if(startover || is.null(done_runs)) {
      cat("", file=ofile)
   }
   cat("Runs (", length(files), ")...", sep="")
   #tmpfile <- tempfile()
   for (i in 1:length(files)) {
     prog <- round(i/length(files)*100, 0)
     if(prog %% 10 == 0) {
        cat(prog, "...", sep="") 
     }

     dat <- read.table(files[i])$V2

     if(!startover) {
        tot <- length(dat)
        base <- tot - lchunk - 1 ## -1 to check pbc more properly
        if(base < 0) {
           base <- 0
        }
        dat <- dat[(base+1):tot]
     } else {
        base <- 0
     }

     ## minimal image
     df <- diff(dat)
     adf <- abs(df)
     if(any(adf > 180)) { # need to remove PBC
        pbc <- rep(0, length(dat))
        pbc[c(FALSE, adf>180)] <- sign(df[adf>180])
        pbc <- cumsum(pbc)
        dat <- dat - pbc*360
       
        # check if pbc is removed
        adf <- abs(diff(dat))
        if(any(adf > 180)) {
           stop(paste("PBC removal failed.", files[i], ", which:", which(adf>180)[1], ", base:",  base))
        }
     } # end remove PBC

     # start from cis (0)
     if(init == "cis") {
        if(abs(dat[1]) > 90) {
           stop(paste("Data processing failed.", files[i], ", base:",  base))
        }
        n <- which(dat > chkpoint[para])  #90)
        m <- which(dat < chkpoint2[para]-360) #-90)
     }

     # start from trans (180)
     if(init == "trans") {
        if(abs(dat[1]-180) > 90) {
           stop(paste("Data processing failed.", files[i], ", base:",  base))
        }
        n <- which(dat < chkpoint[para])  #90)
        m <- which(dat > chkpoint2[para]) #270)
     }

     if(length(n) > 0) {
        segn <- bounds(n)
        segn <- segn[segn[, "length"] >= threslife, ,drop=FALSE]
     }
     else {
        segn <- matrix(integer(0), ncol=3)
     }

     if(length(m) > 0) {
        segm <- bounds(m)
        segm <- segm[segm[, "length"] >= threslife, ,drop=FALSE]
     }
     else {
        segm <- matrix(integer(0), ncol=3)
     }
 
     if(nrow(segn)>0) {
        if(nrow(segm)==0 || (nrow(segm)>0 && segm[1, "start"]>segn[1, "start"])) {
           cat(runs[i], " ", ifelse(startover, 0, base) + segn[1, "start"], " 1\n", file=ofile, append=TRUE)
   #        x11()
   #        par(mfrow = c(2, 1))
   #        plot(dat, typ='p', pch=20, cex=0.3)
   #        abline(v=c(n[1], m[1]), lty=1, col='red')
   #        plot(dat, typ='p', pch=20, xlim=c(n[1]-100, n[1]+100), cex=0.3)
   ##        title(main=paste('life=',life, ", type=", type, sep=""))
   #        abline(v=c(n[1]), lty=1, col='red')
        }
        else {
           cat(runs[i], " ", ifelse(startover, 0, base) + segm[1, "start"], " 0\n", file=ofile, append=TRUE)
        }
     }
     else if(nrow(segm)>0) {
        cat(runs[i], " ", ifelse(startover, 0, base) + segm[1, "start"], " 0\n", file=ofile, append=TRUE)
     }
     else {
        cat(runs[i], " ", ifelse(startover, 0, base) + length(dat), " -1\n", file=ofile, append=TRUE)
     }
   }

   ## sort the table
   tbl <- read.table(ofile)
   tbl <- tbl[order(tbl$V1, tbl$V2, decreasing=TRUE), ]
   tbl <- tbl[!duplicated(tbl$V1), ]
   tbl <- tbl[order(tbl$V1), ]
   write.table(tbl, file=ofile, row.names=FALSE, col.names=FALSE) 

   cat("Done.\n")
#if(file.exists(ofile)) {
#   system(paste("paste", ofile, tmpfile, "> tmp"))
#   system(paste("mv tmp", ofile))
#}
#else {
#   file.copy(tmpfile, ofile)
#}
}
