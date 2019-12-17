files <- list.files("../../restrt", "production.*.restrt_[0-9]+$", recursive=TRUE, full.name=TRUE)

chunks <- sub("production.([0-9]+).*", "\\1", basename(files))

lfiles <- tapply(files, chunks, function(x) {
   frame <- as.numeric(sub("production\\.[0-9]+\\.restrt_([0-9]+)", "\\1", basename(x)))
   x[order(frame)]
})

if(!dir.exists("water")) {
   dir.create("water")
}
tmp <- lapply(names(lfiles), function(x) {
  ofile <- paste("water/ptraj_", x, ".in", sep="")
  cat("parm ../../../prep/Rcatcher_lig_box.prmtop\n", file=ofile)
  for(i in lfiles[[x]]) {
     cat("trajin ", i, "\n", file=ofile, sep="", append=TRUE) 
  }
  cat("
center :1-217
image center familiar
hbond :1-61,63-217 donormask :62@N=,O= solventacceptor :CA|:WAT@O|@Na+,Cl- solventdonor :CA|:WAT|@Na+,Cl- bridgeout water/ch6_donor_bridge_", x, ".dat
hbond :1-61,63-217 acceptormask :62@N=,O= solventdonor :CA|:WAT|@Na+,Cl- solventacceptor :CA|:WAT@O|@Na+,Cl- bridgeout water/ch6_acceptor_bridge_", x, ".dat
go", file=ofile, sep="", append=TRUE)
})

## cpptraj -i water/ptraj_XXX.in
system('for i in water/ptraj_[0-9]*.in; do cpptraj -i $i &> water/err_$(basename $i .in).log; done')
