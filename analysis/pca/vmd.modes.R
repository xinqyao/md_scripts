vmd.nma <- function(...)
  vmd.modes(...)

vmd.pca <- function(...)
  vmd.modes(...)

vmd.modes <- function(modes, mode=NULL, pdb=NULL, file=NULL, scale=50, dual=FALSE,
                      color=c('red', 'green', 'blue'), ca.only=FALSE,
                      type="script", exefile = "vmd", ...) {

  if(!( (inherits(modes, "nma") || inherits(modes,"pca")) ))
    stop("must supply a 'nma' or 'pca' object, i.e. from 'nma()' or 'pca.xyz()'")
 
  color <- match.arg(color)
 
  allowed <- c("script", "launch")
  if(!type %in% allowed) {
    stop(paste("input argument 'type' must be either of:",
               paste(allowed, collapse=", ")))
  }
    
  ## Check if the program is executable
  if(type %in% c("launch")) {
      
      ## determine path to exefile
      exefile1 <- .get.exepath(exefile)
      
      ## Check if the program is executable
      success <- .test.exefile(exefile1)
      
      if(!success) {
          stop(paste("Launching external program failed\n",
                     "  make sure '", exefile, "' is in your search path", sep=""))
      }
      exefile <- exefile1
  }
    
  if(inherits(modes, "nma")) {
    if(is.null(mode))
      mode <- 7
    xyz <- modes$xyz
    mode.vecs <- matrix(modes$modes[,mode], ncol=3, byrow=T)
  }
  else {
    if(is.null(mode))
      mode <- 1
    xyz <- modes$mean
    mode.vecs <- matrix(modes$U[,mode], ncol=3, byrow=T)
  }

  xyz0 <- xyz 
  if(ca.only && !is.null(pdb)) {
     ## merge vectors for each residue
     mode.vecs <- tapply(1:nrow(mode.vecs), pdb$atom$resno, function(i) {
        colMeans(mode.vecs[i, ])
     })
     mode.vecs <- do.call(rbind, mode.vecs)
     xyz <- xyz[atom.select(pdb, "calpha")$xyz]
#     pdb <- atom.select(pdb, "calpha", value=TRUE)
  }
 
  ## calc all vec lengths (for coloring later)
  all.lens <- apply(mode.vecs, 1, function(x) sqrt(sum(x**2)))

  ## output file name
  if(is.null(file)) {
    if(type=="script")
      file <- "R.vmd"
  }
  
  ## use temp-dir unless we output a VMD script
  if(type %in% c("launch"))
    tdir <- tempdir()
  else
    tdir <- "."

  vmdfile <- tempfile(tmpdir=tdir, fileext=".vmd")
  pdbfile <- tempfile(tmpdir=tdir, fileext=".pdb")
  
  ## start building VMD script
  scr <- c("proc vmd_draw_arrow {mol start end} {
   # an arrow is made of a cylinder and a cone
   set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
   graphics $mol cylinder $start $middle radius 0.15
   graphics $mol cone $middle $end radius 0.25
}")
  scr <- c(scr, paste("mol new ", 
    normalizePath(pdbfile, winslash='/', mustWork=FALSE),
    " type pdb", sep=""))
  scr <- c(scr, "mol delrep 0 top
mol representation NewCartoon
mol color ColorID 8
mol addrep top

display update off")
      
  ## define color range 
#  blues <- colorRamp(c("white", color))

  ## change color scale in VMD
  scr <- c(scr, paste("
set color_start [colorinfo num]
for {set i 0} {$i < 1024} {incr i} {
  color change rgb [expr $i + $color_start]", 
  switch(color, 
    "red" = "1 [expr 1.0 - 1.0 / 1023 * $i] [expr 1.0 - 1.0 / 1023 * $i]",
    "green" = "[expr 1.0 - 1.0 / 1023 * $i] 1 [expr 1.0 - 1.0 / 1023 * $i]",
    "blue" = "[expr 1.0 - 1.0 / 1023 * $i] [expr 1.0 - 1.0 / 1023 * $i] 1"
  )), "}" ) 


  ## Arrow widths
#  w.body <- 0.15; w.head <- 0.2

  for ( i in 1:nrow(mode.vecs)) {
    inds    <- atom2xyz(i) 
    coords  <- xyz[inds]
    
    ## For coloring (longest vec has length=1)
    tmp.len <- (sqrt(sum(mode.vecs[i,]**2))-min(all.lens)) / max(all.lens)

    if(tmp.len>1) {
      tmp.len <- 1
    } 
    else if(tmp.len<0) {
      tmp.len <- 0
    }
 
    col <- round(tmp.len * 1023, 0)
 
    ## Main vector
    tmp.vec <- mode.vecs[i,] * scale
    
#    ## For arrow head
#    if(sqrt(sum(tmp.vec**2))<1)
#      norm.vec <-  tmp.vec
#    else
#      norm.vec <- normalize.vector(mode.vecs[i,])
    
    ## Set vectors
    arrow.vec.a <- (coords + tmp.vec)
#    head.vec.a  <- (arrow.vec.a + (norm.vec))
    
    arrow.vec.b <- (coords - tmp.vec)
#    head.vec.b <- (arrow.vec.b - (norm.vec))
    
    a  <- paste("{", paste(coords,      collapse=" "), "}", sep="")
    b1 <- paste("{", paste(arrow.vec.a, collapse=" "), "}", sep="")
#    c1 <- paste(head.vec.a,  collapse=",")
    b2 <- paste("{", paste(arrow.vec.b, collapse=" "), "}", sep="")
#    c2 <- paste(head.vec.b,  collapse=",")
    
    
    ## Arrows
    scr <- c(scr, paste("draw color [expr", col, "+ $color_start]"))
    scr <- c(scr, paste("draw arrow", a, b1))
    if(dual)
      scr <- c(scr, paste("draw arrow", a, b2))
  }
  scr <- c(scr, "display update on")
 
  ## Write PDB structure file
  write.pdb(pdb=pdb, xyz=xyz0, file=pdbfile)
  
  ## Write VMD script or PDB with conect records
  write(scr, file=vmdfile, sep="\n")

  if(type %in% c("launch")) {
    args <- ""
    
    ## Open VMD 
    cmd <- paste(exefile, args, vmdfile)
    
    os1 <- Sys.info()["sysname"]
    if (os1 == "Windows") {
        status <- shell(paste(shQuote(exefile), args, vmdfile))
    }
    else {
        status <- system(cmd)
    }
    
    if(!(status %in% c(0,1))) {
        stop(paste("An error occurred while running command\n '",
                   exefile, "'", sep=""))
    }

  }

  if(type == "script") {
    file.copy(vmdfile, file, overwrite=TRUE)
    unlink(vmdfile)
    message(paste("VMD script written to file", file))
    invisible(file)
  }
  
}