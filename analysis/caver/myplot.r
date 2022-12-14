myplot <- 
function(x,
         key = TRUE,
         resnum.1 = c(1:ncol(x)),
         resnum.2 = c(1:nrow(x)),
         axis.tick.space = 20,
         zlim = range(x, finite = TRUE),
         nlevels = 20,
         levels = pretty(zlim, nlevels),
         color.palette = terrain.colors,
         col = rev(color.palette(length(levels) -1)),
         axes = TRUE,
         key.axes,
         xaxs = "i",
         yaxs = "i",
         las = 1,
         grid = FALSE,
         grid.col = "yellow",
         grid.nx = floor(ncol(x)/30),
         grid.ny = grid.nx,
         center.zero = FALSE,
         flip=TRUE,
         ...) {
  
  if (missing(x)) {
    stop("no 'x' distance matrix specified")
  }

  if(center.zero) {
    if(zlim[1]<0) {
      ## make levels equidistant around 0
      levels = pretty(c(-max(abs(zlim)),max(abs(zlim))), nlevels)
    }
  }
    
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))

  # Color key
  if(key) {
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
    
    par(las = las)
    mar <- mar.orig
    mar[4] <- mar[2]
    mar[2] <- 1
    par(mar = mar)

    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels),
                xaxs = "i", yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1], col = col)
    if (missing(key.axes)) {
      if (axes)
        axis(4)
    }
    else key.axes
    box()
  }    

  # Matrix plot
  mar <- mar.orig
  mar[4] <- 1
  par(mar = mar)
  class(x)=NULL
  
  z <- as.matrix(as.data.frame(t(x)))
  nums <- seq(1,ncol(x),by=axis.tick.space)
  nums2 <- seq(1,nrow(x),by=1)
  a2 <- resnum.2[nums2]
  
  if(flip) {
    z=as.matrix(rev(as.data.frame(t(x)))); a2 <- rev(resnum.2[nums2])
  }
  image(x=1:(ncol(x)),
        y=1:(nrow(x)),
        z=z,
        zlim=zlim,
        col=col, yaxt="n", xaxt="n", ...)
        #xlab="Residue Number", ylab="Residue Number")

#  axis(side=1, at=nums, labels=resnum.1[nums])
  axis(side=1, at=c(seq(1, 2000, 400), ncol(x)), labels=seq(0, 4, 0.8))
  axis(side=2, at=nums2, labels=a2)
  if(grid)
    grid(grid.nx ,grid.ny, col=grid.col)
  box()
}
