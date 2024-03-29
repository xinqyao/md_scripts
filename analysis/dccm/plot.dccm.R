plot.dccm <- function (x, resno = NULL, sse = NULL, colorkey = TRUE, at = c(-1, 
    -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1), main = "Residue Cross Correlation", 
    helix.col = "gray20", sheet.col = "gray80", inner.box = TRUE, 
    outer.box = FALSE, xlab = "Residue No.", ylab = "Residue No.", 
    margin.segments = NULL, segment.col = vmd_colors(), segment.min = 1, 
    ...) 
{
    requireNamespace("lattice", quietly = TRUE)
    colnames(x) = NULL
    rownames(x) = NULL
    if (!is.null(resno)) {
        if (is.pdb(resno)) {
            ca.inds <- atom.select(resno, "calpha", verbose = FALSE)
            resno <- resno$atom$resno[ca.inds$atom]
        }
        if (length(resno) != nrow(x)) {
            warning("Length of input 'resno' does not equal the length of input 'x'; Ignoring 'resno'")
            resno = NULL
        }
    }
    scales <- NULL
    dots <- list(...)
    if ("scales" %in% names(dots)) 
        scales <- dots$scales
    if(!"at" %in% names(scales)) {
       xy.at <- pretty(1:ncol(x))
       xy.at <- xy.at[xy.at <= ncol(x)]
       xy.at[1] <- 1
       if (is.null(resno)) {
           scales$at <- xy.at
           scales$labels <- xy.at
       }
       else {
           labs <- resno[xy.at]
           labs[is.na(labs)] <- ""
           scales$at <- xy.at
           scales$labels <- labs
       }
    }
    dots$scales <- scales
    draw.segment <- function(start, length, xymin, xymax, fill.col = "gray", 
        side = 1) {
        if (side == 1) {
            grid.rect(x = unit(start - 0.5, "native"), y = 0, 
                gp = gpar(fill = fill.col, col = NA), just = c("left", 
                  "bottom"), width = unit(length - 0.5, "native"), 
                height = xymin, vp = vpPath("plot_01.toplevel.vp", 
                  "plot_01.panel.1.1.vp"))
        }
        if (side == 2) {
            grid.rect(x = 0, y = unit(start - 0.5, "native"), 
                gp = gpar(fill = fill.col, col = NA), just = c("left", 
                  "bottom"), width = xymin, height = unit(length - 
                  0.5, "native"), vp = vpPath("plot_01.toplevel.vp", 
                  "plot_01.panel.1.1.vp"))
        }
        if (side == 3) {
            grid.rect(x = unit(start - 0.5, "native"), y = xymax, 
                gp = gpar(fill = fill.col, col = NA), just = c("left", 
                  "bottom"), width = unit(length - 0.5, "native"), 
                height = unit(1, "npc"), vp = vpPath("plot_01.toplevel.vp", 
                  "plot_01.panel.1.1.vp"))
        }
        if (side == 4) {
            grid.rect(x = xymax, y = unit(start - 0.5, "native"), 
                gp = gpar(fill = fill.col, col = NA), just = c("left", 
                  "bottom"), width = unit(1, "npc"), height = unit(length - 
                  0.5, "native"), vp = vpPath("plot_01.toplevel.vp", 
                  "plot_01.panel.1.1.vp"))
        }
    }
    p1 <- do.call(lattice::contourplot, c(list(x, region = TRUE, 
        labels = FALSE, col = "gray40", at = at, xlab = xlab, 
        ylab = ylab, colorkey = colorkey, main = main), dots))
    if (is.pdb(sse)) {
        sse <- pdb2sse(sse)
        sse <- bounds.sse(unname(sse))
    }
    if (length(sse$helix$start) == 0 && length(sse$sheet$start) == 
        0) 
        sse <- NULL
    xymin = 0
    xymax = 1
    if (is.null(sse) && is.null(margin.segments)) {
        print(p1)
    }
    else {
        xlim <- p1$x.limits
        ylim <- p1$y.limits
        uni <- 1/(max(xlim) - min(xlim))
        pad = 0.02
        padref <- pad/uni
        if (!is.null(sse)) {
            xymax <- 1 - (pad)
            p1$x.limits[2] = xlim[2] + padref
            p1$y.limits[2] = ylim[2] + padref
        }
        if (!is.null(margin.segments)) {
            xymin = pad
            p1$x.limits[1] = xlim[1] - padref
            p1$y.limits[1] = ylim[1] - padref
            grps <- table(margin.segments)
            grps = names(grps[grps > segment.min])
            store.grps <- NULL
            for (i in 1:length(grps)) {
                store.grps <- rbind(store.grps, cbind(bounds(which(margin.segments == 
                  grps[i])), grp = as.numeric(grps[i])))
            }
            if (is.null(segment.col)) {
                segment.col <- (store.grps[, "grp"])
            }
            else {
                segment.col <- segment.col[(store.grps[, "grp"])]
            }
        }
        print(p1)
        if (!is.null(sse)) {
            if (length(sse$helix$start) > 0) {
                if (is.null(sse$helix$length)) {
                  sse$helix$length <- (sse$helix$end + 1) - sse$helix$start
                }
                draw.segment(sse$helix$start, sse$helix$length, 
                  xymin = xymin, xymax = xymax, fill.col = helix.col, 
                  side = 3)
                draw.segment(sse$helix$start, sse$helix$length, 
                  xymin = xymin, xymax = xymax, fill.col = helix.col, 
                  side = 4)
            }
            if (length(sse$sheet$start) > 0) {
                if (is.null(sse$sheet$length)) {
                  sse$sheet$length <- (sse$sheet$end + 1) - sse$sheet$start
                }
                draw.segment(sse$sheet$start, sse$sheet$length, 
                  xymin = xymin, xymax = xymax, fill.col = sheet.col, 
                  side = 3)
                draw.segment(sse$sheet$start, sse$sheet$length, 
                  xymin = xymin, xymax = xymax, fill.col = sheet.col, 
                  side = 4)
            }
        }
        if (!is.null(margin.segments)) {
            draw.segment(store.grps[, "start"], store.grps[, 
                "length"], xymin = xymin, xymax = xymax, fill.col = segment.col, 
                side = 1)
            draw.segment(store.grps[, "start"], store.grps[, 
                "length"], xymin = xymin, xymax = xymax, fill.col = segment.col, 
                side = 2)
        }
        if (!outer.box) {
            grid.rect(x = 0, y = 0, gp = gpar(fill = NA, col = "white"), 
                just = c("left", "bottom"), width = 1, height = 1, 
                vp = vpPath("plot_01.toplevel.vp", "plot_01.panel.1.1.vp"))
        }
        if (inner.box) {
            grid.rect(x = xymin, y = xymin, gp = gpar(fill = NA, 
                col = "black"), just = c("left", "bottom"), width = xymax, 
                height = xymax, vp = vpPath("plot_01.toplevel.vp", 
                  "plot_01.panel.1.1.vp"))
        }
    }
}

