dci <- function(x, y) {
   inds <- abs(x) >=0.1 | abs(y) >=0.1
   x <- x[inds]; y <- y[inds]
   inds <- x>=0.05 & y>=0.05 | x<=-0.05 & y<=-0.05
   sum(inds)/length(x)
}

df1 <- read.table("df1.txt")
df2 <- read.table("df2.txt")

dci(df1$V1, df2$V2)
