library(bio3d)
library(ggplot2)
source("myplot.r")

n <- 3   ## The number of tunnels to check for each system (simulation)
hfile <- c("wt/results", "n370s/results", "l444p/results")  ## Header of paths to results
#############

files <- file.path(hfile, "analysis/bottleneck_heat_maps/bottleneck_heat_map.csv")

data <- NULL
for(i in 1:length(files)) {
   dat <- as.matrix(read.csv(files[i], header=FALSE, sep=","))
   dat <- dat[, -ncol(dat)]
   dat <- dat[c(1:n), ]
   data <- rbind(data, dat)
}


## !! IMPORTANT: In the following plots, you might need to change some parameters !!

#myplot(data, zlim=range(data[data>0], finite=TRUE), color.palette=terrain.colors,
pdf("bottleneck_dynamics.pdf", width=3.23, height=3, pointsize=7)
myplot(data, zlim=c(1.0, 3.0), color.palette=terrain.colors, nlevels=20,
   xlab="Time (us)", ylab="Tunnel No.", #resnum.1=c(0:(ncol(data)))*0.2,
   axis.tick.space=20 )
dev.off()



pdf(onefile=TRUE, file='radius_probablity_density.pdf', width=5, height=5)
dat <- data.frame(x=c(apply(data[1:3, 151:1050], 2, max), 
                      apply(data[4:6, 151:1050], 2, max), 
                      apply(data[7:9, 501:2000], 2, max)), 
                      grp=factor(rep(c('wt', 'n370s', 'l444p'), c(900, 900, 1500)), levels=c('wt', 'n370s', 'l444p')))
dat$x[dat$x<0] <- 0
ggplot(dat, aes(x)) + geom_line(aes(col=grp), stat='density') +
#   scale_colour_discrete(guide=guide_legend(title=NULL)) +
   scale_colour_discrete(guide=FALSE) +
   xlab("Bottleneck Radius (A)") +
   ylab('Probability Density')
dev.off()



pdf(onefile=TRUE, file='histogram.pdf', width=5, height=5)
dat <- data.frame(x=c(as.numeric(apply(data[1:3, 151:1050], 2, max)>=1.5), 
                      as.numeric(apply(data[4:6, 151:1050], 2, max)>=1.5), 
                      as.numeric(apply(data[7:9, 501:2000], 2, max)>=1.5)), 
                      grp=factor(rep(c('wt', 'n370s', 'l444p'), c(900, 900, 1500)), levels=c('wt', 'n370s', 'l444p')))
dat0 <- tapply(dat$x, dat$grp, function(x) as.numeric(table(x))/length(x))
dat0 <- data.frame(x=unlist(dat0), grp=rep(c('wt', 'n370s', 'l444p'), each=2), state=rep(c('Closed', 'Open'), 3))

#dat$x[dat$x<0] <- 0
ggplot(dat0, aes(x=state, y=x)) + geom_bar(aes(fill=grp), stat='identity', position="dodge") +
#   scale_colour_discrete(guide=guide_legend(title=NULL)) +
#   scale_colour_discrete(guide=FALSE) +
   xlab("Tunnel") +
   ylab('Probability')
dev.off()



## survey of multiple threholds
pdf(onefile=TRUE, file='lines.pdf', width=5, height=5)
cutoffs <- seq(1, 2.5, 0.1)
dat1 <- matrix(0, length(cutoffs), 4)
for(i in 1:length(cutoffs)) { 
  dat <- data.frame(x=c(as.numeric(apply(data[1:3, 151:1050], 2, max)>=cutoffs[i]), 
                        as.numeric(apply(data[4:6, 151:1050], 2, max)>=cutoffs[i]), 
                        as.numeric(apply(data[7:9, 501:2000], 2, max)>=cutoffs[i])),
                        grp=factor(rep(c('wt', 'n370s', 'l444p'), c(900, 900, 1500)), levels=c('wt', 'n370s', 'l444p')))
  dat0 <- tapply(dat$x, dat$grp, function(x) as.numeric(table(x))/length(x))
  dat0 <- data.frame(x=unlist(dat0), grp=rep(c('wt', 'n370s', 'l444p'), each=2), state=rep(c('Closed', 'Open'), 3))
  dat1[i, ] <- dat0$x[dat0$state=="Open"] 
}
dat1 <- data.frame(x=rep(cutoffs, 3), y=as.vector(dat1), grp=factor(rep(c('wt', 'n370s', 'l444p'), 
 rep(length(cutoffs), 3)), levels=c('wt', 'n370s', 'l444p')))

#dat$x[dat$x<0] <- 0
ggplot(dat1, aes(x=x, y=y)) + geom_line(aes(color=grp)) +
#   scale_colour_discrete(guide=guide_legend(title=NULL)) +
#   scale_colour_discrete(guide=FALSE) +
   scale_colour_manual(values=c("blue", "green", "orange")) +
   xlab("Minimal bottleneck radius of open state") +
   ylab('Probability of open state')
dev.off()



profile <- read.csv(file.path(hfile[1], "analysis/tunnel_profiles_last_frame_cluster1.csv"), header=FALSE)
pdf(onefile=TRUE, file='profile.pdf', width=5, height=5)
dat <- data.frame(x=1:26, radius=unlist(profile[14:39]))
ggplot(dat, aes(x, radius)) + geom_line() +
#   scale_colour_discrete(guide=guide_legend(title=NULL)) +
#   scale_colour_discrete(guide=FALSE) +
   xlab("Tunnel") +
   ylab('Radius')
dev.off()

