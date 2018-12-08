#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)

cmap <- c(obs="black", cur="blue", fut="red")


## Call as: Rscript plot.flux.R label obs cur fut png var txt

## provide input files for one var; others are assumed analogous
vars <- c("huss", "prec", "uas", "vas")

## png = name of output file for plotted figure
## txt = name of output file for plain-text table of metrics

#args <- c("flux rcp85 MPI-ESM-LR RegCM4 birmingham",
#          "data/obs/prec.METDATA.44i.birmingham.nc",
#          "data/gyro/prec.hist.MPI-ESM-LR.RegCM4.day.NAM-44i.gyro.birmingham.nc",
#          "data/gyro/prec.rcp85.MPI-ESM-LR.RegCM4.day.NAM-44i.gyro.birmingham.nc",
#          "test.flux.png",
#          "prec",
#          "test.flux.txt"
#          )

#args <- c("flux rcp85 MPI-ESM-LR RegCM4 birmingham",
#          "data/obs/prec.METDATA.44i.birmingham.nc",
#          "data/raw/prec.hist.MPI-ESM-LR.RegCM4.day.NAM-44i.raw.birmingham.nc",
#          "data/raw/prec.rcp85.MPI-ESM-LR.RegCM4.day.NAM-44i.raw.birmingham.nc",
#          "test.flux.png",
#          "prec",
#          "test.flux.txt"
#          )

args <- c("flux rcp85 MPI-ESM-LR RegCM4 birmingham",
          "data/obs/prec.METDATA.44i.birmingham.nc",
          "data/kddm/prec.hist.MPI-ESM-LR.RegCM4.day.NAM-44i.kddm.birmingham.nc",
          "data/kddm/prec.rcp85.MPI-ESM-LR.RegCM4.day.NAM-44i.kddm.birmingham.nc",
          "test.flux.png",
          "prec",
          "test.flux.txt"
          )



## comment out this line for testing
# args <- commandArgs(trailingOnly=TRUE)

label <- args[1]

basefiles <- c(obs=args[2], cur=args[3], fut=args[4])
infiles <- list()
for (v in vars){
    infiles[[v]] <- gsub(args[6], v, basefiles)
}

outfile <- args[5]
txtfile <- args[7]

nc <- lapply(infiles, function(x){lapply(x, nc_ingest)})



## extract variables of interest from the netcdf objects
data <- list()
for(v in names(nc)){
  data[[v]] <- lapply(nc[[v]], `[[`, v)
}



time <- list()
for(v in names(nc)){
  time[[v]] <- lapply(nc[[v]], `[[`, "time")
}

atime <- rapply(time, alignepochs, epoch="days since 1949-01-01", how="replace")

## Todo: code to check that times are identical.


## day of year
doy <- lapply(atime$prec, function(x){ceiling(x %% yearlength(x@calendar))})

## windows
Nwin <- 360
win <- lapply(doy, function(x){ceiling(Nwin * x / 366)})



## construct dataframe

df <- do.call(rbind, lapply(renest(data), as.data.frame))

df$window <- unlist(win)
df$run <- unlist(mapply(rep, names(win), sapply(win, length)))

## calculate moisture flux = huss * uas | vas

df$uflux <- df$huss * df$uas
df$vflux <- df$huss * df$vas


wetday <- function(x, threshold = 1){
    ifelse(x > threshold, "wet", "dry")
}

## There should be a way to do this with a builtin, but I haven't
## figured it out yet.

for(r in c("obs","cur","fut")){
    for(w in seq(Nwin)){
        mask <- df$run == r & df$window == w
        df$wetday[mask] <- wetday(df$prec[mask])
    }
}

## ordered factor to make sorting happen right
df$run <- factor(df$run, levels=c("obs", "cur", "fut"))

uflux <- aggregate(uflux ~ window * run * wetday, df, mean, na.rm=TRUE)
vflux <- aggregate(vflux ~ window * run * wetday, df, mean, na.rm=TRUE)

# uflux$dir <- "u"
# vflux$dir <- "v"
# 
# ## this is awkward enough that I'm thinking there must be a less
# ## roundabout way to approach this whole situation...
# 
# colnames(uflux) <- gsub("uflux","flux",colnames(uflux))
# colnames(vflux) <- gsub("vflux","flux",colnames(vflux))
# 
# flux <- rbind(uflux, vflux)
# 

## convert from long to wide for matplot
library(reshape)

library(zoo)
smooth <- function(x, k=31){
    N <- length(x)
    y <- rollapply(c(x,x), FUN=mean, width=k, fill=NA)
    y[c(N+seq(k),(k+1):N)]
}

u <- lapply(cast(uflux, window ~ wetday + run)[,-1], smooth)
v <- lapply(cast(vflux, window ~ wetday + run)[,-1], smooth)




## plotting

cmap <- c(obs="black", cur="blue", fut="red")
lmap <- rep(c(2,1), each=3)

rg <- c(-1,1)*max(abs(unlist(c(u,v))))

dev.new()
par(mfrow=c(2,2))
matplot(as.matrix(u[4:6]), type="l", col=cmap, lty=1, ylim=rg, lwd=2, main="u wet" )
matplot(as.matrix(v[4:6]), type="l", col=cmap, lty=1, ylim=rg, lwd=2, main="v wet" )
matplot(as.matrix(u[1:3]), type="l", col=cmap, lty=1, ylim=rg, lwd=2, main="u dry" )
matplot(as.matrix(v[1:3]), type="l", col=cmap, lty=1, ylim=rg, lwd=2, main="v dry" )



stop()


# positions for month labels
mons <- seq(12)*30-15
umon <- sapply(u, `[`, mons)
vmon <- sapply(v, `[`, mons)
moname <- c("J","F","M","A","M","J","J","A","S","O","N","D")
pos <- c(4, 4, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3)


matplot(u, v, type="l", col=cmap, lty=lmap, xlim=rg, ylim=rg, lwd=2 )

for(i in 4:6){
    points(x=umon[,i], y=vmon[,i], col=cmap[i-3], pch=16)
    text(x=umon[,i], y=vmon[,i], col=cmap[i-3], labels=moname, pos=pos)
}

abline(v=0, h=0, col="gray")

## metrics?

stop()

## Empty metrics dataframe

metrics <- data.frame(infile = "dummy", period="seas", analysis="joint",
                      ktau=0, q50over=0, modeP=0, modeT=0,
                      stringsAsFactors=FALSE)

## plotting

png(outfile, units="in", res=120, width=10, height=7)

#par(mfrow=c(2,2), oma=c(0,0,3,0), mgp=c(2,1,0), mar=c(3.5,3.5,3,1.5))
#
#for (s in seas){
#
#  sd <- jsdata[[s]]
#
#  mplot(sd, type="n", x="prec", y="temp", main=s, xlab=punits, ylab=tunits)
#
#  ## freezing line
#  abline(h=0, col="gray")
#
#  for(i in names(sd)){
#    d <- sd[[i]]
#  
#    ## gridded 2D density estimate
#    gf <- kde2d(d$prec, d$temp, lims=par("usr"))
#
#    
#    ## find contour levels containing 95% and 50% of total 2D-PDF AUC
#    ## (Ordering outer to inner makes plotting outliers easier.)  
#    plev <- c(0.95, 0.5)
#  
#    sz <- sort(gf$z)
#    auc <- sum(sz)
#    clev <- sz[findInterval((1-plev)*auc, cumsum(sz))]
#
#    
#    ## plot contours
#    contour(gf, levels=clev, drawlabels=FALSE,
#          lwd=seq(length(plev)), col=cmap[i], add=TRUE)
#
#    
#    ## add mode point
#    imode <- which(gf$z == max(sz), arr.ind=TRUE)
#    modept <- list(x=gf$x[imode[1]], y=gf$y[imode[2]])
#    points(modept$x, modept$y, pch=19, col=cmap[i], cex=1.5)
#    ## TODO: change mode point to something like 2D Tukey median
#    
#  
#    ## plot individual points outside 90% contour
#    ## sp::point.in.polygon has nicer syntax, but crashes R
#
#    a90 <- conpolyarr(gf, clev[1])
#      
#    outpts <- !in.out(a90, cbind(d$prec, d$temp))
#    points(d$prec[outpts], d$temp[outpts], col=cmap[i], pch=19, cex=0.2)
#
#    
#    
#    ## Calculate metrics
#
#    ## correlation between P & T  (Kendall's Tau)
#    ktau <- cor(d$prec, d$temp, method="kendall")
#
#    ## % overlap between obs & mod q50 contours
#    if (i == "obs"){
#      ## As long as obs is first, this will carry over to later loops      
#      grid <- as.matrix(expand.grid(gf$x, gf$y, KEEP.OUT.ATTRS=FALSE))
#      oqa50 <- conpolyarr(gf, clev[2])      
#      oq50io <- in.out(oqa50, grid)
#    }
#
#    qa50 <- conpolyarr(gf, clev[2])
#    q50io <- in.out(qa50, grid)
#    q50over <- sum(q50io & oq50io) / sum(q50io)
#
#    ## Add to metrics df
#    metrics <- rbind(metrics,
#                     list(infiles$prec[i], s, "joint",
#                          ktau, q50over,
#                          modept$x, modept$y))
#  }
#}
#
#
#mtext(paste(label, "joint PDFs"), line=1, outer=TRUE)
#
### centered legend outside plots
#
#par(fig=c(0,1,0,1), new=TRUE)
#plot(c(0,1), c(0,1), type="n", bty="n", xaxt="n", yaxt="n", ann=FALSE)
#
#nlev <- length(plev)
#legend(0.5, 0.5, ncol=2, bty="n", seg.len=1, x.intersp=0.5, xjust=0.5, yjust=0.65,
#       c(names(cmap), paste0(plev*100,"%"), "mode"),
#       fill=c(cmap, rep(NA, nlev+1)), border=NA,
#       lwd=c(NA,NA,NA,seq(nlev),NA),
#       pch=c(rep(NA,nlev+3),19)
#       )
# 
#dev.off()
#
#
### write out metrics
#
#metrics <- metrics[-1,]
#
#write.table(format(metrics, trim=TRUE, digits=3),
#            file=txtfile, quote=FALSE, sep="\t", row.names=FALSE)
#
