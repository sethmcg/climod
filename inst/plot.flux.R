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

## ordered factor to make sorting happen correctly
df$run <- factor(df$run, levels=c("obs", "cur", "fut"))

uflux <- aggregate(uflux ~ window * run * wetday, df, mean, na.rm=TRUE)
vflux <- aggregate(vflux ~ window * run * wetday, df, mean, na.rm=TRUE)


## smooth and convert from long to wide for matplot
library(reshape)

library(zoo)
smooth <- function(x, k=31){
    N <- length(x)
    y <- rollapply(c(x,x), FUN=mean, width=k, fill=NA)
    y[c(N+seq(k),(k+1):N)]
}

u <- lapply(cast(uflux, window ~ wetday + run)[,-1], smooth)
v <- lapply(cast(vflux, window ~ wetday + run)[,-1], smooth)


## Ugh, gross.  There must be a better way to do this

wd <- list(wet=list(), dry=list())

flux <- list(u=wd,v=wd)

flux$u$wet <- u[4:6]
flux$v$wet <- v[4:6]
flux$u$dry <- u[1:3]
flux$v$dry <- v[1:3]

flux <- lapply(lapply(lapply(flux, renest), `names<-`, c("obs","cur","fut")), renest)



## plotting

cmap <- c(obs="black", cur="blue", fut="red")
lmap <- rep(c(2,1), each=3)

rg <- c(-1,1)*max(abs(unlist(c(u,v))))

png(outfile, units="in", res=120, width=7, height=7)

par(mfrow=c(2,2), oma=c(0,0,3,0), mgp=c(2,1,0), mar=c(3.5,3.5,3,1.5))

for(dir in c("u","v")){
    for(damp in c("wet","dry"){
    }
}

matplot(as.matrix(u[4:6]), type="l", col=cmap, lty=1, ylim=rg, lwd=2, main="u wet" )
matplot(as.matrix(v[4:6]), type="l", col=cmap, lty=1, ylim=rg, lwd=2, main="v wet" )
matplot(as.matrix(u[1:3]), type="l", col=cmap, lty=1, ylim=rg, lwd=2, main="u dry" )
matplot(as.matrix(v[1:3]), type="l", col=cmap, lty=1, ylim=rg, lwd=2, main="v dry" )



stop()

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
