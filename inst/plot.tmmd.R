#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)

## Call as: Rscript plot.pfit.R label obs cur fut png var txt

## provide input files for tmax; tmin is assumed analogous
## png = name of output file for plotted figure
## txt = name of output file for plain-text table of metrics

# for testing
args <- c("tmmd rcp85 HadGEM2-ES RegCM4 boulder",
          "obs/tmax.obs.livneh.boulder.nc",
          "raw/tmax.hist.HadGEM2-ES.RegCM4.boulder.nc",
          "raw/tmax.rcp85.HadGEM2-ES.RegCM4.boulder.nc",
          "test.tmmd.png",
          "tmax",
          "test.tmmd.txt"
          )

## comment out this line for testing
#args <- commandArgs(trailingOnly=TRUE)




label <- args[1]
        
maxfiles <- c()
maxfiles["obs"] <- args[2]
maxfiles["cur"] <- args[3]
maxfiles["fut"] <- args[4]

outfile <- args[5]

varname <- args[6]
stopifnot(varname == "tmax")

txtfile <- args[7]

minfiles <- gsub("tmax", "tmin", maxfiles)

## pseudo-filenames used for metrics
dtrfiles <- gsub("tmax", "dtr", maxfiles)

ncmax <- lapply(maxfiles, nc_ingest)
ncmin <- lapply(minfiles, nc_ingest)



## extract variables of interest from the netcdf objects
tmax <- lapply(ncmax,"[[","tmax")
tmin <- lapply(ncmin,"[[","tmin")
 dtr <- mapply(`-`, tmax, tmin)

units <- tmax$cur@units
    
time <- lapply(ncmax,"[[","time")
time <- lapply(time, alignepochs, "days since 1950-01-01")



## width of moving window
mwinwidth = 30

## number of windows
N = 365

### generate climatology moving window index arrays
### N ~1-day moving windows w/ 30-day outer pool
dwin <- lapply(time, cslice, outer=mwinwidth, num=N, split=FALSE)

## day of year associated with each window (noleap calendar for simplicity)
day <- seq(1, 365, length=N)


vars <- c("Tmax", "Tmin", "DTR")

rawdata <- list(tmax, tmin, dtr)
names(rawdata) <- vars






## ought to merge this stuff into a windowed cslice apply...

cwindow <- function(x){
    mapply(slice, x, dwin, MoreArgs=list(outer=TRUE), SIMPLIFY=FALSE)
}

cdata <- lapply(rawdata, cwindow)

capply <- function(x, f, ...){
  lapply(rapply(x, f, ..., how="replace"), unlist)
}


traw <- lapply(cdata, capply, mean, na.rm=TRUE)
terr <- lapply(traw, function(x){
    list(cur=x$cur-x$obs, fut=x$fut-x$obs)})

tdata <- list(raw=traw, err=terr)

minval <- rapply(terr, min, how="replace")
maxval <- rapply(terr, max, how="replace")

minday <- rapply(terr, which.min, how="replace")
maxday <- rapply(terr, which.max, how="replace")



## plotting


## color palette for figures
ocf <- c(obs="black", cur="blue", fut="red")
cf  <- ocf[-1]


png(outfile, units="in", res=120, width=7, height=7)

par(mfrow=c(3,2), oma=c(0,0,3,0), mar=c(4,4,3,2), mgp=c(2.5,1,0))

for(v in vars){
    matplot(day, as.matrix(tdata$raw[[v]]), type="l", col=ocf, lty=1,
            xlab="day of year", ylab=units, main=paste("mean daily",v))
    abline(h=0, col="gray")

    
    matplot(day, as.matrix(tdata$err[[v]]), type="l", col=cf, lty=1,
            xlab="day of year", ylab=units, main=paste(v,"difference from obs"))
    abline(v=unlist(maxday[[v]]), col=cf, lty=2)
    abline(v=unlist(minday[[v]]), col=cf, lty=3)
    abline(h=0, col="gray")
}

mtext(label, line=1, outer=TRUE)

### legends in top corners outside plots
par(fig=c(0,1,0,1), mfrow=c(1,1), mar=c(0,0,0,0), oma=rep(0.5,4), new=TRUE)
plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n", ann=FALSE)    
legend("topright", names(ocf), col=ocf, lty=1,
       horiz=TRUE, cex=0.65, seg.len=1)
legend("topleft", c("max", "min"), col="gray", lty=c(2,3),
       horiz=TRUE, cex=0.65)

dev.off()


## metrics


mtemp <- lapply(terr, function(x){as.list(as.data.frame(x))})

## Average difference from obs.  Using hand-rolled mean absolute
## deviation / area under curve [auc()] instead of built-in (proper)
## median absolute deviation [mad()] because we also want to know what
## fraction of the area under the curve is positive.

auc <- function(x){mean(abs(x))}

ppos <- function(x){
    tot <- auc(x)
    x[x<0] <- 0
    auc(x) / tot * 100
}

mauc <- rapply(mtemp, auc, how="replace")
mpos <- rapply(mtemp, ppos, how="replace")


infiles <- list(Tmax=as.list(maxfiles[-1]),
                Tmin=as.list(minfiles[-1]),
                DTR =as.list(dtrfiles[-1]))

allmet <- list(dmad=mauc, dppos=mpos,
               dmaxday=maxday, dmaxval=maxval,
               dminday=minday, dminval=minval)

## rearrange to [var][per][metric]
allmet <- lapply(renest(allmet), renest)


metrics <- data.frame(infile="dummy", period="ann", analysis="tmmd",
                      dauc=0, dpctpos=0, dmaxday=0, dmaxval=0,
                      dminday=0, dminval=0, stringsAsFactors=FALSE)

for(v in vars){
    for(p in names(cf)){

        m0 <- list(infile=infiles[[v]][[p]], period="ann", analysis="tmmd")
        m1 <- allmet[[v]][[p]]

        metrics <- rbind(metrics, c(m0, m1))
    }
}



## write out metrics

metrics <- metrics[-1,]

write.table(format(metrics, trim=TRUE, digits=3),
            file=txtfile, quote=FALSE, sep="\t", row.names=FALSE)

