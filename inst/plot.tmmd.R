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
args <- c("rcp85 HadGEM2-ES RegCM4 birmingham",
          "obs/tmax.obs.livneh.birmingham.nc",
          "save-test/tmax.hist.HadGEM2-ES.RegCM4.birmingham.nc",
          "save-test/tmax.rcp85.HadGEM2-ES.RegCM4.birmingham.nc",
          "test.tmmd.png",
          "tmax",
          "test.tmmd.txt"
          )

## comment out this line for testing
args <- commandArgs(trailingOnly=TRUE)


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

ncmax <- lapply(maxfiles, nc_ingest)
ncmin <- lapply(minfiles, nc_ingest)



## extract variables of interest from the netcdf objects
tmax <- lapply(ncmax,"[[","tmax")
tmin <- lapply(ncmin,"[[","tmin")
 dtr <- mapply(`-`, tmax, tmin)

units <- tmax$cur@units
    
time <- lapply(ncmax,"[[","time")
time <- lapply(time, alignepochs, "days since 1950-01-01")



## color palette for bias-correction
ocf <- c(obs="black", cur="blue", fut="red")


## width of moving window
mwinwidth = 30

### generate climatology moving window index arrays
### 360 ~1-day moving windows w/ 30-day outer pool
dwin <- lapply(time, cslice, outer=mwinwidth, num=360, split=FALSE)


data <- namelist(tmax, tmin, dtr)

cdata <- lapply(data, function(x){ mapply(slice, x, dwin,
                                          MoreArgs=list(outer=TRUE),
                                          SIMPLIFY=FALSE)})

capply <- function(x, f, ...){
  lapply(rapply(x, f, ..., how="replace"), unlist)
}


tdata <- lapply(lapply(cdata, capply, mean, na.rm=TRUE), as.matrix)

#sdata <- lapply(lapply(cdata, capply, sd), as.matrix)
#upper <- mapply(`+`, tdata, sdata, SIMPLIFY=FALSE)
#lower <- mapply(`-`, tdata, sdata, SIMPLIFY=FALSE)


## plotting

png(outfile, units="in", res=120, width=7, height=7)

par(mfrow=c(3,1), oma=c(0,0,3,0))

yr <- range(unlist(tdata[c("tmax","tmin")]))

x <- seq(1, 365, length=nrow(tdata$tmax))

matplot(x, tdata$tmax, type="l", col=ocf, lty=1, ylim=yr,
        xlab="day of year", ylab=units, main="mean daily Tmax")
abline(v=which.max(tdata$tmax[,"obs"]), col="gray")
abline(v=which.min(tdata$tmax[,"obs"]), col="gray", lty=2)
abline(h=0, col="gray")

matplot(x, tdata$tmin, type="l", col=ocf, lty=1, ylim=yr,
        xlab="day of year", ylab=units, main="mean daily Tmin")
abline(v=which.max(tdata$tmin[,"obs"]), col="gray")
abline(v=which.min(tdata$tmin[,"obs"]), col="gray", lty=2)
abline(h=0, col="gray")

legend("topright", names(ocf), col=ocf, lty=1, horiz=TRUE)
legend("topleft", c("max", "min"), col="gray", lty=c(1,2), horiz=TRUE)

matplot(x, tdata$dtr,  type="l", col=ocf, lty=1,
        xlab="day of year", ylab=units, main="mean daily DTR")
abline(v=which.max(tdata$dtr[,"obs"]), col="gray")
abline(v=which.min(tdata$dtr[,"obs"]), col="gray", lty=2)

mtext(label, line=1, outer=TRUE)
    
dev.off()



## metrics


mtemp <- lapply(tdata, function(x){as.list(as.data.frame(x))})


mcor <- lapply(mtemp, function(x){lapply(x, function(y){cor(y, x$obs)})})
mcor <- lapply(renest(mcor), unlist)
mmad <- lapply(mtemp, function(x){lapply(x, function(y){mad(y-x$obs, center=0)})})
mmad <- lapply(renest(mmad), unlist)


metrics <- data.frame(infile="dummy", period="ann", analysis="tmmd",
                      tmaxcor=0, tmincor=0, dtrcor=0,
                      tmaxmad=0, tminmad=0, dtrmad=0,
                      stringsAsFactors=FALSE)

for(p in names(maxfiles)){

  m0 <- list(infile=maxfiles[p], period="ann", analysis="tmmd")
 
  m1 <- as.list(mcor[[p]])
  names(m1) <- paste0(names(m1),"cor")

  m2 <- as.list(mmad[[p]])
  names(m2) <- paste0(names(m2),"mad")

  metrics <- rbind(metrics, c(m0, m1, m2))
}


## write out metrics

metrics <- metrics[-1,]

write.table(format(metrics, trim=TRUE, digits=3),
            file=txtfile, quote=FALSE, sep="\t", row.names=FALSE)

