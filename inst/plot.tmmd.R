#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)

## Call as: Rscript --vanilla tmmd.R label tmax.obs tmax.cur tmax.fut out

## tmin files are assumed analogous to tmax files.  dtr is calculated
## as tmax-tmin

args <- commandArgs(trailingOnly=TRUE)

# # for testing
# args <- c("rcp85 HadGEM2-ES WRF ftlogan",
#           "obs/tmax.obs.livneh.ftlogan.nc",
#           "raw/tmax.hist.HadGEM2-ES.WRF.ftlogan.nc",
#           "raw/tmax.rcp85.HadGEM2-ES.WRF.ftlogan.nc",
#           "tmax.test.tmmd.png")
#
# args <- c("rcp85 WRF GFDL-ESM2M armsite",
#          "obs/tmax.obs.livneh.armsite.nc",
#          "save-test/tmax.hist.GFDL-ESM2M.WRF.armsite.nc",
#          "save-test/tmax.rcp85.GFDL-ESM2M.WRF.armsite.nc",
#          "fig-save/tmmd/tmmd.rcp85.WRF.GFDL-ESM2M.armsite.png",
#          "tmax")


label <- args[1]
        
maxfiles <- c()
maxfiles["obs"] <- args[2]
maxfiles["cur"] <- args[3]
maxfiles["fut"] <- args[4]

outfile <- args[5]



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

matplot(x, tdata$tmin, type="l", col=ocf, lty=1, ylim=yr,
        xlab="day of year", ylab=units, main="mean daily Tmin")
abline(v=which.max(tdata$tmin[,"obs"]), col="gray")
abline(v=which.min(tdata$tmin[,"obs"]), col="gray", lty=2)

legend("topright", names(ocf), col=ocf, lty=1, horiz=TRUE)
legend("topleft", c("max", "min"), col="gray", lty=c(1,2), horiz=TRUE)

matplot(x, tdata$dtr,  type="l", col=ocf, lty=1,
        xlab="day of year", ylab=units, main="mean daily DTR")
abline(v=which.max(tdata$dtr[,"obs"]), col="gray")
abline(v=which.min(tdata$dtr[,"obs"]), col="gray", lty=2)

mtext(label, line=1, outer=TRUE)
    
dev.off()


