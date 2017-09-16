#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)

## Call as: Rscript --vanilla pfit.R label obs cur fut out

args <- commandArgs(trailingOnly=TRUE)

## for testing
#args <- c("rcp85 WRF MPI-ESM-LR elcentro",
#          "obs/prec.obs.livneh.elcentro.nc",
#          "save-test/prec.hist.MPI-ESM-LR.WRF.elcentro.nc",
#          "save-test/prec.rcp85.MPI-ESM-LR.WRF.elcentro.nc",
#          "fig-save/pfit/pfit.rcp85.WRF.MPI-ESM-LR.elcentro.png",
#          "prec"          
#          )


label <- args[1]
        
infiles <- c()
infiles["obs"] <- args[2]
infiles["cur"] <- args[3]
infiles["fut"] <- args[4]

outfile <- args[5]


v <- "prec"


## color palette for bias-correction
ocf <- c(obs="black", cur="blue", fut="red")


unzintensity <- function(y){
  if(all(is.na(y))){ 0 } else { mean(y, na.rm=TRUE) }
}

unzfrequency <- function(y){
  sum(!is.na(y))/length(y)*100
}

## Precip frequency, intensity, total
pfit <- function(ocfdata){

    wet <- lapply(ocfdata, function(x){lapply(x, unlist)})         
    zwet <- rapply(wet, rezero, how="replace")
    
    result <- list()
    
    result$tot  <- sapply(zwet, function(x){sapply(x, mean)} )
    result$int  <- sapply(wet, function(x){sapply(x, unzintensity)} )
    result$freq <- sapply(wet, function(x){sapply(x, unzfrequency)})
    return(result)
}


## pfit plotting
pfitplot <- function(ocfdata, fname, title=""){
    png(fname, units="in", res=120, width=7, height=7)
    par(mfrow=c(3,1), oma=c(0,0,3,0))

    prec <- pfit(ocfdata)
    yint <- max(unlist(prec$int), na.rm=TRUE)
    ytot <- max(unlist(prec$tot))

    x <- seq(1, 365, length=length(ocfdata$obs))

    matplot(x, prec$freq, type="l", col=ocf, lty=1, ylim=c(0,100),
            xlab="day of year", ylab="% wet", main="mean precip frequency")
    abline(v=which.max(prec$freq[,"obs"]), col="gray")
    abline(v=which.min(prec$freq[,"obs"]), col="gray", lty=2)
    legend("bottomright", names(ocf), col=ocf, lty=1, horiz=TRUE)
    legend("bottomleft", c("max", "min"), col="gray", lty=c(1,2), horiz=TRUE)
    matplot(x, prec$int,  type="l", col=ocf, lty=1, ylim=c(0,yint),
            xlab="day of year", ylab="mm/day", main="mean precip intensity")
    abline(v=which.max(prec$int[,"obs"]), col="gray")
    abline(v=which.min(prec$int[,"obs"]), col="gray", lty=2)
    matplot(x, prec$tot,  type="l", col=ocf, lty=1, ylim=c(0,ytot),
            xlab="day of year", ylab="mm/day", main="mean total precip")
    abline(v=which.max(prec$tot[,"obs"]), col="gray")
    abline(v=which.min(prec$tot[,"obs"]), col="gray", lty=2)

    mtext(title, line=1, outer=TRUE)
    
    dev.off()
}


        
nc <- lapply(infiles, nc_ingest)

## extract variables of interest from the netcdf objects
data <- lapply(nc,"[[",v)
units <- data$cur@units
    
time <- lapply(nc,"[[","time")

time <- lapply(time, alignepochs, "days since 1950-01-01")


## width of moving window
mwinwidth = 30

### generate climatology moving window index arrays
### 360 ~1-day moving windows w/ 30-day outer pool
dwin <- lapply(time, cslice, outer=mwinwidth, num=360, split=FALSE)

zdata <- lapply(data, unzero)
                
pdata <- mapply(slice, zdata, dwin, MoreArgs=list(outer=TRUE),
                SIMPLIFY=FALSE)

pfitplot(pdata, outfile, label)

