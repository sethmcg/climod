#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)

## Call as: Rscript plot.pfit.R label obs cur fut png var txt

## png = name of output file for plotted figure
## txt = name of output file for plain-text table of metrics

# for testing
args <- c("rcp85 HadGEM2-ES RegCM4 birmingham",
          "obs/prec.obs.livneh.birmingham.nc",
          "save-test/prec.hist.HadGEM2-ES.RegCM4.birmingham.nc",
          "save-test/prec.rcp85.HadGEM2-ES.RegCM4.birmingham.nc",
          "test.pfit.png",
          "prec",
          "test.pfit.txt"
          )

## comment out this line for testing
args <- commandArgs(trailingOnly=TRUE)


label <- args[1]
        
infiles <- c()
infiles["obs"] <- args[2]
infiles["cur"] <- args[3]
infiles["fut"] <- args[4]

outfile <- args[5]

v <- args[6]
stopifnot(v == "prec")

txtfile <- args[7]


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
    class(result) <- "pfit"
    return(result)
}


## pfit plotting
plot.pfit <- function(prec, fname, title=""){
    png(fname, units="in", res=120, width=7, height=7)
    par(mfrow=c(3,1), oma=c(0,0,3,0), mar=c(4,4,3,2), mgp=c(2.5,1,0))

    yint <- max(unlist(prec$int), na.rm=TRUE)
    ytot <- max(unlist(prec$tot))

    x <- seq(1, 365, length=length(prec$tot[,"obs"]))

    matplot(x, prec$freq, type="l", col=ocf, lty=1, ylim=c(0,100),
            xlab="day of year", ylab="% wet", main="mean precip frequency")
    abline(v=which.max(prec$freq[,"obs"]), col="gray")
    abline(v=which.min(prec$freq[,"obs"]), col="gray", lty=2)
    matplot(x, prec$int,  type="l", col=ocf, lty=1, ylim=c(0,yint),
            xlab="day of year", ylab="mm/day", main="mean precip intensity")
    abline(v=which.max(prec$int[,"obs"]), col="gray")
    abline(v=which.min(prec$int[,"obs"]), col="gray", lty=2)
    matplot(x, prec$tot,  type="l", col=ocf, lty=1, ylim=c(0,ytot),
            xlab="day of year", ylab="mm/day", main="mean total precip")
    abline(v=which.max(prec$tot[,"obs"]), col="gray")
    abline(v=which.min(prec$tot[,"obs"]), col="gray", lty=2)

    mtext(title, line=1, outer=TRUE)

    ## legends in top corners outside plots
    par(fig=c(0,1,0,1), mfrow=c(1,1), mar=c(0,0,0,0), oma=rep(0.5,4), new=TRUE)
    plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n", ann=FALSE)    
    legend("topright", names(ocf), col=ocf, lty=1,
           horiz=TRUE, cex=0.65, seg.len=1)
    legend("topleft", c("max", "min"), col="gray", lty=c(1,2),
           horiz=TRUE, cex=0.65)

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

prec <- pfit(pdata)

plot(prec, outfile, label)



## metrics

## omit RMSE
## rmse <- function(x){ sqrt(mean(x^2)) }

## Omit MASE
## mase <- function(x){ mean(abs(x)) / mean(abs(diff(x))) }


mprec <- lapply(prec, function(x){as.list(as.data.frame(x))})

mcor <- lapply(mprec, function(x){lapply(x, function(y){cor(y, x$obs)})})
mcor <- lapply(renest(mcor), unlist)
mmad <- lapply(mprec, function(x){lapply(x, function(y){mad(y-x$obs, center=0)})})
mmad <- lapply(renest(mmad), unlist)


metrics <- data.frame(infile="dummy", period="ann", analysis="pfit",
                      freqcor=0, intcor=0, totcor=0,
                      freqmad=0, intmad=0, totmad=0,
                      stringsAsFactors=FALSE)

for(p in names(infiles)){

  m0 <- list(infile=infiles[p], period="ann", analysis="pfit")
 
  m1 <- rev(as.list(mcor[[p]]))
  names(m1) <- paste0(names(m1),"cor")

  m2 <- rev(as.list(mmad[[p]]))
  names(m2) <- paste0(names(m2),"mad")

  metrics <- rbind(metrics, c(m0, m1, m2))
}


## write out metrics

metrics <- metrics[-1,]

write.table(format(metrics, trim=TRUE, digits=3),
            file=txtfile, quote=FALSE, sep="\t", row.names=FALSE)
