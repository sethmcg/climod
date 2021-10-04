#library(climod)
library(devtools)
load_all("~/climod")
library(ncdf4)
library(colorspace)


## Call as: Rscript plot.trend.R label obs cur fut png var txt

## png = name of output file for plotted figure
## txt = name of output file for plain-text table of metrics

# for testing
args <- c("tmin rcp85 MPI-ESM-MR CRCM5-UQAM NAM-44i mbcn-gridMET",
          "test/tmin.obs.gridMET.day.NAM-44i.nc",
          "test/tmin.hist.MPI-ESM-MR.CRCM5-UQAM.day.NAM-44i.raw.nc",
          "test/tmin.rcp85.MPI-ESM-MR.CRCM5-UQAM.day.NAM-44i.raw.nc",
          "figs/test.trend.png",
          "tmin",
          "figs/test.trend.txt"
          )

## comment out this line for testing
#args <- commandArgs(trailingOnly=TRUE)

label <- args[1]
        
infiles <- c()
infiles["obs"] <- args[2]
infiles["cur"] <- args[3]
infiles["fut"] <- args[4]

outfile <- args[5]

v <- args[6]

txtfile <- args[7]


nc <- lapply(infiles, nc_ingest)

## extract variables of interest from the netcdf objects
Data <- lapply(nc,"[[", v)

units <- Data$cur@units


## 4-th root precip
if(v == "prec"){
  suppressWarnings(Data <- lapply(Data, '^', 1/4))
  units <- expression(sqrt(bquote(units), 4))
}


time <- lapply(nc,"[[","time")
epoch <- "days since 1900-01-01"
time <- lapply(time, alignepochs, epoch)


## Not proper months, just 1/12 of the year length

## decade = 1980 == [1976,1985]
## that way we don't get single-year decades at ends

## We'd want to alignepochs() in this function if we hadn't *just*
## done it.  Note that it may be a smidge off; need to add 0.001 to
## time to get fyear values that match NCL date function results.

ydm <- function(T){

    ylen   <- yearlength(T)
    fyear  <- T / ylen + 1900
    year   <- c(floor(fyear))
    decade <- signif(fyear, 3)
    month  <- ceiling(fyear %% 1 * 12 + 1e-10) ## 1e-10 to avoid month=0
    data.frame(year, decade, month)
}



## concatenate model current & future

tcat <- function(x){
    x$mod <- copyatts(x$fut, c(x$cur, x$fut))
    return(x)}

Data <- tcat(Data)
time <- tcat(time)

df <- mapply(cbind, lapply(time, ydm), data=Data, SIMPLIFY=FALSE)


agg <- function(x){aggregate(data ~ decade + month, x, mean)}
xt  <- function(x){xtabs(data ~ decade + month, x)}
trend <- lapply(lapply(df, agg), xt)

dx <- lapply(df, function(x){sort(unique(x$decade))})


## trend in units per century
slope <- mapply(function(x,y){sprintf("%3.2e",lm(y~x)$coeff[2,]*100)},
                dx, trend, SIMPLIFY=FALSE)


## make room for legend & slope annotations
extend <- function(x, low=0, high=0){
    w = diff(x)
    x + c(-w*low, w*high)
}
xr <- extend(range(unlist(dx)), lo=1/8, hi=1/8)
yr <- extend(range(unlist(trend)), hi=1/12)


cmap <- c(qualitative_hcl(12, "Dark3")[(12:1-6)%%12+1], "black")



png(outfile, units="in", res=120, height=6, width=9)
#dev.new(height=6, width=9)

par(mar=c(3,4.2,3,2), oma=c(0,0,3,0))



matplot(dx$mod, trend$mod, type='l', lty=1, lwd=2, col=cmap,
        xlim=xr, ylim=yr, axes=FALSE, main=v, xlab='', ylab=units)

matpoints(dx$obs, trend$obs, pch=1, col=cmap)

abline(v=min(df$fut$year), lty=3)

axis(1)
axis(2)



llwd <- c(rep(3,12),NA)
lpch <- c(rep(NA,12),1)

legend("top", horiz=TRUE, cex=0.7, seg.len=4/3, x.i=2/3,
       c(month.abb,"obs"), lwd=llwd, col=cmap, pch=lpch,
       text.width=diff(xr)/13*0.5)


text(max(unlist(dx)), tail(trend$fut,1), slope$fut,
     pos=4, col=cmap, cex=0.6)

text(min(unlist(dx)), tail(trend$cur,1), slope$cur,
     pos=2, col=cmap, cex=0.6)


mtext(label, outer=TRUE, cex=1.2)

dev.off()


## metrics


metrics <- data.frame(infile="dummy", period="mon", analysis="trend",
                      slope = 0, stringsAsFactors=FALSE)

for(p in names(infiles)){
  metrics <- rbind(metrics, data.frame(infile = infiles[p],
                                 period = month.abb,
                                 analysis = "trend",
                                 slope = slope[[p]],
                                 row.names=NULL))
}



## write out metrics

metrics <- metrics[-1,]

write.table(metrics, file=txtfile, quote=FALSE, sep="\t",
            row.names=FALSE)

