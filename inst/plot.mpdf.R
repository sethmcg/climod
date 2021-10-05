#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)

## Call as: Rscript plot.mpdf.R label obs cur fut png var txt

## png = name of output file for plotted figure
## txt = name of output file for plain-text table of metrics

# for testing
args <- c("tmin rcp85 MPI-ESM-MR CRCM5-UQAM NAM-44i mbcn-gridMET",
          "test/tmin.obs.gridMET.day.NAM-44i.nc",
          "test/tmin.hist.MPI-ESM-MR.CRCM5-UQAM.day.NAM-44i.mbcn-gridMET.nc",
          "test/tmin.rcp85.MPI-ESM-MR.CRCM5-UQAM.day.NAM-44i.mbcn-gridMET.nc",
          "figs/test.mpdf.png",
          "tmin",
          "figs/test.mpdf.txt"
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


## extract time & make sure everything is aligned

time <- lapply(nc,"[[","time")
epoch <- "days since 1900-01-01"
time <- lapply(time, alignepochs, epoch)


## concatenate model current & future

tcat <- function(x){
    x$mod <- copyatts(x$fut, c(x$cur, x$fut))
    x$cur <- x$fut <- NULL
    return(x)
}

Data <- tcat(Data)
time <- tcat(time)


## Calculate month & year
## (Not proper months, just 1/12 of the year length)

ym <- function(T){
    ylen   <- yearlength(T)
    fyear  <- T / ylen + 1900
    year   <- c(floor(fyear))
    month  <- ceiling(fyear %% 1 * 12 +1e-10) ## 1e-10 to avoid month = 0
    data.frame(time=fyear, year, month)
}

df <- mapply(cbind, lapply(time, ym), data=Data, SIMPLIFY=FALSE)


## Split model data into 30-year periods, plus a current period
## matching the obs period.

period <- list(obs  = c(1979, 2017),
               cur  = c(1979, 2017),
               hist = c(1950, 1980),
               turn = c(1980, 2010),
               near = c(2010, 2040),
               mid  = c(2040, 2070),
               far  = c(2070, 2100)
               )

df$cur <- df$mod[df$mod$time >= 1979 & df$mod$time < 2017,]

df <- c(df, split(df$mod, cut(df$mod$time, seq(1950,2100,by=30), right=FALSE)))
names(df)[1:5+3] <- tail(names(period), 5)

df$mod <- NULL


## Split all data frames by month

mdf <- renest(lapply(df, function(x){
    y <- split(x, x$month)
    names(y) <- month.abb
    y})
)


md <- lapply(mdf, function(x){lapply(x, `[[`, "data")})


## find common x-range
extend <- 3*stats::bw.nrd(unlist(md))
    xr <- extend*c(-1,1)+range(unlist(md))


## calculate distributions
mpdf <- lapply(md, function(ld){
    lapply(ld, akde, min=min(xr), max=max(xr))})


## find common y-range for pdfs
yr <- c(0, max(unlist(renest(lapply(mpdf, renest))$y)))

## calculate skill scores
skill <- list()
for(s in c("pdf","tail")){
    skill[[s]] <- list()
    for(m in month.abb){
        skill[[s]][[m]] <- sapply(md[[m]][-1], paste0(s,"skill"), md[[m]]$obs)
    }
}

## skill scores reformatted as matrices for plotting
skillmat <- lapply(skill, do.call, what="rbind")


## calculate monthly means & error in mean
mu <- lapply(renest(rapply(md, mean, how="replace")),unlist)
delta <- as.matrix(mu)[,-1] - mu$obs


## And plot!

yeartext <- paste0('(',sapply(period, paste, collapse='-'),')')
legendtext <- paste(names(period), yeartext)
                    

cmap <- c(obs  = "gray",
          cur  = "black", 
          hist = "darkgreen",
          turn = "blue",
          near = "purple",
          mid  = "red",
          far  = "orange")

lmap <- c(4,rep(2,6))
names(lmap) <- names(period)



## plotting

png(outfile, units="in", res=120, width=10, height=9)
#dev.new(width=10,height=9)

par(mfrow=c(4,4), oma=c(0,0,3,0), mar=c(2,2,2,1)+0.6)

ty <- max(yr) * (1 - 1:6/9)

for(m in month.abb){
    mplot(mpdf[[m]], type='l', col=cmap, lwd=lmap, lty=1,
          xlim=xr, ylim=yr, main=m, xlab=units, ylab='')

    text(min(xr), ty, sprintf("%4.3f", skill$pdf[[m]]),
         pos=4, col=cmap[-1])

    text(max(xr), ty, sprintf("%4.3f", skill$tail[[m]]),
         pos=2, col=cmap[-1])

    mtext("pdf skill",  adj=0, cex=0.7)
    mtext("tail skill", adj=1, cex=0.7)
}

matplot(skillmat$pdf, main="PDF skill", xaxt='n', ylim=c(0,1),
        type='b', col=cmap[-1], pch=1, lty=1)
axis(side=1, at=1:12, labels=substr(month.abb,1,1), cex.axis=0.9)

matplot(skillmat$tail, main="tail skill", xaxt='n', ylim=c(0,1),
        type='b', col=cmap[-1], pch=1, lty=1)
axis(side=1, at=1:12, labels=substr(month.abb,1,1), cex.axis=0.9)

syr <- c(-1,1)*max(abs(delta))
matplot(delta, main="mean delta from obs", xaxt='n', ylim=syr,
        type='b', col=cmap[-1], pch=1, lty=1)
abline(h=0, col='gray')
axis(side=1, at=1:12, labels=substr(month.abb,1,1), cex.axis=0.9)

#mtext(paste(label, "PDFs (",units,")"), line=3, outer=TRUE)
mtext(paste(label, "PDFs (",units,")"), outer=TRUE)


## legend

np <- length(period)

plot(0,0, type='n', xaxt='n', yaxt='n', bty='n', ylab='', xlab='')

legend("center", c(names(period), yeartext), title="Period", ncol=2,
       lwd=lmap, bty='n', col=c(cmap,rep(NA,np)),
       seg.len=rep(c(2,-4), each=np))


dev.off()

## metrics

N <- length(delta)
metrics <- cbind(rep(infiles["fut"],N),
                 rep("mpdf",N),
                 as.data.frame.table(skillmat$pdf,
                 stringsAsFactors=FALSE),
                 as.data.frame.table(skillmat$tail)[,3],
                 as.data.frame.table(delta)[,3]
                 )

colnames(metrics) <- c("infile","analysis","period", "era", "pdfsk", "tailsk", "delta") 
metrics <- metrics[,c(1,3,2,4:ncol(metrics))]


## write out metrics

write.table(metrics,
            file=txtfile, quote=FALSE, sep="\t", row.names=FALSE)
