#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)

## Call as: Rscript plot.spells.R label obs cur fut png var txt

## png = name of output file for plotted figure
## txt = name of output file for plain-text table of metrics

# for testing
args <- c("spells rcp85 HadGEM2-ES RegCM4 birmingham",
          "big7-obs/prec.METDATA.22i.birmingham.nc",
          "big7-raw/prec.hist.HadGEM2-ES.RegCM4.day.NAM-22i.raw.birmingham.nc",
          "big7-raw/prec.rcp85.HadGEM2-ES.RegCM4.day.NAM-22i.raw.birmingham.nc",
          "test.spells.png",
          "prec",
          "test.spells.txt"
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
stopifnot(v == "prec")

txtfile <- args[7]


## color palette for bias-correction
ocf <- c(obs="black", cur="blue", fut="red")


        
nc <- lapply(infiles, nc_ingest)

## extract variables of interest from the netcdf objects
data <- lapply(nc,"[[",v)
units <- data$cur@units
    
time <- lapply(nc,"[[","time")

time <- lapply(time, alignepochs, "days since 1950-01-01")

## calculate wet/dry spells

wetdry <- lapply(data, function(x){c(x >= 1)})


## calculate spell size and year/day from boolean timeseries x & time
syd <- function(x, t){
    runs <- rle(x)

    size <- runs$length
    type <- factor(runs$values+1, labels=c("dry","wet"))
    
    istop  <- cumsum(runs$length)
    istart <- c(1, head(istop+1,-1))
    tmid <- c((t[istop] + t[istart] ) / 2)

    ylen <- yearlength(t@calendar)
    doy  <- tmid %% ylen
    year <- 1950 + floor(tmid / ylen)

    namelist(doy, year, size, type)
}

spells <- mapply(syd, wetdry, time, SIMPLIFY=FALSE)

ocf <- c("obs","cur","fut")
nspell <- sapply(renest(spells)$size, length)

for(i in ocf){
    spells[[i]]$run <- factor(rep(i, nspell[i]), levels=ocf)
}


#wdmap <- adjustcolor(c(dry="brown", wet="green"), alpha=0.5)
#wdmap <- c(dry=rgb(140,81,10,128,max=255), wet=rgb(1,102,95,128,max=255))
#wdmap <- c(dry=rgb(140,81,10,max=255), wet=rgb(1,102,95,max=255))

wdmap <- c(dry=rgb(191,129,45,max=255), wet=rgb(53,151,143,max=255))

bubble <- function(x, cmap=wdmap, size=0.075, pch=19, ...){
    plot(x$doy, x$year, xlab="day of year", ylab="year", pch=pch,
         col=cmap[x$type], cex=x$size*size, ...)    
}

#bubble(spells$obs, main="obs")


#png(fname, units="in", res=120, width=7, height=11)

## Use panels sized proportionally to number of years in each dataset

nyr <-sapply(lapply(renest(spells)$year, range), diff)

yrat <- c(0, cumsum(nyr) / sum(nyr))

smat <- cbind(rep(0,3), rep(1,3), yrat[1:3], yrat[2:4])





dev.new(width=9, height=12)

split.screen(smat)

## mar: bottom, left, top, right
## mgp: axis title, axis label, axis line  (3, 1, 0)

for(i in 1:3){
    screen(i)
    if(i==1) {par(oma=c(0,0,2,0))}
    par(mar=c(3,3,2,1), mgp=c(1.5,0.5,0))
    bubble(spells[[i]])
    title(names(spells)[i], line=0.5) 
    if(i==1) {mtext(label, line=0, outer=TRUE)}
}
close.screen(all.screens=TRUE)


df <- do.call(rbind, lapply(spells, as.data.frame))

df$week <- pmin(round(df$doy / 7)+1 , 52)
df$month <- pmax(round(df$doy / 30.5 + 0.5), 1)

par(mfrow=c(3,2))

for(i in levels(df$run)){
    for(j in levels(df$type)){
        boxplot(size ~ week, data=subset(df, run == i & type == j),
                main=paste(i,j), xlab="week", ylim=c(1,max(df$size)),
                log="y")
    }
}



# ####################################
#     png(fname, units="in", res=120, width=7, height=7)
#     par(mfrow=c(3,1), oma=c(0,0,3,0), mar=c(4,4,3,2), mgp=c(2.5,1,0))
# 
# 
#     mtext(title, line=1, outer=TRUE)
# 
#     ## legends in top corners outside plots
#     par(fig=c(0,1,0,1), mfrow=c(1,1), mar=c(0,0,0,0), oma=rep(0.5,4), new=TRUE)
#     plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n", ann=FALSE)    
#     legend("topright", names(ocf), col=ocf, lty=1,
#            horiz=TRUE, cex=0.65, seg.len=1)
#     legend("topleft", c("max", "min"), col="gray", lty=c(1,2),
#            horiz=TRUE, cex=0.65)
# 
#     dev.off()
# 
# 
# 
# plot(prec, outfile, label)
# 
# 
# 
# ## metrics
# 
# mprec <- lapply(prec, function(x){as.list(as.data.frame(x))})
# 
# mcor <- lapply(mprec, function(x){lapply(x, function(y){cor(y, x$obs)})})
# mcor <- lapply(renest(mcor), unlist)
# mmad <- lapply(mprec, function(x){lapply(x, function(y){mad(y-x$obs, center=0)})})
# mmad <- lapply(renest(mmad), unlist)
# 
# 
# metrics <- data.frame(infile="dummy", period="ann", analysis="pfit",
#                       freqcor=0, intcor=0, totcor=0,
#                       freqmad=0, intmad=0, totmad=0,
#                       stringsAsFactors=FALSE)
# 
# for(p in names(infiles)){
# 
#   m0 <- list(infile=infiles[p], period="ann", analysis="pfit")
#  
#   m1 <- rev(as.list(mcor[[p]]))
#   names(m1) <- paste0(names(m1),"cor")
# 
#   m2 <- rev(as.list(mmad[[p]]))
#   names(m2) <- paste0(names(m2),"mad")
# 
#   metrics <- rbind(metrics, c(m0, m1, m2))
# }
# 
# 
# ## write out metrics
# 
# metrics <- metrics[-1,]
# 
# write.table(format(metrics, trim=TRUE, digits=3),
#             file=txtfile, quote=FALSE, sep="\t", row.names=FALSE)
# 
# 
