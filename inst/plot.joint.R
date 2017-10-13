#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)
library(MASS)
library(mgcv)
#library(RColorBrewer)
#library(colorspace)

seas <- c("DJF","MAM","JJA","SON")
cmap <- c(obs="black", cur="blue", fut="red")

## Recursive subset x[y] on identically structured nested lists
rsub <- function(x,y){
  if(is.atomic(x) && is.atomic(y)){
    copyatts(x, x[y])
 } else {
   Map(rsub, x, y)
 }
}


## Call as: Rscript plot.joint.R label obs cur fut out

args <- commandArgs(trailingOnly=TRUE)

# # for testing
# args <- c("rcp85 HadGEM2-ES RegCM4 birmingham",
#           "obs/prec.obs.livneh.birmingham.nc",
#           "save-test/prec.hist.HadGEM2-ES.RegCM4.birmingham.nc",
#           "save-test/prec.rcp85.HadGEM2-ES.RegCM4.birmingham.nc",
#           "test.joint.png"
#           )


label <- args[1]

infiles <- list()
infiles$prec <- c(obs=args[2], cur=args[3], fut=args[4])
infiles$tmax <- gsub("prec", "tmax", infiles$prec)
infiles$tmin <- gsub("prec", "tmin", infiles$prec)

# precfiles <- c()
# precfiles["obs"] <- args[2]
# precfiles["cur"] <- args[3]
# precfiles["fut"] <- args[4]
# 
# tmaxfiles <- gsub("prec","tmax",precfiles)
# tminfiles <- gsub("prec","tmin",precfiles)
# infiles <- c(precfiles, tmaxfiles, tminfiles)
# names(infiles) <- outer(names(precfiles),c("prec","tmax","tmin"),paste0)

outfile <- args[5]

#nc <- lapply(infiles, nc_ingest)

nc <- lapply(infiles, function(x){lapply(x, nc_ingest)})



## extract variables of interest from the netcdf objects
#v <- rep(c("prec","tmax","tmin"), each=3)
#data <- mapply(`[[`, nc, v, SIMPLIFY=FALSE)

data <- list()
for(v in names(nc)){
  data[[v]] <- lapply(nc[[v]], `[[`, v)
}


#punits <- data$obsprec@units
#tunits <- data$obstmax@units

punits <- data$prec$obs@units
tunits <- data$tmax$obs@units


#time <- lapply(nc,"[[","time")

time <- list()
for(v in names(nc)){
  time[[v]] <- lapply(nc[[v]], `[[`, "time")
}

#time <- lapply(time, alignepochs, "days since 1949-12-01")

atime <- rapply(time, alignepochs, epoch="days since 1949-12-01", how="replace")


## Aligning on 12-01 instead of 01-01 makes the cslices be
## DJF/MAM/JJA/SON instead of JFM/AMJ/JAS/OND.


## trim data to common coverage periods

rtime <- renest(atime)
rdata <- renest(data)

####  The basic problem that I am in the middle of fixing here: we
####  need a 1:1 mapping between precip & temp.  So for obs/cur/fut,
####  the time coordinates of the three variables need to match.  It
####  may commonly be the case that the endpoints will not match, in
####  which case we can trim things to the minimal coverage period.
####  But if the coordinates are off (e.g., precip is on 0.5, while
####  temp is on 0, as is the case with the WRF data currently), then
####  I should be fixing it in the data, rather than hackily fixing it
####  here.  So this next little bit is a hack while I wait for the
####  data to be fixed properly.

htime <- rtime
for(p in names(htime)){
  offset <- htime[[p]]$prec[1] %% 1
  for(v in c("tmax","tmin")){
    adjust <- offset - htime[[p]][[v]][1] %% 1
    if(adjust != 0){ warning(paste("Manually adjusting", p, v, "time to match prec time!"))}
    htime[[p]][[v]] <- htime[[p]][[v]] + adjust
  }
}


## Data could be missing in the middle (WRF), so we use %in% rather
## than just checking the beginning and ending.

matching <- lapply(htime, function(x){
  list(prec = ((x$prec %in% x$tmax) & (x$prec %in% x$tmin)),
       tmax = ((x$tmax %in% x$tmin) & (x$tmax %in% x$prec)),
       tmin = ((x$tmin %in% x$prec) & (x$tmin %in% x$tmax)))})


ttime <- rsub(htime, matching)
ddata <- rsub(rdata, matching)



## calculate tmrg = T_midrange = (Tmax + Tmin)/2
for(p in names(ddata)){
  stopifnot(identical(ttime[[p]]$tmax, ttime[[p]]$tmin))
#  ttime[[p]]$tmrg <- ttime[[p]]$tmax
  ddata[[p]]$tmrg <- (ddata[[p]]$tmax + ddata[[p]]$tmin)/2
}



#for(i in c("obs","cur","fut")){
#  ti <- paste0(i, c("tmin","tmax","tmrg"))  
#  stopifnot(identical(time[[ti[1]]], time[[ti[2]]]))
# time[[ti[3]]] <- time[[ti[1]]]
#  data[[ti[3]]] <- (data[[ti[1]]] + data[[ti[2]]])/2
#}



## generate seasonal climatology window index arrays
#cwin <- lapply(time, cslice, ratio=1, num=4, split=FALSE)
cwin <- lapply(renest(ttime)$prec, cslice, ratio=1, num=4, split=FALSE)

## slice data
#sdata <- mapply(slice, data, cwin, SIMPLIFY=FALSE)
sdata <- mapply(function(d,s){lapply(d, slice, s)}, ddata, cwin, SIMPLIFY=FALSE)


## subset
#pdata <- sdata[grep("prec", names(sdata))]
#tdata <- sdata[grep("tmrg", names(sdata))]
#
#names(pdata) <- names(cmap)
#names(tdata) <- names(cmap)

pdata <- lapply(sdata, `[[`, "prec")
tdata <- lapply(sdata, `[[`, "tmrg")


## Log precip
punits <- paste0("log(",punits,")")
suppressWarnings(pdata <- rapply(pdata, log10, how="replace"))


## drop non-finite values & sub-trace precip

## This doesn't work with WRF because current test-suite prec & temp
## have different length (plus other time coordinate issues)

## Needed: once times are aligned, trim based on max(min(time)) and
## min(max(time)) to get common coverage periods.


ptrace <- -4   # log scale; threshold based on gridded obs range of values
 
pgood <- rapply(pdata, how="replace", function(x){is.finite(x) & x > ptrace})
tgood <- rapply(tdata, how="replace", function(x){is.finite(x)})
isgood <- mapply(function(x,y){mapply(`&`, x, y, SIMPLIFY=FALSE)}, pgood, tgood, SIMPLIFY=FALSE)


psub <- rsub(pdata, isgood)
tsub <- rsub(tdata, isgood)


## combine, push p & t down to bottom of nesting & pull seas up to top
jsdata <- renest(lapply(renest(list(temp=tsub, prec=psub)), renest))
names(jsdata) <- seas


## plotting

png(outfile, units="in", res=120, width=10, height=7)

par(mfrow=c(2,2), oma=c(0,0,3,0), mgp=c(2,1,0), mar=c(3.5,3.5,3,1.5))

for (s in seas){

  sd <- jsdata[[s]]

  mplot(sd, type="n", x="prec", y="temp", main=s, xlab=punits, ylab=tunits)

  ## freezing line
  abline(h=0, col="gray")

  for(i in names(sd)){
    d <- sd[[i]]
  
    ## gridded 2D density estimate
    gf <- kde2d(d$prec, d$temp, lims=par("usr"))
  
    ## find contour levels containing 95% and 50% of total 2D-PDF AUC
    ## (Ordering outer to inner makes plotting outliers easier.)  
    plev <- c(0.95, 0.5)
  
    sz <- sort(gf$z)
    auc <- sum(sz)
    clev <- sz[findInterval((1-plev)*auc, cumsum(sz))]
  
    ## plot contours
    contour(gf, levels=clev, drawlabels=FALSE,
          lwd=seq(length(plev)), col=cmap[i], add=TRUE)
    
    ## add mode point
    modept <- which(gf$z == max(sz), arr.ind=TRUE)
    points(gf$x[modept[1]], gf$y[modept[2]], pch=19, col=cmap[i], cex=1.5)
  
  
    ## plot individual points outside 90% contour
    ## sp::point.in.polygon has nicer syntax, but crashes R
  
    pout <- contourLines(x=gf$x, y=gf$y, z=gf$z, levels=clev[1])
  
    ## Need to convert lists to 2-column arrays.  contourLines returns a
    ## list because there could be multiple polygons.  If there's more
    ## than 1 polygon, we need to merge them into one big array and
    ## stick rows of NA between them.
  
    la90 <- lapply(pout, function(p){cbind(p$x, p$y)})
  
    a90 <- la90[[1]]
    if(length(la90) > 1){
      for(i in 2:length(la90)){
        a90 <- rbind(a90, NA, la90[[i]])
      }
    }
      
    outpts <- !in.out(a90, cbind(d$prec, d$temp))
    points(d$prec[outpts], d$temp[outpts], col=cmap[i], pch=19, cex=0.2)
  
  }

}


mtext(paste(label, "joint PDFs"), line=1, outer=TRUE)

## centered legend outside plots

par(fig=c(0,1,0,1), new=TRUE)
plot(c(0,1), c(0,1), type="n", bty="n", xaxt="n", yaxt="n", ann=FALSE)

nlev <- length(plev)
legend(0.5, 0.5, ncol=2, bty="n", seg.len=1, x.intersp=0.5, xjust=0.5, yjust=0.65,
       c(names(cmap), paste0(plev*100,"%"), "mode"),
       fill=c(cmap, rep(NA, nlev+1)), border=NA,
       lwd=c(NA,NA,NA,seq(nlev),NA),
       pch=c(rep(NA,nlev+3),19)
       )
 
dev.off()

