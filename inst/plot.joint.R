#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)
library(MASS)
suppressMessages(library(mgcv))


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


## Call as: Rscript plot.joint.R label obs cur fut png var txt

## provide input files for prec; tmax & tmin are assumed analogous
## png = name of output file for plotted figure
## txt = name of output file for plain-text table of metrics

args <- c("rcp85 HadGEM2-ES WRF birmingham",
          "obs/prec.obs.livneh.birmingham.nc",
          "save-test/prec.hist.HadGEM2-ES.WRF.birmingham.nc",
          "save-test/prec.rcp85.HadGEM2-ES.WRF.birmingham.nc",
          "test.joint.png",
          "prec",
          "test.joint.txt"
          )

## comment out this line for testing
#args <- commandArgs(trailingOnly=TRUE)

label <- args[1]

infiles <- list()
infiles$prec <- c(obs=args[2], cur=args[3], fut=args[4])
infiles$tmax <- gsub("prec", "tmax", infiles$prec)
infiles$tmin <- gsub("prec", "tmin", infiles$prec)

outfile <- args[5]
txtfile <- args[7]

nc <- lapply(infiles, function(x){lapply(x, nc_ingest)})



## extract variables of interest from the netcdf objects
data <- list()
for(v in names(nc)){
  data[[v]] <- lapply(nc[[v]], `[[`, v)
}


punits <- data$prec$obs@units
tunits <- data$tmax$obs@units


time <- list()
for(v in names(nc)){
  time[[v]] <- lapply(nc[[v]], `[[`, "time")
}


atime <- rapply(time, alignepochs, epoch="days since 1949-12-01", how="replace")
## Aligning on 12-01 instead of 01-01 makes the cslices be
## DJF/MAM/JJA/SON instead of JFM/AMJ/JAS/OND.


## Subset data to common times across variables

rtime <- renest(atime)
rdata <- renest(data)

##  We need a 1:1 mapping between precip & temp.  The later stages
##  deal with differences in start and end times, missing data int he
##  middle, etc.  However, it depends on the time corodinates of the
##  three variables being otherwise the same.  The code chunk below
##  deals with the special case of one variable having its coordinates
##  at a different point in the interval (e.g., precip is at 12:00
##  while temp is at 00:00).  The right place to fix this is in the
##  datafiles, but that's unfinished while this code is being
##  developed, so we fix it manually and emit a warning about it.

htime <- rtime
for(p in names(htime)){
  offset <- htime[[p]]$prec[1] %% 1
  for(v in c("tmax","tmin")){
    adjust <- offset - htime[[p]][[v]][1] %% 1
    if(adjust != 0){
      warning(paste("Manually adjusting", p, v,
                    "time to match prec time!"))}
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
  ddata[[p]]$tmrg <- (ddata[[p]]$tmax + ddata[[p]]$tmin)/2
}


## generate seasonal climatology window index arrays
cwin <- lapply(renest(ttime)$prec, cslice, ratio=1, num=4, split=FALSE)


## slice data
sdata <- mapply(function(d,s){lapply(d, slice, s)}, ddata, cwin, SIMPLIFY=FALSE)


## subset
pdata <- lapply(sdata, `[[`, "prec")
tdata <- lapply(sdata, `[[`, "tmrg")


## Log precip
punits <- paste0("log(",punits,")")
suppressWarnings(pdata <- rapply(pdata, log10, how="replace"))


## drop non-finite values & sub-trace precip

ptrace <- -4   # log scale; threshold based on gridded obs range of values
 
pgood <- rapply(pdata, how="replace", function(x){is.finite(x) & x > ptrace})
tgood <- rapply(tdata, how="replace", function(x){is.finite(x)})
isgood <- mapply(function(x,y){mapply(`&`, x, y, SIMPLIFY=FALSE)}, pgood, tgood, SIMPLIFY=FALSE)


psub <- rsub(pdata, isgood)
tsub <- rsub(tdata, isgood)


## combine, push p & t down to bottom of nesting & pull seas up to top
jsdata <- renest(lapply(renest(list(temp=tsub, prec=psub)), renest))
names(jsdata) <- seas


## Function: gets polygons for contours of field as 2-column array
## contourLines

conpolyarr <- function(g2d, contourlevel){
  
  pout <- contourLines(x=g2d$x, y=g2d$y, z=g2d$z, levels=contourlevel)
  
  ## Need to convert lists to 2-column arrays.  contourLines returns a
  ## list because there could be multiple polygons.  If there's more
  ## than 1 polygon, we need to merge them into one big array and
  ## stick rows of NA between them.
  
  la <- lapply(pout, function(p){cbind(p$x, p$y)})
  
  a <- la[[1]]
  if(length(la) > 1){
    for(i in 2:length(la)){
      a <- rbind(a, NA, la[[i]])
    }
  }
  return(a)
}


## Empty metrics dataframe

metrics <- data.frame(infile = "dummy", period="seas", analysis="joint",
                      ktau=0, q50over=0, modeP=0, modeT=0,
                      stringsAsFactors=FALSE)

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
    imode <- which(gf$z == max(sz), arr.ind=TRUE)
    modept <- list(x=gf$x[imode[1]], y=gf$y[imode[2]])
    points(modept$x, modept$y, pch=19, col=cmap[i], cex=1.5)
    ## TODO: change mode point to something like 2D Tukey median
    
  
    ## plot individual points outside 90% contour
    ## sp::point.in.polygon has nicer syntax, but crashes R

    a90 <- conpolyarr(gf, clev[1])
      
    outpts <- !in.out(a90, cbind(d$prec, d$temp))
    points(d$prec[outpts], d$temp[outpts], col=cmap[i], pch=19, cex=0.2)

    
    
    ## Calculate metrics

    ## correlation between P & T  (Kendall's Tau)
    ktau <- cor(d$prec, d$temp, method="kendall")

    ## % overlap between obs & mod q50 contours
    if (i == "obs"){
      ## As long as obs is first, this will carry over to later loops      
      grid <- as.matrix(expand.grid(gf$x, gf$y, KEEP.OUT.ATTRS=FALSE))
      oqa50 <- conpolyarr(gf, clev[2])      
      oq50io <- in.out(oqa50, grid)
    }

    qa50 <- conpolyarr(gf, clev[2])
    q50io <- in.out(qa50, grid)
    q50over <- sum(q50io & oq50io) / sum(q50io)

    ## Add to metrics df
    metrics <- rbind(metrics,
                     list(infiles$prec[i], s, "joint",
                          ktau, q50over,
                          modept$x, modept$y))
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


## write out metrics

metrics <- metrics[-1,]

write.table(format(metrics, trim=TRUE, digits=3),
            file=txtfile, quote=FALSE, sep="\t", row.names=FALSE)
