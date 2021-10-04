#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)

## Call as: Rscript plot.mpdf.R label obs cur fut png var txt

## png = name of output file for plotted figure
## txt = name of output file for plain-text table of metrics

# for testing
args <- c("prec rcp85 HadGEM2-ES RegCM4 birmingham",
          "obs/prec.obs.livneh.birmingham.nc",
          "save-test/prec.hist.HadGEM2-ES.RegCM4.birmingham.nc",
          "save-test/prec.rcp85.HadGEM2-ES.RegCM4.birmingham.nc",
          "test.mpdf.png",
          "prec",
          "test.mpdf.txt"
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
data <- lapply(nc,"[[", v)

units <- data$cur@units
    
time <- lapply(nc,"[[","time")
time <- lapply(time, alignepochs, "days since 1950-01-01")


## Log precip
if(v == "prec"){
  units <- paste0("log(",units,")")
  suppressWarnings(data <- lapply(data, log10))
}


## generate monthly climatology window index arrays
dwin <- lapply(time, cslice, ratio=1, num=12, split=FALSE)

## slice data
cdata <- mapply(slice, data, dwin, SIMPLIFY=FALSE)


## drop non-finite & sub-trace values
trace <- -7   # log scale
cdata <- rapply(cdata, function(x){x[is.finite(x) & x>trace]}, how="replace")

## rearrange
mdata <- renest(cdata)

## color palette for bias-correction
ocf <- c(obs="black", cur="blue", fut="red")

## plotting limits
extend <- 3*stats::bw.nrd(unlist(mdata))
    xr <- extend*c(-1,1)+range(unlist(mdata))

## plotting

png(outfile, units="in", res=120, width=10, height=7)

par(mfrow=c(3,4), oma=c(0,0,3,0), mar=c(2.5,1.5,4.5,1.5)+0.1)

mdata <- renest(cdata)
names(mdata) <- month.abb


skill <- list(obs=list(pdf=c(), tail=c()), cur=list(), fut=list())
for(p in names(ocf)){
  skill[[p]]$pdf  <- mapply( pdfskill, cdata$obs, cdata[[p]])
  skill[[p]]$tail <- mapply(tailskill, cdata$obs, cdata[[p]])
}
skill <- rapply(skill, how="replace", function(x){names(x)<- month.abb;x})
skill <- rapply(skill, how="replace", format, trim=TRUE, digits=3)


for(m in month.abb){
  mplot(lapply(mdata[[m]], akde, min.x=min(xr), max.x=max(xr)),
        yaxt="n", col=ocf, type='l', lty=1, main=m, xlab=units, ylab="")

  ## Skill annotations above plot
  mtext("pdf skill", adj=0, line=2.5, cex=0.7)
  mtext(skill$cur$pdf[m], adj=0, line=1.5, col=ocf[2], cex=0.7)
  mtext(skill$fut$pdf[m], adj=0, line=0.5, col=ocf[3], cex=0.7)

  mtext("tail skill", adj=1, line=2.5, cex=0.7)
  mtext(skill$cur$tail[m], adj=1, line=1.5, col=ocf[2], cex=0.7)
  mtext(skill$fut$tail[m], adj=1, line=0.5, col=ocf[3], cex=0.7)
}

mtext(paste(label, "PDFs (",units,")"), line=1, outer=TRUE)


## upper-margin legend outside plots

par(fig=c(0,1,0,1), oma=c(0,0,1,1), mar=c(0,0,0,0), new=TRUE)
plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n", ann=FALSE)
legend("topright", names(ocf), col=ocf, lty=1, seg.len=1, horiz=TRUE)


dev.off()


## metrics

metrics <- data.frame(infile="dummy", period="mon", analysis="mpdf",
                      pdfsk = 0, tailsk = 0, stringsAsFactors=FALSE)

for(p in names(infiles)){
  metrics <- rbind(metrics, data.frame(infile = infiles[p],
                                 period = month.abb,
                                 analysis = "mpdf",
                                 pdfsk = skill[[p]]$pdf,
                                 tailsk = skill[[p]]$tail,
                                 row.names=NULL))
}


## write out metrics

metrics <- metrics[-1,]

write.table(metrics,
            file=txtfile, quote=FALSE, sep="\t", row.names=FALSE)
