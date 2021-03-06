#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
suppressMessages(library(quantreg))

## Call as: Rscript plot.quants.R label rsave out var
## rsave = name of saved Rdata file with peak slices

## for testing
args <- c("tmin HadGEM2-ES RegCM4 ftlogan",
          "v1.lof.bc/dmaps/tmin.rcp85.HadGEM2-ES.RegCM4.ftlogan.Rdata",
          "test.quants.png",
          "tmin")

## comment out this line for testing
args <- commandArgs(trailingOnly=TRUE)


label <- args[1]

load(args[2])

outfile <- args[3]

v <- args[4]


cmap <- c(obs="black", cur="blue", fut="red")
lmap <- c(3,2,5,1,5,2,3)

## Quantiles - equally spaced if data is normal
tau <- c(0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98)


png(outfile, units="in", res=120, width=7, height=7)

par(mfrow=c(3,2), oma=c(0,0,3,0), mar=c(3,3,2,1))

for (p in names(peaks)){

  pyear <- lapply(tslice[[p]], unlist)
  pfix  <- lapply(outerdata$fix[[p]], unlist)
  praw  <- lapply(outerdata$raw[[p]], unlist)

  df <- mapply(data.frame, pyear, pfix, praw, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  df <- lapply(df, setNames, c("year","fix","raw"))
  names(df) <- names(pyear)

  fraw <- as.formula("raw ~ year")
  ffix <- as.formula("fix ~ year")

  qfix <- lapply(df, rq, formula=ffix, tau=tau)
  qraw <- lapply(df, rq, formula=fraw, tau=tau)
    
  mplot(df, x="year", y="fix", type="n", xlab="", ylab="")

  for(i in names(qfix)){
    matplot(pyear[[i]], predict(qfix[[i]]),
            type="l", col=cmap[i], add=TRUE, lty=lmap)
  }

  title(paste0(p ,": ", dates[p]), line=1)

}
  
mtext(paste(label,"quantile regression (",units,")"), line=1, outer=TRUE)

dev.off()

