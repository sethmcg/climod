#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
#library(ncdf4)
library(quantreg)

## Call as: Rscript --vanilla box.R label obs cur fut out varname slice

args <- commandArgs(trailingOnly=TRUE)

#args <- c("tmax HadGEM2-ES RegCM4 ftlogan",
#          "tmax.test.Rdata",
#          "test.quants.png",
#          "tmax")

label <- args[1]

load(args[2])

outfile <- args[3]

v <- args[4]


cmap <- c(obs="black", cur="blue", fut="red")
lmap <- c(3,2,5,1,5,2,3)

## Quantiles - equally spaced if data is normal
tau <- c(0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98)


png(outfile, units="in", res=120, width=7, height=7)

par(mfrow=c(4,1), oma=c(0,0,2,0), mar=c(3,3,2,1))

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
  
mtext(paste(label,"quantile regression (",units,")"), outer=TRUE)

dev.off()

