library(DT)
library(RColorBrewer)
library(colorspace)
library(devtools)
load_all("~/Desktop/climod")

## Call as: Rscript dtable.R analysis indir outname

## analysis = type of analysis (pfit, gev, etc.)
## indir    = directory where metrics files (type.*.txt) are foudn
## outname  = name of output files (out.html & out_files/)

## Note: if you only want some of the files in the directory, just run
## it for everything and filter them out in the datatable.

## Note also: if outname is html/dirname, saving generates an error
## message and saves to html/html/dirname instead.

args <- c("pfit",
          "~/Desktop/fig/raw/all/metrics",
          "test.pfit")

## Comment out this line for testing
#args <- commandArgs(trailingOnly=TRUE)

analysis <- args[1]
indir    <- args[2]
outname  <- args[3]


outfile <- paste0(outname, ".html")
outdir  <- paste0(outname, "_files")
caption <- gsub(".", " ", basename(outname), fixed=TRUE)


## Read in all the analysis.*.txt files in indir and bind them
## together into a big dataframe

infiles <- list.files(indir, paste0(analysis,".*.txt"))

rawdata <- lapply(paste0(indir,"/",infiles), read.table,
                  header=TRUE, stringsAsFactors=FALSE)

rawdf <- do.call(rbind, rawdata)


## The first three columns are filename, period, and analysis.

## Remove duplicate rows based on those columns.  (Typically the obs
## row, which is duplicated across analyses of different models.)
## Have to based redundancy on only the first three (identifier)
## columns, because some analyses (e.g., GEV) don't produce
## numerically identical results from run to run.

dups <- duplicated(rawdf[,1:3])
ddf <- rawdf[-which(dups),]


## Reform the dataframe: drop the analysis column and convert the
## filename column to categorical variables.

## Model filenames look like: dir/var.scen.GCM.RCM.loc.nc
## Obs data filenames look like: dir/var.obs.dataset.loc.nc

fields <- strsplit(basename(ddf$infile), ".", fixed=TRUE)

catlist <- list()

## Obs: drop last field (.nc)
obsind <- grep("obs", ddf$infile)
catlist[obsind] <- lapply(fields[obsind], `[`, -5)

## Model: drop last field, combine RCM & GCM into dataset
modind <- seq(length(fields))[-obsind]
catlist[modind] <- lapply(fields[modind], function(x){
    c(x[1], x[2], paste0(x[4],'.',x[3]), x[5])
})

## Convert into a dataframe
categoricals <- as.data.frame(t(as.matrix(catlist)))
colnames(categoricals) <- c("variable", "scenario", "dataset", "location")


dframe <- cbind(categoricals, ddf[,c(-1,-3)])

## relevel scenario to put obs first in sorting
dframe$scenario <- relevel(dframe$scenario, ref="obs")


## drop nodata data points (used in test suite)
dframe <- dframe[!dframe$location == "nodata",]


## initial sort by scenario, location, and dataset
sdf <- dframe[with(dframe, order(scenario, location, dataset)), ]


## Drop obs for analyses where all metrics are relative to obs
if(analysis == "pfit"){
    sdf <- sdf[!sdf$scenario == "obs",]
}



## Create the datatable
    
## Note: colors need to be either hex (which RColorBrewer and
## colorspace emit) or names for HTML; R color names don't work.


## Basic colorization for categoricals first

html <- datatable(sdf,
                  options = list(paging=FALSE),
                  rownames=FALSE,
                  caption=caption,
                  filter="top"
                  ) %>%
    formatStyle("scenario",
                backgroundColor=styleEqual(
                    c("obs","hist","rcp45","rcp85"),
                    c("#FFFFFF","#E4E4FF","#FFFFE4","#FFE4E4")
                    )) %>%
    formatStyle("dataset",
                backgroundColor=styleEqual(
                    levels(sdf$dataset),
                    brewer.pal(nlevels(sdf$dataset), "Pastel1")
                    )) %>%
    formatStyle("location",
                color=styleEqual(
                    levels(sdf$location),
                    rainbow_hcl(nlevels(sdf$location))
                    ))


## Red to blue bg color style for correlations, etc.
redblue <- styleInterval(seq(-0.95,0.95,0.1),
                         c(hsv(1,   (10:1)/10, 1),
                           "#FFFFFF",
                           hsv(2/3, (1:10)/10, 1)
                           ))

if(analysis == "gev"){
    html <- html %>%
    formatStyle("rlevlo",
              background=styleColorBar(sdf$rlevlo, "#CCCCCC")
              ) %>%
    formatStyle("rlevmid",
              background=styleColorBar(sdf$rlevmid, "#CCCCCC")
              ) %>%
    formatStyle("rlevhi",
              background=styleColorBar(sdf$rlevhi, "#CCCCCC")
              )
}


if(analysis == "pfit"){

    ## Max values of MAD metrics for colorbar backgrounds
    madmax <- lapply(sdf[,c("freqmad","intmad","totmad")], max, na.rm=TRUE)

    html <- html %>%
      formatStyle("freqcor", backgroundColor=redblue) %>%
      formatStyle("intcor",  backgroundColor=redblue) %>%
      formatStyle("totcor",  backgroundColor=redblue) %>%
      formatStyle("freqmad",
                  background=styleColorBar(c(0, madmax$freqmad), "#CCCCCC")
                  ) %>%
      formatStyle("intmad",
                  background=styleColorBar(c(0, madmax$intmad), "#CCCCCC")
                  ) %>%
      formatStyle("totmad",
                  background=styleColorBar(c(0,madmax$totmad), "#CCCCCC")
                  )    
}


## Save to file

htmltools::save_html(html, file=outfile, lib=outdir, background="white")

