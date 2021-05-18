library(DT)
library(RColorBrewer)
library(colorspace)
library(devtools)
load_all("~/Desktop/climod")

## Call as: Rscript dtable.R analysis indir outname

## analysis = type of analysis (pfit, gev, etc.)
## indir    = directory where metrics files (type.*.txt) are found
## outname  = name of output files (out.html & out_files/)

## Note: if you only want some of the files in the directory, just run
## it for everything and filter them out in the datatable.

## Note also: if outname is html/dirname, saving generates an error
## message and saves to html/html/dirname instead.

#args <- c("pfit",
#          "~/Desktop/fig/raw2/all/metrics",
#          "test.pfit")
args <- c("mpdf",
          "~/Desktop/DCA/upper/lowbar/mpdf/metrics",
          "upper-atm.SGP.mpdf",
          "variable.level.scenario.GCM.RCM.drop.years.drop.lon.lat.drop",
          FALSE,
          TRUE)

## Comment out this line for testing
#args <- commandArgs(trailingOnly=TRUE)

analysis <- args[1]
indir    <- args[2]
outname  <- args[3]


## pattern for spitting filename to columns.  Components are used as
## colum names, 'drop' columns are removed.

if(length(args)<4){
    defaultpattern <- TRUE
    colpattern <- "variable.scenario.GCM.RCM.location.drop"
} else {
    defaultpattern <- FALSE
    colpattern <- args[4]
}


## If true, then obs filenames have a dataset component instead of
## RCM/GCM analogs; RCM & GCM are merged to dataset to match.

if(length(args)<5){
    merge <- TRUE
} else {
    merge <- args[5]
}


## Read in all the analysis.*.txt files in indir and bind them
## together into a big dataframe

infiles <- list.files(indir, paste0(analysis,".*.txt"))

rawdata <- lapply(paste0(indir,"/",infiles), read.table,
                  header=TRUE, stringsAsFactors=FALSE)

rawdf <- do.call(rbind, rawdata)


## The first three columns are filename, period, and analysis.

## Remove duplicate rows based on those columns.  (Typically the obs
## row, which is duplicated across analyses of different models.)
## Have to base redundancy on only the first three (identifier)
## columns, because some analyses (e.g., GEV) don't produce
## numerically identical results from run to run.

ddf <- rawdf[!duplicated(rawdf[,1:3]),]

## Reform the dataframe: drop the analysis column and convert the
## filename column to categorical variables.


## extract filename components
fields <- strsplit(basename(ddf$infile), ".", fixed=TRUE)

cnames <- unlist(strsplit(basename(colpattern), ".", fixed=TRUE))

keep <- !grepl("drop", cnames)

categoricals <- as.data.frame(t(as.matrix(fields)),stringsAsFactors=TRUE)[,keep]
colnames(categoricals) <- cnames[keep]

## relevel scenario to put obs first in sorting
categoricals$scenario <- relevel(categoricals$scenario, ref="obs")


# ### OLD:
# 
# ## Model filenames look like: dir/var.scen.GCM.RCM.loc.nc
# ## Obs data filenames look like: dir/var.obs.dataset.loc.nc
# 
# fields <- strsplit(basename(ddf$infile), ".", fixed=TRUE)
# 
# ## Obs: drop last field (.nc)
# obsfields <- lapply(fields, `[`, -5)
# 
# ## Model: drop last field, combine RCM & GCM into dataset
# modfields <- lapply(fields, function(x){c(x[1], x[2], paste0(x[4],'.',x[3]), x[5])})
# 
# 
# ## Use grepl instead of grep (returns logical vector instead of
# ## indices) her because there may be no obs at all, in which case grep
# ## would return 0, and array[-0] = array[0] = nothing.
# obsmask <- grepl("obs", ddf$infile)
# modmask <- !obsmask
# 
# catlist <- list()
# catlist[obsmask] <- obsfields[obsmask]
# catlist[modmask] <- modfields[modmask]
# 
# 
# ## Convert into a dataframe
# categoricals <- as.data.frame(t(as.matrix(catlist)))
# colnames(categoricals) <- c("variable", "scenario", "dataset", "location")


dframe <- cbind(categoricals, ddf[,c(-1,-3)])

## drop nodata data points (used in test suite)
if("location" %in% colnames(dframe)){
    dframe <- dframe[!dframe$location == "nodata",]
}

## get month-based periods to sort correctly
if(all(dframe$period %in% month.abb)){
    dframe$period <- factor(dframe$period, levels=month.abb)
}


## Drop obs for analyses where all metrics are relative to obs
if(analysis %in% c("mpdf", "pfit", "tmmd")){
    dframe <- dframe[dframe$scenario != "obs",]
}


# ## initial sort by categoricals
# #sdf <- dframe[with(dframe, order(variable, scenario, location, dataset)), ]
# sdf <- dframe[with(dframe, order(c(colnames(categoricals),"period"))),]


#for (v in c("variable", "scenario", "dataset", "location", "period")){
#    sdf[[v]] <- factor(sdf[[v]])
#}

## Create the datatable
    
## Note: colors need to be either hex (which RColorBrewer and
## colorspace emit) or names for HTML; R color names don't work.


## sdf has 56640 rows, which is apparently too many for DataTable, so
## we have to split this up by month.

for (m in month.abb){


    outfile <- paste0(outname, ".", m, ".html")
    outdir  <- paste0(outname, ".", m, "_files")
    caption <- paste0(gsub(".", " ", basename(outname), fixed=TRUE), " ", m)

    
    sdf <- dframe[dframe$period == m,]


    ## Basic colorization for categoricals first

    html <- datatable(sdf,
                      options = list(paging=FALSE),
                      rownames=FALSE,
                      caption=caption,
                      filter="top"
                      )

    ## scenario has special custom coloring
    
    if("scenario" %in% colnames(sdf)){
        html <- html %>%
            formatStyle("scenario",
                        backgroundColor=styleEqual(
                            c("obs","hist","rcp45","rcp85"),
                            c("#FFFFFF","#E4E4FF","#FFFFE4","#FFE4E4")
                    ))
    }

    pal <- function(N, pname){
        if(pname %in% c("rainbow_hcl","heat_hcl","terrain_hcl","diverging_hcl")){
            do.call(pname, list(N))
        } else {
            brewer.pal(N, pname)
        }
    }

    catcolor <- function(cname, pname){
        if(cname %in% colnames(sdf)){
            html <- html %>%
                formatStyle(cname,
                            backgroundColor=styleEqual(
                                levels(sdf[[cname]]),
                                pal(nlevels(sdf[[cname]]), pname)
                                )
                            )
            
        }
        return(html)
    }

    html <- catcolor("dataset", "Pastel1")
    html <- catcolor("location", "rainbow_hcl")
    html <- catcolor("RCM", "Pastel1")
    html <- catcolor("GCM", "Pastel2")
    html <- catcolor("lat", "PuOr")
    html <- catcolor("lon", "Spectral")

# 
# if("dataset" %in% colnames(sdf)){
#     html <- html %>%
#         formatStyle("dataset",
#                     backgroundColor=styleEqual(
#                         levels(sdf$dataset),
#                         brewer.pal(nlevels(sdf$dataset), "Pastel1")
#                         ))
# }
# 
# if("location" %in% colnames(sdf)){
#  html <- html %>%
#      formatStyle("location",
#                  color=styleEqual(
#                      levels(sdf$location),
#                      rainbow_hcl(nlevels(sdf$location))
#                      ))
# }
# 
# if("RCM" %in% colnames(sdf)){
#  html <- html %>%
#      formatStyle("RCM",
#                  color=styleEqual(
#                      levels(sdf[["RCM"]),
#                      brewer.pal(nlevels(sdf[["RCM"]]), "Pastel1")
#                      ))
# }
# 


#--#
#--#
#--### Red to blue bg color style for correlations, etc.
#--#redblue <- styleInterval(seq(-0.95,0.95,0.1),
#--#                         c(hsv(1,   (10:1)/10, 1),
#--#                           "#FFFFFF",
#--#                           hsv(2/3, (1:10)/10, 1)
#--#                           ))
#--#
#--### Red to blue bg color style for percentages
#--#pctredblue <- styleInterval(seq(2.5,97.5,5),
#--#                         c(hsv(1,   (10:1)/10, 1),
#--#                           "#FFFFFF",
#--#                           hsv(2/3, (1:10)/10, 1)
#--#                           ))
#--#
#--### Rainbow bg color style for days of the year
#--#doyrainbow <- styleInterval(1:365, rainbow_hcl(366, c=50, l=100))
#--#
#--#
#--#if(analysis == "gev"){
#--#    ## background bars for the different return levels
#--#    for(m in c("rlevlo", "rlevmid", "rlevhi")){
#--#        html <- html %>% formatStyle(m, background=styleColorBar(sdf[[m]], "#CCCCCC"))
#--#    }
#--#}
#--#
#--#
#--#
#--#if(analysis == "pfit"){
#--#
#--#    ## range of MAD metrics for background bars
#--#    vrange <- lapply(sdf[,c("freqmad","intmad","totmad")], range, na.rm=TRUE)
#--#    vrange <- lapply(vrange, function(x){x[1]<-0;x})  ## set min for MAD to zero
#--#
#--#    ## red-blue background for correlations
#--#    for(m in c("freqcor", "intcor", "totcor")){
#--#        html <- html %>% formatStyle(m, backgroundColor=redblue)
#--#    }
#--#
#--#    ## background bar from 0 to max for MADs
#--#    for(m in c("freqmad","intmad","totmad")){
#--#        html <- html %>% formatStyle(m, background=styleColorBar(vrange[[m]], "#CCCCCC"))
#--#    }
#--#}
#--#
#--#
#--#
#--#if(analysis == "tmmd"){
#--#
#--#    vrange <- lapply(sdf[,c("dmad", "dmaxval","dminval")], range)
#--#    vrange$dmad[1] <- 0  ## set min for MAD to zero
#--#
#--#    ## dmad, dminval, dmaxval: colorbars from min to max
#--#    for(d in c("dmad", "dminval", "dmaxval")){
#--#        html <- html %>%
#--#            formatStyle(d, background=styleColorBar(vrange[[d]], "#CCCCCC"))
#--#    }
#--#    
#--#    ## dpctpos: red-blue background 0:100
#--#    html <- html %>% formatStyle("dpctpos", backgroundColor=pctredblue)
#--#
#--#    ## minday & maxday: colorwheel
#--#    for(d in c("dminday", "dmaxday")){
#--#        html <- html %>% formatStyle(d, backgroundColor=doyrainbow)
#--#    }
#--#}
#--#    


    ## Save to file

    htmltools::save_html(html, file=outfile, lib=outdir, background="white")
}
