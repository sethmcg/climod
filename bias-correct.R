library(ncdf4)
library(KernSmooth)
library(devtools)
load_all()

## width of moving window
mwinwidth = 30

ocf <- c("obs","cur","fut")

for(varname in c("prec","tmax","tmin")){
#varname <- "prec"

norm <- ifelse(varname == "prec", "power", "zscore")

print("#################")
print(paste(varname))


indir <- "tests/raw"
infiles <- paste0(indir, "/", varname, ".", ocf, ".nc")
names(infiles) <- ocf

#infile.o <- paste0("tests/raw/",varname,".obs.nc")
#infile.c <- paste0("tests/raw/",varname,".cur.nc")
#infile.f <- paste0("tests/raw/",varname,".fut.nc")
#
#fin.o <- nc_open(infile.o)
#fin.c <- nc_open(infile.c)
#fin.f <- nc_open(infile.f)
        

nc <- lapply(infiles, nc_ingest)


## get netcdf data in friendly form
#nc.o <- nc_ingest(fin.o)
#nc.c <- nc_ingest(fin.c)
#nc.f <- nc_ingest(fin.f)
#>#nc <- lapply(fin, nc_ingest)


time <- lapply(nc,"[[","time")

indata <- lapply(nc,"[[",varname)

outdata <- lapply(indata,"+",NA)

#
### set up arrays for storing results
#save.c <- array(dim=dim(nc.c[[varname]]))
#save.f <- array(dim=dim(nc.f[[varname]]))
#


cwin <- lapply(time, cslice, outer=mwinwidth, num=360)

### generate climatology moving window index arrays
### 360 1-ish-day moving windows w/ 30-day outer pool
#cwin.o <- cslice(nc.o$time, outer=30, num=360)
#cwin.c <- cslice(nc.c$time, outer=30, num=360)
#cwin.f <- cslice(nc.f$time, outer=30, num=360)
#

    
## iterate over space
#xname <- names(dimnames(nc.o[[varname]]))[1]
#yname <- names(dimnames(nc.o[[varname]]))[2]

#nx <- length(nc.o[[xname]])
#ny <- length(nc.o[[yname]])

### Assert all 3 datasets have same spatial dimensions and are ordered
### lon, lat, time...

nxy <- dim(indata$obs)[1:2]
nx <- nxy[1]
ny <- nxy[2]


#x = 2
#y = 1


for(x in 1:nx){
    for(y in 1:ny){

        print(paste0("x:",x,", y:",y))

        data <- lapply(indata, function(a){a[x,y,]})
        
#        
#        ## get data from nc object 
#        data.o <- nc.o[[varname]][x,y,]
#        data.c <- nc.c[[varname]][x,y,]
#        data.f <- nc.f[[varname]][x,y,]



        ##  For stability, it's best to dedrizzle precip data here        
        if(varname == "prec"){

            ## floor data at zero
            data <- lapply(data, pmax, 0)

            ### MOVE TO DEDRIZZLE FUNCTION
            
            ## calculate wet/dry fraction
            pwet <- sum(data$obs > 0) / length(data$obs)

            ## find threshold equalizing model wet/dry with obs
            threshold <- sort(data$cur)[(1-pwet)*length(data$cur)]

            ## set values below threshold to zero
            data$cur[data$cur < threshold] <- 0
            data$fut[data$fut < threshold] <- 0

            ### END DEDRIZZLE FUNCTION

            ### MOVE TO UNZERO FUNCTION

            if(any(is.na(unlist(data)))){
                warning("Some values already NA")
            }
            data <- lapply(data, function(a){a[a == 0] <- NA; a})

            ### END UNZERO FUNCTION

        }
            
#        ## fix up precip
#        if(varname == "prec"){
#            ## floor all datasets at zero
#            data.o <- pmax(data.o, 0)
#            data.c <- pmax(data.c, 0)
#            data.f <- pmax(data.f, 0)
#            #>#data <- lapply(data,pmax,0)
#
#            ## calculate wet/dry fraction
#            nt <- length(data.o)
#            pwet <- sum(data.o > 0) / nt
#
#            ## find threshold equalizing model wet/dry with obs
#            threshold <- sort(data.c)[(1-pwet)*nt]
#
#            ## set values below threshold to zero
#            data.c[data.c < threshold] <- 0
#            data.f[data.f < threshold] <- 0
#
#            ## set zeros to NA for carrying through calculations
#            if(any(is.na(c(data.c, data.f, data.o)))){
#                warning("Some values already NA when setting zero -> NA")
#            }
#            data.o[data.o == 0] <- NA
#            data.c[data.c == 0] <- NA
#            data.f[data.f == 0] <- NA
#        }

        ## If all data for one input dataset is NA, result is NA.  No
        ## warnings or errors are needed here; all NA is expected over
        ## oceans, outside domain, etc.
        if(any(sapply(data, function(a){all(is.na(a))}))){
            next
        }
        
#        if(all(is.na(data.o)) || all(is.na(data.c)) || all(is.na(data.f))){
#            next
#        }

        
#        ## window data using outer window
        wind <- mapply(slice, data, cwin, MoreArgs=list(outer=TRUE), SIMPLIFY=FALSE)
        
#        wind.o <- slice(data.o, cwin.o, outer=TRUE)
#        wind.c <- slice(data.c, cwin.c, outer=TRUE)
#        wind.f <- slice(data.f, cwin.f, outer=TRUE)
        
        
        ## normalize each window separately
        nwd <- rapply(wind, normalize, norm=norm, how="replace")
        
#        nwd.o   <- lapply(wind.o, normalize, norm=norm)
#        nwd.c   <- lapply(wind.c, normalize, norm=norm)
#        nwd.f   <- lapply(wind.f, normalize, norm=norm)


        ## construct distribution mappings
        dmaps <- mapply(distmap, nwd$cur, nwd$obs, SIMPLIFY=FALSE)
        
#        dmaps <- mapply(distmap, nwd.c, nwd.o, SIMPLIFY=FALSE)

        
        ## get normalized inner window for bias correction
        fixme <- mapply(subslice, nwd, cwin, SIMPLIFY=FALSE)

        ## drop obs from fixme

        fixme <- fixme[-(names(fixme)=="obs")]
        
#        fixme.c <- subslice(nwd$cur, cwin$cur)
#        fixme.f <- subslice(nwd$fut, cwin$fut)
        
        ## bias-correct innner window
        fixed <- lapply(fixme, function(x){mapply(predict, dmaps, x)})

        
#        fixed.c <- mapply(predict, dmaps, fixme.c)
#        fixed.f <- mapply(predict, dmaps, fixme.f)


        ## copy atts from nwd to fixed

        fixed$cur <- mapply(copyatts, nwd$cur, fixed$cur, SIMPLIFY=FALSE)
        fixed$fut <- mapply(copyatts, nwd$fut, fixed$fut, SIMPLIFY=FALSE)
        
#        fixed.c <- mapply(copyatts, nwd$cur, fixed.c, SIMPLIFY=FALSE)
        ## whoops! nwd$cur here is an error; should be $fut
#        fixed.f <- mapply(copyatts, nwd$cur, fixed.f, SIMPLIFY=FALSE)

        

        if(varname == "prec"){

            ## MOVE TO REZERO FUNCTION
            
            ## fixed values may be negative; set to zero.
            ## (Don't warn; this is normal)
            fixed <- rapply(fixed, function(x){pmax(x,0)}, how="replace")
            
#            neg.c <- sum(unlist(fixed.c) < 0, na.rm=TRUE)
#            if(neg.c > 0){
#                warning(paste(neg.c,"negative values in current after bias correction"))
#                fixed.c <- lapply(fixed.c, pmax, 0)
#            }
#            neg.f <- sum(unlist(fixed.f) < 0, na.rm=TRUE)
#            if(neg.f > 0){
#                warning(paste(neg.f,"negative values in future after bias correction"))
#                fixed.f <- lapply(fixed.f, pmax, 0)
#            }

          
            ## convert NA back to zero
            fixed <- rapply(fixed, function(x){x[is.na(x)]<-0;return(x)}, how="replace")
            
#            fixed.c <- lapply(fixed.c, function(x){x[is.na(x)]<-0;return(x)})
#            fixed.f <- lapply(fixed.f, function(x){x[is.na(x)]<-0;return(x)})
            
            ## END REZERO FUNCTION
          
        }

        ## denormalize bias-corrected data
        bc <- list()
        bc$obs <- subslice(wind$obs, cwin$obs)
        bc$cur <- mapply(denormalize, fixed$cur, match=nwd$obs)
        bc$fut <- mapply(denormalize, fixed$fut, match=nwd$obs, adjust=nwd$cur)
        
#        bc.c <- mapply(denormalize, fixed.c, nwd$obs)
#        bc.f <- mapply(denormalize, fixed.f, nwd$obs, nwd$cur)

                
        ## collate BC inner windows back into timeseries
#        result.c <- unslice(bc.c, cwin$cur)
#        result.f <- unslice(bc.f, cwin$fut)

        result <- mapply(unslice, bc, cwin)
        
        
        ## save results in array
        for(i in names(result)){
            outdata[[i]][x,y,] <- result[[i]]
        }
#        outdata$cur[x,y,] <- result.c
#        outdata$fut[x,y,] <- result.f
    }
}
    
###################



outdir <- "tests"
outfiles <- paste0(outdir, "/", varname, ".", ocf, ".bc.nc")
names(outfiles) <- ocf

#outfile.c <- paste0("tests/",varname,".cur.bc.nc")
#outfile.f <- paste0("tests/",varname,".fut.bc.nc")


for(i in c("cur","fut")){

### BEGIN NC_REPLACE FUNCTION
    

    file.copy(infiles[i], outfiles[i], overwrite=TRUE, copy.mode=FALSE)
#file.copy(infile.c, outfile.c, overwrite=TRUE, copy.mode=FALSE)
#file.copy(infile.f, outfile.f, overwrite=TRUE, copy.mode=FALSE)
    
    fout <- nc_open(outfiles[i], write=TRUE)
#fout.c <- nc_open(outfile.c, write=TRUE)
#fout.f <- nc_open(outfile.f, write=TRUE)
    
    ncvar_put(fout, varname, outdata[[i]])
#ncvar_put(fout.c, varname, save.c)
#ncvar_put(fout.f, varname, save.f)

    ncatt_put(fout, varname, "bias_correction", "KDDM")
#ncatt_put(fout.c, varname, "bias_correction", "KDDM")
#ncatt_put(fout.f, varname, "bias_correction", "KDDM")


## BEGIN APPEND HISTORY FUNCTION

    historyatt <- ncatt_get(fout,0,"history")$value
#history.c = ncatt_get(fout.c,0,"history")$value
#history.f = ncatt_get(fout.f,0,"history")$value

    historyatt = paste0(date(),
        ": bias-corrected using R package 'climod'\n",
        "    see bias_correction attribute for details\n",
        historyatt)   
#history.c = paste0(date(),
#    ": bias-corrected using R package 'climod'\n",
#    "    see bias_correction attribute for details\n",
#    history.c)
#history.f = paste0(date(),
#    ": bias-corrected using R package 'climod'\n",
#    "    see bias_correction attribute for details\n",
#    history.f)

    ncatt_put(fout, 0, "history", historyatt)    
#ncatt_put(fout.c,0,"history",history.c)
#ncatt_put(fout.f,0,"history",history.f)

## END APPEND HISTORY FUNCTION

#bcnote = paste(
#    "Bias-corrected using KDDM (kernel density distribution mapping)\n",
#    "calculated over a 30-day moving window across years\n",
#    "observed data from",infile.o,"\n",
#    "current data from",infile.c,"\n",
#    "future data from",infile.f)
bcnote = paste(
    "Bias-corrected using kernel density distribution mapping\n",
    "calculated over a", mwinwidth, "day moving window across years\n",
    "observed data from",infiles["obs"],"\n",
    "current data from",infiles["cur"],"\n",
    "future data from",infiles["fut"])


    ncatt_put(fout, 0, "bias_correction", bcnote)
#ncatt_put(fout.c,0,"bias_correction",bcnote)
#ncatt_put(fout.f,0,"bias_correction",bcnote)
    
    nc_close(fout)
#nc_close(fout.c)
#nc_close(fout.f)

    ## END NC REPLACE FUNCTION
}


}
