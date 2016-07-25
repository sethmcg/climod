library(ncdf4)
library(KernSmooth)
library(devtools)
load_all()

## width of moving window
mwinwidth = 30

ocf <- c("obs","cur","fut")

for(varname in c("prec","tmax","tmin")){
    
    norm <- ifelse(varname == "prec", "power", "zscore")
    
    print("#################")
    print(paste(varname))
    
    
    indir <- "tests/raw"
    infiles <- paste0(indir, "/", varname, ".", ocf, ".nc")
    names(infiles) <- ocf
            
    
    nc <- lapply(infiles, nc_ingest)
    
    
    ## get netcdf data in friendly form
    
    time <- lapply(nc,"[[","time")
    
    indata <- lapply(nc,"[[",varname)
    
    
    ### set up storage data structures
    
    outdata <- lapply(indata,"+",NA)
    
    
    ### generate climatology moving window index arrays
    ### 360 1-ish-day moving windows w/ 30-day outer pool
    cwin <- lapply(time, cslice, outer=mwinwidth, num=360)
        
    
    ### Assert all 3 datasets have same spatial dimensions and are ordered
    ### lon, lat, time...
    
    nxy <- dim(indata$obs)[1:2]
    nx <- nxy[1]
    ny <- nxy[2]
    
    
    for(x in 1:nx){
        for(y in 1:ny){
            
            print(paste0("x:",x,", y:",y))
            
            ## extract data for this gridcell
            data <- lapply(indata, function(a){a[x,y,]})
            
            if(varname == "prec"){
                ## Remove excess drizzle and set zero values to NA
                ## For stability, it's best to dedrizzle all at once, before slicing
                data <- dedrizzle(data)
                data <- lapply(data, unzero)
            }
            
            ## If all data for one input dataset is NA, result is NA.
            ## No warnings or errors are needed here; all NA is
            ## expected over oceans, outside domain, etc.
            
            if(any(sapply(data, function(a){all(is.na(a))}))){
                for(i in names(result)){
                    outdata[[i]][x,y,] <- NA
                }                
                next
            }
                        
            ## window data using outer window
            wind <- mapply(slice, data, cwin, MoreArgs=list(outer=TRUE), SIMPLIFY=FALSE)
            
            ## invert list nesting
            datatobc <- renest(wind)
            
            ## bias-correct each window
            fixdata <- lapply(datatobc, biascorrect, norm)
            
            ## re-invert the list
            bc <- renest(fixdata)
            
            ## collate inner windows back into timeseries
            result <- mapply(unslice, bc, cwin, SIMPLIFY=FALSE)
            
            ## rezero precipitation
            if(varname == "prec"){
                result <- rapply(result, rezero, how="replace")
            }
            
            ## save results in array
            for(i in names(result)){
                outdata[[i]][x,y,] <- result[[i]]
            }
        }
    }
    
    ###################
    
    outdir <- "tests"
    outfiles <- paste0(outdir, "/", varname, ".", ocf, ".bc.nc")
    names(outfiles) <- ocf
    
    
    for(i in c("cur","fut")){
        
        file.copy(infiles[i], outfiles[i], overwrite=TRUE, copy.mode=FALSE)
        
        fout <- nc_open(outfiles[i], write=TRUE)
        
        ncvar_put(fout, varname, outdata[[i]])
        
        ncatt_put(fout, varname, "bias_correction", "KDDM")
        
        nc_history(fout, paste("bias-corrected using R package 'climod'\n",
                               "    see bias_correction attribute for details"))
            
        bcnote = paste(
            "Bias-corrected using kernel density distribution mapping\n",
            "calculated over a", mwinwidth, "day moving window across years\n",
            "observed data from",infiles["obs"],"\n",
            "current data from",infiles["cur"],"\n",
            "future data from",infiles["fut"])
            
        ncatt_put(fout, 0, "bias_correction", bcnote)
        
        nc_close(fout)
    } 
}
