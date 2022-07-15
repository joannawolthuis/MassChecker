# ================= GROUP PEAKS =================

proj_short = "example"
grouper = "maldiq" #alternatively, mzclust

library(xcms)
library(MassSpecWavelet)
library(MALDIquant)
library(MSnbase)
library(stringr)
library(Biobase)
library(data.table)

### PRIOR TO RUNNING THIS SCRIPT ###
# In the RAW file folder, place a csv named "samples.csv" indicating both which replicates belong to which sample and the sample injection order. 
# Can be with or without replicates (i.e. comma for each row: A, A, A, B, B, B, C, C, C OR A, B, C)

# indicate which folder has all the result folders used in grouping
basefolder = ""
foldergroups = list(all = c("folder_res_1",  # examples, add your own
                            "folder_res_2", 
                            "folder_res_3"))

folders = file.path(basefolder, foldergroups[[proj_short]])

maxrows = 2147483647

repl = 3 # technical replicate amount
ppm = 2 # ppm for peak grouping
regroup_only = F # if repeating a dataset done earlier, but only change is post peak collection

print(paste0("Grouping using ", grouper))

outdir = "" #output directory

centroided = F
projname = paste0(proj_short, "_", grouper, "_", ppm, "ppm")
projname_noppm = proj_short
process_replicates_mode = "both"

mapper = data.table::rbindlist(lapply(folders, function(folder){
  print(folder)
  res=data.table::data.table()
  try({
    files=list.files(file.path(folder, "converted"), pattern=".mzML")
    files_peaks=list.files(file.path(folder, "peaks"))
    print(paste("Files without peaks called:", length(files) - (length(files_peaks)/2)))
    prefix=lcPrefixC(files, ignore.case=FALSE)
    prefix=gsub("0$", "", prefix)
    csvloc=dirname(readLines(file.path(folder,
                                       "rawfiles.txt"),
                             n = 1))
    csv = file.path(csvloc, 
                    "samples.csv")
    if(length(csv)==0) stop("we'll need a csv with the sample names (1 per row) for this, please save it as .csv in the results folder!")
    if(length(csv)>1) csv = csv[1]
    info = data.table::fread(file = csv,header=F, blank.lines.skip=T, fill=TRUE)
    info = info[, sapply(1:ncol(info), function(i) !all(is.na(info[[i]]))), with=F]
    info = info[info[[1]] != ""]
    print(head(info))
    info = info[complete.cases(info),]
    if(FALSE){#ncol(info) == 2){
      res = data.table::data.table(filename=info[[1]],
                                   sample=info[[2]],
                                   sampgroup=c(1),
                                   injection=1:(length(filenames))) 
    }else{
      if(info[1,1]=="V1") info = data.table::fread(file = csv,sep=";",header=T)
      if(length(unique(unlist(info[,1][1:repl,]))) != 1){
        info = info[rep(1:nrow(info), each=repl)]
      }
      filenames = info[,1]
      filenames=gsub(unlist(filenames), pattern=" ", replacement="")
      filenames=gsub(unlist(filenames), pattern = "QC-", replacement = "QC")
      
      topad = nchar(stringr::str_match(string = files[1], ".*_(\\d+?)\\.mzML")[,2])
      numbers=str_pad(1:(length(filenames)), 
                      topad, 
                      pad = "0")
      
      res = data.table::data.table(filename=paste0(prefix,numbers),
                                   sample=filenames,
                                   sampgroup=if(ncol(info)==2) unlist(info[,2]) else c(1),
                                   injection=1:(length(filenames))) 
    }
  })
  res
}),use.names = T)

mapper <- mapper[complete.cases(mapper),]
print(mapper)
print(dim(mapper))

qcpos = grep("^QC", mapper$sample)

if(length(qcpos) > 0){
  i = 1
  split_qc = split(qcpos, ceiling(seq_along(qcpos)/repl))
  for(replgroup in split_qc){
    mapper$sample[replgroup] <- paste0("QC", i)
    i <<- i + 1
  }  
}

thr = 1

mapper$batch = gsub(mapper$filename, pattern = "_\\d+$", replacement = "")

# == DUPLICATES ==

uniqs = unique(mapper[,c("sample", "batch")])
dupliBatch <- c(which(duplicated(uniqs$sample)),
                which(duplicated(uniqs$sample,
                                 fromLast = T)))
duplis = uniqs[dupliBatch,]

removeMe <- unlist(lapply(split(duplis, duplis$sample), function(l){
  mergy = merge(l, mapper, by = c("sample", "batch"))
  splitty = split(mergy, mergy$batch)
  dates = names(splitty)
  conv = as.Date(gsub(names(splitty), pattern = "RES_|RES_BSP_|RES-", replacement = ""), format="%Y%m%d") # this checks for date based on part of our sample naming (RES_DATE_INJECTION)
  names(conv) = dates
  dumpBatch = names(which(conv == min(conv)))
  if(length(dumpBatch) > 0){
    splitty[[dumpBatch]]$filename
  }else{
    c()
  }
}))


if(length(removeMe) > 0){
  mapper = mapper[!(filename %in% removeMe)]
  print(removeMe)
}

switch(process_replicates_mode,
       noreps = NULL,
       repsonly = {
         mapper$sample = paste0(mapper$sample, "_", "REP", 1:repl)
       },
       both = {
         mapper_reps = mapper
         mapper_reps$sample = paste0(mapper_reps$sample, "_", "REP", 1:repl)
         mapper = rbind(mapper, 
                        mapper_reps)
       })
batch_mapper = mapper[,c(1,2,4,5)]

split.by.batch = split(batch_mapper, batch_mapper$batch)

batch_mapper_final = rbindlist(lapply(split.by.batch, function(bmap){
  switch(process_replicates_mode,
         noreps = {
           bmap_subs = unique(bmap[,c("sample", "batch")])
           bmap_subs$injection = 1:nrow(bmap_subs)
         },
         sepreps = {bmap_subs = bmap},
         both = {
           bmap_avged = unique(bmap[!grepl("REP\\d", sample),c("sample", "batch")])
           bmap_avged$injection = 1:nrow(bmap_avged)
           bmap_sep = unique(bmap[grepl("REP\\d", sample),c("sample", "batch")])
           bmap_sep$injection = 1:nrow(bmap_sep)
           bmap_subs = rbind(bmap_avged,
                             bmap_sep)
         })
  # - - - -
  bmap_subs
}))

extra = paste0("replhandling_", process_replicates_mode)
fn = file.path(outdir, paste0(projname, "_batches", extra, ".csv"))
fwrite(batch_mapper_final, file = fn)
print(paste("Saved batch info to", fn))

# ===========
filter_by_replicates=F

lapply(c("pos","neg"), function(mode){
  try({
    print(paste("Grouping peaks for", mode, "mode..."))
    peakfiles = list.files(file.path(folders, if(centroided) "sum" else "peaks"), 
                           pattern = mode, 
                           full.names = T)
    
    pkfiles_no_ext = gsub(x=basename(peakfiles), 
                          pattern=paste0("_",mode,"\\..*$"), 
                          replacement="",perl = T)
    
    remove = which(is.na(match(pkfiles_no_ext, mapper$filename)))
    
    if(length(remove) > 0){
      peakfiles = peakfiles[-remove]
    }
    
    if(!regroup_only){
      pklists <- pbapply::pblapply(1:length(peakfiles), function(i){
        
        load(peakfiles[i])
        pks = as.data.frame(chromPeaks(summed_peaks))[,c("sample", "mz", "into")]
        
        pks <- pks[pks$into > thr,]
        pks$mzmin <- pks$mz
        pks$mzmax <- pks$mz
        pks$rt <- c(NA)
        pks$rtmin <- c(NA)
        pks$rtmax <- c(NA)
        pks$sample <- gsub(basename(peakfiles[i]), 
                           pattern = "_(pos|neg)\\.mzML", 
                           replacement="")
        pks
      })
      allPeaks <- data.table::rbindlist(pklists, fill=T)
      data.table::fwrite(allPeaks, file = file.path(outdir, 
                                                    paste0(projname_noppm,
                                                           "_",
                                                           mode,
                                                           "_ungrouped.csv")))
    }else{
      allPeaks <- fread(file.path(outdir, paste0(projname_noppm,"_",mode,"_ungrouped.csv")))
    }
    
    print(dim(allPeaks))
    
    if(grouper == "mzclust"){
      
      data.table::setkey(allPeaks, mz)
      peaks_only <- data.table(mz = unique(allPeaks$mz))
      
      print(dim(peaks_only))
      
      minsamp = ceiling(0.01 * length(unique(allPeaks$sample)))
      
      # skip peaks that are only in 1% of samples to cut down a bit
      groups <- xcms:::mzClustGeneric(peaks_only,
                                      mzppm = ppm,
                                      minsamp = minsamp
      )
      
      grouptbl_rows <- pbapply::pblapply(1:length(groups$idx),
                                         function(i){
                                           res = data.table(mz = as.character(peaks_only[groups$idx[[i]]]$mz),
                                                            mzmed = as.character(groups$mat[i,"mzmed"]))
                                           res
                                         })  
      
      grouptbl <- rbindlist(grouptbl_rows)
      grouptbl$mz <- as.character(grouptbl$mz)
      setkey(grouptbl, "mz")
      
      allPeaks_dt <- as.data.table(allPeaks[,c("sample","mz","into")])
      
      allPeaks_dt$mz <- as.character(allPeaks_dt$mz)
      setkey(allPeaks_dt, "mz")
      
      merged = as.data.table(unique(merge(allPeaks_dt, 
                                          grouptbl,
                                          by = "mz")))
      
      merged[ , into := sum(as.numeric(into)), 
              by = list(sample, mzmed)]
      
      merged <- unique(merged[,c("sample",
                                 "mzmed",
                                 "into")])
      
    }else if(grouper == "maldiq"){
      
      require(MALDIquant)
      data.table::setkey(allPeaks, sample)
      
      samples = unique(allPeaks$sample)
      
      masspeaks = pbapply::pblapply(samples, function(samp){
        l = allPeaks[sample == samp]
        MALDIquant::createMassPeaks(l$mz, l$into, metaData = list(sample = samp))
      })
      
      binned = MALDIquant::binPeaks(l = masspeaks,
                                    method = "relaxed",
                                    tolerance = ppm * 1e-6)
      
      binned_allPeaks <- pbapply::pblapply(1:length(binned), function(i){
        l = binned[[i]]
        data.table::data.table(mzmed = l@mass,
                               sample = samples[i],
                               into = l@intensity)
      })
      
      merged = unique(data.table::rbindlist(binned_allPeaks))
    }
    
    merged$label=c("REPLACEME")
    mapper_essential <- unique(mapper[,c("sample","filename")])
    colnames(mapper_essential) <- c("sample_final", "sample")
    
    merged$sample <- gsub(x=merged$sample, 
                          pattern=paste0("_",mode,".*$"), 
                          replacement="", perl=T)
    
    setkey(mapper_essential, "sample")
    setkey(merged, "sample")
    
    merged_plus_mapper =  merge(merged, mapper_essential, by = "sample", allow.cartesian=T)
    merged_plus_mapper = unique(merged_plus_mapper[complete.cases(merged_plus_mapper),])
    
    keycols = c("sample_final", "mzmed")
    data.table::setkeyv(merged_plus_mapper, keycols)
  
    majcount=1 # how many replicates need to have a given m/z value
    
    samps = unique(merged_plus_mapper$sample_final)
    switch(process_replicates_mode,
           sepreps = {
             merged = merged_plus_mapper[,-"sample"]
           },
           noreps = {
             keep.per.sample = pbapply::pblapply(samps, function(samp){
               merged_ss = merged_plus_mapper[sample_final == samp,]
               mz_prev = table(merged_ss$mzmed)
               keep = names(mz_prev)[mz_prev >= majcount]
               filtered = merged_ss[mzmed %in% keep]
               filtered = filtered[,-"sample"]
               filtered[, .(into = mean(into)), by = c("mzmed", "sample_final", "label")]
             })  
             merged = data.table::rbindlist(keep.per.sample)
           },
           both = {
             keep.per.sample = pbapply::pblapply(samps, function(samp){
               merged_ss = merged_plus_mapper[sample_final == samp,]
               if(grepl("REP\\d", samp)){
                 merged = merged_ss[,-"sample"]
               }else{
                 mz_prev = table(merged_ss$mzmed)
                 keep = names(mz_prev)[mz_prev >= majcount]
                 filtered = merged_ss[mzmed %in% keep]
                 filtered = filtered[,-"sample"]
                 filtered[, .(into = mean(into)), by = c("mzmed", "sample_final", "label")]  
               }
             })  
             merged = data.table::rbindlist(keep.per.sample, use.names = T)
           })
    
    # ==========================
    
    merged[ , cast_cat := findInterval( mzmed , seq( 60 , 700 , 30 ) ) ]
    merged_list <- split(merged, by = 'cast_cat' )
    merged_list <- lapply(merged_list, function(x) x[, cast_cat := NULL])
    print(head(merged_list[[1]]))
    print("Merging peaklist parts...")
    merged_list <- pbapply::pblapply(merged_list, 
                                     function(z) data.table::dcast(z, 
                                                                   formula = sample_final + label ~ mzmed, 
                                                                   value.var = "into",
                                                                   fun.aggregate = sum))
    
    metTable <- Reduce(function( ... ) merge( ... , 
                                              by = c("sample_final", "label") , 
                                              all = TRUE),
                       merged_list)  
    
    peaktbl <- data.table::as.data.table(metTable)
    
    peaktbl <- peaktbl[!is.na(sample_final)]
    colnames(peaktbl)[colnames(peaktbl)=="sample_final"] <- "sample"
    
    print(dim(peaktbl))
    print(peaktbl[1:5,1:5])
    
    extra = paste0("_replhandling_", process_replicates_mode)
    
    fn = file.path(outdir, paste0(projname, "_", mode, extra, ".csv"))
    fwrite(peaktbl, file = fn)
    
    print(paste("Saved to", fn))
  })                                                                                                                                                                                                                                           
})

print("Done!")