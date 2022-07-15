filedir="result directory of folder with only internal standards run through this pipeline"

library(xcms)
library(MassSpecWavelet)
library(MALDIquant)
library(MSnbase)
library(data.table)

name = "reference_peaks_bsp_is"
ppm = 1.5 

tblpermode = lapply(c("pos","neg"), function(mode){
  # LOAD IN ALL RELEVANT PEAK FILES
  peakfolder = file.path(filedir, "results", "peaks")
  peakfiles = list.files(peakfolder, pattern=paste0("_", mode),full.names = T)
  
  peaklist = pbapply::pblapply(peakfiles, function(pfile){
    load(pfile)  
    peaks_formatted = chromPeaks(summed_peaks)
    peaks_table = as.data.frame(peaks_formatted)
    print(peaks_table)
    MALDIquant::createMassPeaks(mass = peaks_table$mz,
                                intensity = peaks_table$into)
  })
  
  print(peaklist)
  # 'REFERENCEPEAKS'
  refpks = MALDIquant::referencePeaks(peaklist,
                                      tolerance = ppm*1e-6)
  
  print(refpks)
  # SAVE AS CSV  
  data.frame(mz = refpks@mass,
             mode = c(mode))  
})

refpks_tbl = rbindlist(tblpermode, fill=T)

fwrite(refpks_tbl, file = file.path(filedir, paste0(name, ".csv")))

