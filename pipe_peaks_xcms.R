args <- commandArgs(TRUE)

print(args)
options(digits=22)

filedir =  gsub(x = args[1],
                pattern = "\\\\", 
                replacement = "")
filename = args[2]
calib_path = args[3]
align=as.logical(args[4])
callmethod = args[5]
#scriptdir = args[6]
isCalled = if(args[6] == "no") FALSE else TRUE

# ---------------------

library(xcms)
library(MassSpecWavelet)
library(MALDIquant)
library(MSnbase)
library(data.table)
library(stringr)

if(calib_path != ""){
  refmz = list(pos = c(221.157514788, 235.173164852, 235.203357216, 255.226475392, 
                       291.235765108, 122.003584108, 100.02163948, 78.039694852, 381.367325968, 
                       403.36096562, 177.05405788, 126.10090575, 142.09582037, 138.043951166, 
                       116.062006538, 94.08006191, 216.07027431, 232.06518893, 171.098560022, 
                       187.093474642, 210.083244302, 172.106385054, 188.101299674, 171.074750568, 
                       80.130989718, 181.02750366, 159.045559032, 137.063614404, 157.091652216, 
                       135.109707588, 170.100358278, 148.11841365, 126.136469022, 133.068699784, 
                       195.043153724, 173.061209096, 136.050606916, 131.056146292, 197.041045278, 
                       151.079264468, 118.083158482, 175.05910065, 117.099142904, 179.073596844, 
                       89.115256052, 153.077156022, 202.129358792, 180.147414164, 224.11130342, 
                       107.070105306, 222.079410506, 157.102679984, 200.097465878, 135.120735356, 
                       178.11552125, 225.108729236, 161.088972144, 179.084624612, 153.158395832, 
                       193.150905144, 171.168960516, 229.123809352, 207.141864724, 235.130254234),
               neg = c(265.148441188, 112.00181969, 170.09183215, 186.08674677, 302.00399151, 
                       206.068509892, 222.063424512, 193.00768387, 250.96630624, 124.121916118, 
                       149.064711564, 207.023333934, 264.981956304, 148.058434818, 151.062603118, 
                       187.03928086, 133.095154684, 214.109539002, 178.13286126, 176.100968346, 
                       234.059590716, 292.018213086, 169.082860194, 212.077646088, 133.106182452
               ))
}else{
  refmz = list(pos=c(), neg=c())
  align=F
}

i=1
cwtJW <- function (ms, scales = 1, wavelet = "mexh"){
  if(i==1){
    print("Running Joanna-adjusted version...")
  }
  i <<- i + 1
  if (wavelet == "mexh") {
    psi_xval <- seq(-8, 8, length = 1024)
    psi <- (2/sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) *
      exp(-psi_xval^2/2)
  }
  else if (is.matrix(wavelet)) {
    if (nrow(wavelet) == 2) {
      psi_xval <- wavelet[1, ]
      psi <- wavelet[2, ]
    }
    else if (ncol(wavelet) == 2) {
      psi_xval <- wavelet[, 1]
      psi <- wavelet[, 2]
    }
    else {
      stop("Unsupported wavelet format!")
    }
  }
  else {
    stop("Unsupported wavelet!")
  }
  oldLen <- length(ms)
  ms <- MassSpecWavelet:::extendNBase(ms, nLevel = NULL, base = 2)
  len <- length(ms)
  nbscales <- length(scales)
  wCoefs <- NULL
  psi_xval <- psi_xval - psi_xval[1]
  dxval <- psi_xval[2]
  xmax <- psi_xval[length(psi_xval)]
  scales = scales[scales > 0]
  for (i in 1:length(scales)) {
    scale.i <- scales[i]
    f <- rep(0, len)
    j <- 1 + floor((0:(scale.i * xmax))/(scale.i * dxval))
    if (length(j) == 1)
      j <- c(1, 1)
    lenWave <- length(j)
    f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
    if (length(f) > len)
      stop(paste("scale", scale.i, "is too large!"))
    wCoefs.i <- 1/sqrt(scale.i) * convolve(ms, f)
    wCoefs.i <- c(wCoefs.i[(len - floor(lenWave/2) + 1):len],
                  wCoefs.i[1:(len - floor(lenWave/2))])
    wCoefs <- cbind(wCoefs, wCoefs.i)
  }
  if (length(scales) == 1)
    wCoefs <- matrix(wCoefs, ncol = 1)
  colnames(wCoefs) <- scales
  wCoefs <- wCoefs[1:oldLen, , drop = FALSE]
  return(wCoefs)
}

assignInNamespace("cwt", cwtJW, ns="MassSpecWavelet")

for(mode in c("pos", "neg")){
  
  print(paste("finding peaks in", mode, "mode"))
  
  mzlist = refmz[[mode]]

  mzdata_file = gsub(file.path(filedir, "sum", basename(filename)),
                     pattern = "\\.raw",
                     replacement=paste0("_",mode,".mzML"))
  
  if(!file.exists(mzdata_file)){ 
    stop("File doesn't exist!")
  }
  
  summed_file <- readMSData(files = mzdata_file,
                            mode = "onDisk")
  
  print("loaded file...")
   
  if(isCalled){
    print("already called")
  }else{
    if(callmethod == "msnbase"){
      msw <- MSWParam(scales=c(1,5,7))
      summed_peaks <- findChromPeaks(summed_file, 
                                     param = msw) 
      smoothed = smooth(summed_file, method = "SavitzkyGolay", halfWindowSize = 2)
      peaked = pickPeaks(smoothed, refineMz = "descendPeak")
      peaks = data.table(mz = peaked[[1]]@mz,
                         int = peaked[[1]]@intensity)
      unlockBinding(sym = "chromPeakData", env = summed_peaks@msFeatureData@.xData)
      unlockBinding(sym = "chromPeaks", env = summed_peaks@msFeatureData@.xData)
      reformatted = data.frame(row.names = paste0("CP", str_pad(1:nrow(peaks), 3, pad = "0")),
                               mz = peaks$mz,
                               mzmin = peaks$mz,
                               mzmax = peaks$mz,
                               rt = c(-1),
                               rtmin = c(-1),
                               rtmax = c(-1),
                               into = peaks$int,
                               maxo = peaks$int,
                               sn = c(1337),
                               intf = c(NA),
                               maxf = c(NA),
                               sample = c(1))
      summed_peaks@msFeatureData@.xData$chromPeakData <- DataFrame(ms_level=as.integer(rep(1, nrow(peaks))),
                                                                   is_filled = c(FALSE),
                                                                   row.names = rownames(reformatted))
      summed_peaks@msFeatureData@.xData$chromPeaks <- as.matrix(reformatted)
    }else if(callmethod == "msw"){
      msw <- MSWParam(snthresh = 0, 
                      verboseColumns = FALSE, 
                      scales = seq(1, 7, .5),
                      nearbyPeak = T, 
                      peakScaleRange = 0.5,
                      minNoiseLevel = 0, 
                      ridgeLength = 24, 
                      peakThr = 100,
                      ampTh = 1e-40,
                      tuneIn = T)
      
      summed_peaks <- findChromPeaks(summed_file, 
                                     param = msw)  
    }
  
  if(align){
    try({
      calib = CalibrantMassParam(mz = mzlist, mzppm = 3)
      summed_peaks = calibrate(summed_peaks, param=calib)
    }) 
  }
  
  fn = file.path(filedir, 
                 "peaks", 
                 gsub(basename(mzdata_file), 
                      pattern="\\.mzML",
                      replacement=".RData"))
  
  print("=====================================")
  
  print("writing...")
  
  if(file.exists(fn)) file.remove(fn)
  save(summed_peaks, 
       file = fn)
}

