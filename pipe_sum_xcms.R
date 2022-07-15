args <- commandArgs(TRUE)

options(digits=22)

filedir = gsub(x = args[1], pattern = "\\\\", replacement = "")
filename = args[2]
calib_path = args[3]
align=as.logical(args[4])
binme=as.logical(args[5])

print(dput(args))

if(calib_path != ""){
  calib_tbl = data.table::fread(calib_path)
  # LIST: FOUND IN >80% OF PLASMA SAMPLE SCANS MATCHING INTERNAL STANDARDS
  # IF USING OTHER SAMPLE TYPES PLEASE RENEW
  refmz = list(pos = c(214.107384422, 129.065853968, 240.120629178, 218.13868455, 
                       88.0390422366667, 258.15606597, 257.147178284, 129.089663422, 
                       235.165233656, 299.109375116, 118.086255054, 303.142913608, 223.10771878, 
                       267.13152822, 136.073679, 254.136279242, 245.149583592, 309.129242072, 
                       123.078430022, 121.086589412, 137.072624777, 249.18088372, 268.151929306, 
                       259.165233656, 130.086255054, 150.094945764, 90.05549407, 287.170319022, 
                       285.178478412, 389.154540932, 373.180603692, 263.196533784, 351.198659064, 
                       132.101905118, 387.186509776, 371.212572536, 349.230627908, 368.22520817, 
                       159.117752144, 273.18088372, 137.094080086, 287.196533784, 144.101905118, 
                       189.128316828, 145.086589412, 153.084046722, 164.104979128, 363.155062896, 
                       407.23879252, 231.138881512, 245.154531576, 394.292779882, 132.102444262, 
                       195.172330406, 463.301392776, 259.17018164, 400.342135382, 207.172330406, 
                       134.118896095333, 209.18798047, 273.185831704, 221.18798047, 
                       120.003194118, 98.02124949, 76.039304862, 138.043951166, 116.062006538, 
                       94.08006191, 88.039304862, 72.08077575, 151.071333272, 168.097882378, 
                       178.008673422, 156.026728794, 102.054954926, 134.044784166, 139.047798596, 
                       100.07569037, 116.07060499, 162.05014431, 194.02221504, 148.06043423, 
                       140.068199682, 172.040270412, 115.086589412, 177.061043352, 86.096425814, 
                       150.058325784, 155.079098724, 133.097154096, 165.086983336, 159.086589398, 
                       83.047129894, 135.11280416, 157.088555644, 135.106611016, 130.049869546, 
                       192.024323486, 170.042378858, 101.059705948, 159.076418652, 175.071333272, 
                       115.075356012, 114.091340434, 176.065794374, 177.127346446, 154.083849746, 
                       175.095142726, 198.084912386, 219.082841436, 180.086648978, 176.102967758, 
                       192.097882378, 174.111127148, 213.074834048, 197.100896808, 175.11895218, 
                       150.112469802, 182.084540532, 165.123368844, 193.129516864, 192.145501286, 
                       225.108729236, 189.086983336, 206.076359058, 181.094748788, 213.066819518, 
                       95.047129894, 184.09441443, 191.08487489, 196.10564783, 162.112469802, 
                       195.092365564, 227.064436294, 164.128119866, 208.129182506, 179.139018908, 
                       207.145166928, 120.08077575, 136.07569037, 216.097882378, 267.08229188, 
                       173.128454224, 252.12506269, 161.060631192, 210.05014431, 226.04505893, 
                       188.068199682, 220.037051542, 204.063114302, 166.086255054, 182.081169674, 
                       199.10771878, 226.104979114, 204.123034486, 76.039844006, 203.139018908, 
                       108.07876438, 221.149583592, 111.078430022, 259.150252294, 149.059705948, 
                       165.054620568),
               neg = c(233.150680752, 309.17795814, 369.198019632, 307.187460184, 
                       395.189585464, 381.224239616, 92.065509006, 91.05923226, 132.030231262, 
                       167.059278934, 184.001763194, 168.027825954, 146.045881326, 130.087352214, 
                       166.064029956, 209.077909736, 173.104399276, 209.081077018, 206.067010694, 
                       160.097916898, 208.064902248, 164.07170215, 200.048379892, 216.043294512
               ))
}else{
  refmz = list(pos=c(), neg=c())
  align=F
}

library(xcms)
library(MassSpecWavelet)
library(MALDIquant)
library(MSnbase)
library(pbapply)
library(R.utils)
library(data.table)

# binning?
lowMZ = 60
highMZ = 700
resol = 140000

nsegment = 2 * (highMZ-lowMZ)
segment = seq(from=lowMZ, to=highMZ, length.out=nsegment+1)
breaks.fwhm=NULL
breaks.fwhm.avg=NULL

breaks.list <- pbapply::pblapply(1:nsegment, function(i){
  startsegm <- segment[i]
  endsegm <- segment[i + 1]
  resol.mz <- resol*(1/sqrt(2)^(log2(startsegm/200)))
  fwhmsegm <- startsegm/resol.mz
  break.fwhm <- seq(from = (startsegm + fwhmsegm),to = endsegm, by = 0.2*fwhmsegm)
  # --- avg ---
  range = seq(from = (startsegm + fwhmsegm),to=endsegm, by=0.2*fwhmsegm)
  deltaMZ = range[2] - range[1]
  # --- return ---
  list(st=break.fwhm)
})

breaks.fwhm <- unlist(lapply(breaks.list, function(x) x$st))

binSum <- function(x){ # x is a list of spectra
  # z needs to be aligned?
  if(align){
    mzlist = refmz[[mode]]
    require(MALDIquant)
    maldi_spec = lapply(x, function(spec){
      createMassSpectrum(mass = spec@mz, intensity = spec@intensity)
    })
    refpeaks = createMassPeaks(mass = mzlist, 
                               intensity = rep(9999999999, length(mzlist)))
    
    aligned = alignSpectra(maldi_spec, reference = refpeaks,
                           tolerance = 5e-6, allowNoMatches = T, emptyNoMatches = T)
    x <- lapply(1:length(aligned), function(i){
      maldi_spec = aligned[[i]]
      xcms_spec = x[[i]]
      xcms_spec@mz <- maldi_spec@mass
      xcms_spec@intensity <- maldi_spec@intensity
      # - - - - 
      xcms_spec
    })  
  }
  # - - - - - - - - - - - 
  if(binme){
    binbreaks = breaks.fwhm
  }else{
    peaks_per_scan = lapply(x, function(spec){
      data.table(mass = spec@mz, 
                 intensity = spec@intensity)
    })
    peaks_all_scans = rbindlist(peaks_per_scan)
    peaks_reordered = unique(peaks_all_scans[order(mass, decreasing=F)])
    binbreaks = peaks_reordered$mass
  }
  binned <- lapply(x, function(z) bin(z, breaks = binbreaks))
  intensities = lapply(binned, function(z) z@intensity)
  summedInt = Reduce(`+`, intensities)
  new_sp <- binned[[1]] 
  new_sp@intensity <- summedInt  
  new_sp@peaksCount <- length(new_sp@mz)
  clean(new_sp)
}

for(mode in c("pos", "neg")){
  
  print(paste("Summing", mode, "mode spectra!"))
  
  mzdata_file = gsub(file.path(filedir, 
                               "split", 
                               basename(filename)),
                     pattern = "\\.raw",
                     replacement=paste0("_", mode, ".mzML"))
  
  if(!file.exists(mzdata_file)){ 
    stop("File doesn't exist!")
  }
  
  my_file <- readMSData(files = mzdata_file,
                        mode = "onDisk")

  print("loaded file")
  
  summed_file <- combineSpectra(Spectra(spectra(my_file)), 
                                method = binSum)

  if(summed_file[[1]]@peaksCount == 0){
    if(align){
      print("no spectra matched reference peaks! returning non-aligned spectra")
      align <<- F
      summed_file <- combineSpectra(Spectra(spectra(my_file)), fun = binSum)
      if(summed_file[[1]]@peaksCount == 0){
        stop("no peaks here!")
      }
    }else{
      stop("no peaks here!")
    }
  }
  
  print("summed spectra...")
  
  summed_spec = as(summed_file, "MSnExp")

  fn = file.path(filedir, "sum", basename(mzdata_file))
  print(fn)
  
  # write negative mode
  tries = 15
  currtry = 0
  res = NULL
  
  while(is.null(res) & currtry < tries){
    currtry <<- currtry + 1
    try({
      withTimeout({
        writeMSData(summed_spec, file = fn)
      }, timeout = 60, onTimeout = "error")
      res = TRUE
    })
  }
  
  if(currtry == tries){
    print("failed to write")
  }else{
    print("writing success!")
  }
}
