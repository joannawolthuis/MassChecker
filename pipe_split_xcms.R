args <- commandArgs(TRUE)

print(args)
options(digits=22)

filedir = gsub(x = args[1], pattern = "\\\\", replacement = "")
filename = args[2]

library(data.table)

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

library(xcms)
library(MassSpecWavelet)
library(MALDIquant)
library(MSnbase)
library(pbapply)

mzdata_file = gsub(file.path(filedir, "converted", basename(filename)),
                   pattern = "\\.raw",
                   replacement=".mzML")

if(!file.exists(mzdata_file)){ 
  stop("file doesn't exist!")
}

raw_file <- readMSData(files = mzdata_file,
                       mode = "onDisk")

print("read file!")

pos_spec = raw_file[polarity(raw_file) == 1]
neg_spec = raw_file[polarity(raw_file) == 0]

outlier <- function (x,method="mean",addthres=FALSE){
  if (method=="mean") {
    avrg <- mean(x)
    stdev <-sd(x)
    dtf <<- data.frame(ID=seq.int(length(x)), obs=x, outlier=abs(x-avrg)>2*stdev)
    midp <<- avrg
    lower <<- avrg-2*stdev
    upper <<- avrg+2*stdev
    outliern <<- length(which(dtf=="TRUE"))
  } else {}
  if (method=="median") {
    med <- median(x)
    MAD <-median(abs(med-x))
    dtf <<- data.frame(ID=seq.int(length(x)), obs=x, outlier=abs(x-med)>2*(MAD/0.6745))
    midp <<- med
    lower <<- med-2*(MAD/0.6745)
    upper <<- med+2*(MAD/0.6745)
    outliern <<- length(which(dtf=="TRUE"))
  } else {}
  if (method=="boxplot") {
    Q1 <- quantile(x, 0.25)
    Q3 <- quantile(x, 0.75)
    IntQ <-Q3-Q1
    dtf <<- data.frame(ID=seq.int(length(x)), obs=x, outlier=x<Q1-1.5*IntQ | x>Q3+1.5*IntQ)
    midp <<- median(x)
    lower <<- Q1-1.5*IntQ
    upper <<- Q3+1.5*IntQ
    outliern <<- length(which(dtf=="TRUE"))
  } else {}
  return(dtf)
}

# --- remove outlier scans ---

noOuts = F

if(noOuts){
  
  tikkie = tic(pos_spec)
  outliers <- outlier(tikkie, method="median")
  keep_pos <- which(!outliers$outlier)
  
  tikkie = tic(neg_spec)
  outliers <- outlier(tikkie, method="median")
  keep_neg <- which(!outliers$outlier)
  
}else{
  
  keep_pos = 1:length(pos_spec)
  keep_neg = 1:length(neg_spec) 
  
}

min_scan_frac_with_signal = 0.7

if(length(keep_pos)>1){
  if((length(keep_pos)/length(pos_spec)) >= min_scan_frac_with_signal){
    pos_spec <- pos_spec[keep_pos]
  }else{
    stop("Too much pos spectra are outliers!")
  }
}
if(length(keep_neg)>1){
  if((length(keep_neg)/length(neg_spec)) >= min_scan_frac_with_signal){
    neg_spec <- neg_spec[keep_neg]
  }else{
    stop("Too much neg spectra are outliers!")
  }}

print("=====================================")
print("writing...")

require(R.utils)

# write positive mode
tries=3
currtry=0
res = NULL
while(is.null(res) & currtry < tries){
  currtry <<- currtry + 1
  try({
    withTimeout({
      fn_pos = gsub(file.path(filedir, "split", basename(mzdata_file)),
                    pattern = "\\.mzML",
                    replacement = paste0("_pos.mzML"))
      writeMSData(pos_spec, file = fn_pos)
    }, timeout = 300, onTimeout = "error")   
    res = TRUE
  })
}  
if(currtry == 3){
  print("failed to write + mode")
}else{
  print("writing success")
}

currtry=0
res = NULL
while(is.null(res) & currtry < tries){
  currtry <<- currtry + 1
  try({
    withTimeout({
      fn_neg = gsub(file.path(filedir,"split", basename(mzdata_file)),
                    pattern = "\\.mzML",
                    replacement=paste0("_neg.mzML"))
      writeMSData(neg_spec, file = fn_neg)
    }, timeout = 300, onTimeout = "error")   
    res = TRUE
  })
}  
if(currtry == 3){
  print("failed to write - mode")
}else{
  print("writing success")
}

