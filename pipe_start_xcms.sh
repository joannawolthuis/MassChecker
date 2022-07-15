#SBATCH -J START_PIPELINE
#SBATCH -o $HOME/start."%j".out
#SBATCH -e $HOME/start."%j".err
# Default in slurm
# Request 24 hours run time
#SBATCH -t 0:30:0
#SBATCH --mem=20G

# specify folder with RAW files
filedir= ""
# specify folder with the pipe_ scripts
scriptdir=""
# specify output directory
outdir=""
# specify folder that has the Thermo parser and internal standard csv
pkgdir=""

# ===========================================

thermotool="$pkgdir/ThermoRawFileParser.exe"
calibrants="$pkgdir/80perc_intstds_100iso_5ppm.csv"
useThermoCalling="no"
callmethod="msnbase" #alternative is 'wavelet'
autokill=no # kill if dependency is invalid

echo "writing output to $outdir"

# create directories and submit SLURM jobs
rm "$filedir/*.out"
echo "$filedir"
jobsuffix=`dbus-uuidgen`
#echo $jobsuffix

foldername=`basename "$filedir"`
resdir="$outdir/$foldername"
mkdir "$resdir"

filelist="$resdir"/rawfiles.txt
      
if test -f "$filelist"; then
  rm "$filelist"
fi
      
ls "$filedir"/*.raw > "$filelist"
filecount=`wc -l "$filelist"`
maxfiles=`echo $filecount | cut -f1 -d" "`
      
# ===========================================

rm -r "$resdir"/converted
mkdir -p "$resdir"/converted
convid=$(sbatch -J convert_$jobsuffix --output=convert-%j.out --parsable --array=1-$maxfiles "$scriptdir"/pipe_convert_xcms.sh "$filedir" "$thermotool" "$useThermoCalling" "$resdir")
echo $convid

sleep 1

rm -r "$resdir"/split
mkdir -p "$resdir"/split
splitid=$(sbatch --parsable --kill-on-invalid-dep=$autokill --output=split-%j.out --dependency=afterany:$convid -J split_$jobsuffix --array=1-$maxfiles "$scriptdir"/pipe_split_xcms.sh "$resdir" "$calibrants")
echo $splitid

sleep 1

rm -r "$resdir"/sum
mkdir -p "$resdir"/sum
sumid=$(sbatch --parsable --kill-on-invalid-dep=$autokill --output=sum-%j.out --dependency=afterany:$splitid -J sum_$jobsuffix --array=1-$maxfiles "$scriptdir"/pipe_sum_xcms.sh "$resdir" "$calibrants" "TRUE" "TRUE")
echo $sumid
      
sleep 1
      
rm -r "$resdir"/peaks
mkdir -p "$resdir"/peaks
peakid=$(sbatch --parsable --kill-on-invalid-dep=$autokill --output=peaks-%j.out --dependency=afterany:$sumid -J peaks_$jobsuffix --array=1-$maxfiles "$scriptdir"/pipe_peaks_xcms.sh "$resdir" "$calibrants" "FALSE" "$callmethod" "$useThermoCalling")
echo $peakid

# next run pipe_group.R to create final output document