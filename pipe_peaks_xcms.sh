#!/bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -l h_rt=00:45:00
#$ -l h_vmem=20G

index="$2/rawfiles.txt"
filename=`sed -n "${SLURM_ARRAY_TASK_ID}p" < "$index"`
basename=`basename "$filename"`
echo $basename

guixr load-profile $HOME/.guix-profile/ -- <<EOF
  echo "Running inside a protected environment";
  Rscript --no-save $1/pipe_peaks_xcms.R "$2" "$basename" "$3" "$4" "$5" "$6" 
  exit
EOF

