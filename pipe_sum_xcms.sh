#!/bin/bash -l
# Default in slurm
# Request 24 hours run time
#SBATCH -t 0:30:0
#SBATCH --mem=20G

index="$1/rawfiles.txt"
echo $index
filename=`sed -n "${SLURM_ARRAY_TASK_ID}p" < "$index"`
echo $filename
basename=`basename "$filename"`
echo $basename

guixr load-profile $HOME/.guix-profile/ -- <<EOF
  echo "Running inside a protected environment";
  Rscript --no-save $1/pipe_sum_xcms.R "$2" "$basename" "$3" "$4" "$5"
  exit
EOF
