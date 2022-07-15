#!/bin/bash -l
#SBATCH -J GROUP_PEAKS
#SBATCH -o /home/cog/jwolthuis/group_peaks."%j".out
#SBATCH -e /home/cog/jwolthuis/group_peaks."%j".err
# Default in slurm
# Request 24 hours run time
#SBATCH -t 24:0:0
#SBATCH --mem=300G

guixr load-profile /home/cog/jwolthuis/.guix-profile/ -- <<EOF
  echo "Running inside a protected environment";
  Rscript --no-save /hpc/compgen/users/jwolthuis/MassChecker/pipeline_new/pipe_group_xcms.R
  exit
EOF

ß