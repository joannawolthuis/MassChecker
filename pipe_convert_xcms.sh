#!/bin/bash -l
# Default in slurm
# Request 24 hours run time
#SBATCH -t 0:30:0
#SBATCH --mem=20G

index="$5/rawfiles.txt"
filename=`sed -n "${SLURM_ARRAY_TASK_ID}p" < "$index"`
basename=`basename "$filename"`
echo $basename

if [ $5 == "yes" ]
then
  peakCall=""
else
  peakCall="-p"
fi

echo "writing output to $5"

guixr load-profile $HOME/.guix-profile/ -- <<EOF
  echo "Running inside a protected ecelnvironment";
  mono $2 -i="$filename" -o="$5/converted/" $peakCall
  exit
EOF

