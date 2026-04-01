#!/bin/bash
#SBATCH --job-name=angsd_structure	                  # Job name
#SBATCH --partition=batch	                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=120gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/home/%u/logs/%x.%j.out		# Standard output and error log, make sure this folder exists
#SBATCH --error=/home/%u/logs/%x.%j.err		# Standard output and error log, make sure this folder exists

WORKING_DIR=$1
INPUT_BEAGLE=$2
OUT_PREFIX=$3

printf "Script called with parameters $1 $2 $3 \n"

ml angsd/0.940-GCC-12.3.0
ml PCAngsd/1.10

cd $WORKING_DIR
pcangsd -b $INPUT_BEAGLE -t 8 -o $OUT_PREFIX --inbreedSamples
