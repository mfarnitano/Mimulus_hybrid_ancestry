#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p batch --mem=96000
#SBATCH --time=72:00:00
#SBATCH --array=1-4%4

###load modules
ml ancestryinfer/20230111-foss-2022a-Python-2.7.18
ml BEDTools/2.30.0-GCC-11.3.0

###variables
WORKING_DIR=/scratch/mcf96392/ancestry_IM767_v2
READLIST=~/ancestry_pipeline_scripts/North_south_moms_readlist.txt
GENOME1=/home/mcf96392/Genomes/JGI_finals_fixed/Mguttatusvar_IM767_887_v2.0.fixed.fa
GENOME2=/home/mcf96392/Genomes/JGI_finals_fixed/Mnasutusvar_SF_822_v2.0.fixed.fa

cd $WORKING_DIR

###variables that probably don't need modification
READ_LENGTH=150
SAVE_FILES=1
MAX_ALIGN=2000000
FOCAL_CHROMS=0
REC_RATE=0.0000000387
PROGRAM_PATH=/apps/eb/ancestryinfer/20230111-foss-2022a-Python-2.7.18/bin
MIN_QUALITY=30

###pull sample row from array input readlist: readlist should be R1file(tab)R2file(newline)
READLIST_LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${READLIST}) #pulls one line from input list for each array task
READLIST_X="readlist_""$SLURM_ARRAY_TASK_ID"".txt"
printf "%s\n" "$READLIST_LINE" > $READLIST_X

###execute ancestryinfer script
perl run_samtools_to_hmm_v10_mod.pl $READLIST_X $GENOME1 $GENOME2 $READ_LENGTH $SAVE_FILES $MAX_ALIGN $FOCAL_CHROMS $REC_RATE $PROGRAM_PATH $MIN_QUALITY

