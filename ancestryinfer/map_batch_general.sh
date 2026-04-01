#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=batch
#SBATCH --mem=64000
#SBATCH --time=7-00:00:00
#SBATCH --array=1-4%4

ml ancestryinfer/20230111-foss-2022a-Python-2.7.18

WORKING_DIR=/scratch/mcf96392/ancestry_IM767_v2
#READLIST=~/ancestry_pipeline_scripts/CAC_LM_moms_2019-2022_readlist.txt
READLIST=~/ancestry_pipeline_scripts/North_south_moms_readlist.txt
GENOME1=/home/mcf96392/Genomes/JGI_finals_fixed/Mguttatusvar_IM767_887_v2.0.fixed.fa
GENOME2=/home/mcf96392/Genomes/JGI_finals_fixed/Mnasutusvar_SF_822_v2.0.fixed.fa

READLIST_LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${READLIST}) #pulls one line from input list for each array task

cd $WORKING_DIR

READLIST_X="readlist_""$SLURM_ARRAY_TASK_ID"".txt"
printf "%s\n" "$READLIST_LINE" > $READLIST_X

run_map_v3.pl $READLIST_X $GENOME1 $GENOME2 PE _allchrs
