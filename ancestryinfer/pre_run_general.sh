#!/bin/sh 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1 
#SBATCH -p batch --mem=96000 
#SBATCH --time=05:00:00

AIMS_LIST=~/ancestry_pipeline_scripts/IM767_AIMs.txt	#path to tab-separated file with chr, pos, ref, alt
READLIST=~/ancestry_pipeline_scripts/CAC_LM_moms_2019-2022_readlist.txt	#path to tab-separated file with read1_filename, read2_filename
WORKING_DIR=/scratch/mcf96392/ancestry_IM767_v2

MOD_SAMTOOLS_SCRIPT=~/ancestry_pipeline_scripts/run_samtools_to_hmm_v10_mod.pl

mkdir -p $WORKING_DIR
cd $WORKING_DIR

cp $MOD_SAMTOOLS_SCRIPT $WORKING_DIR
###Note: things to check before starting
###Use a genome that doesnt have underscores in chromosome names! Use `sed 's/Chr_/Chr/g' old > new` and re-index if necessary
###Check that AIMs list also doesnt have underscores in chromosome names
###fastq files must be in the working directory, and must be named PREFIX_read_1.fastq.gz, PREFIX_read_2.fastq.gz
 
###prep some files
sed 's/\t/_/g' $AIMS_LIST > AIMs_list_ready
awk 'BEGIN{OFS="\t"}{print $1 "_" $2, $2, $3, $4}' $AIMS_LIST > AIMs_list_ready.mod
awk 'BEGIN{OFS="\t"}{print $1, $2, $2, $3, $4}' $AIMS_LIST > AIMs_list_ready.mod.bed
printf "AIMs_list_ready\n" > current_aims_file

cut -f1 $READLIST | sed 's/$/_allchrs.sam.hmm.combinedparental.format/' > HMM.parental.files.list_allchrs
cut -f1 $READLIST | sed 's/$/_allchrs.sam.hmm.combined.pass.formatted/' > HMM.hybrid.files.list_allchrs

