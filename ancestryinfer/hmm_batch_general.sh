#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p batch --mem=120000
#SBATCH --time=72:00:00

###load modules
ml ancestryinfer/20230111-foss-2022a-Python-2.7.18

###variables
WORKING_DIR=/scratch/mcf96392/ancestry_IM767_v2
AIMS_COUNTS=~/ancestry_pipeline_scripts/IM767_AIMs_counts_realcMs.txt

cd $WORKING_DIR

###variables that probably won't change
GENOME_PROP=0.5
GEN_ADMIX=0	#set=0 to have hmm estimate time
FOCAL_CHROMS=0
READ_LENGTH=150
ERROR_RATE=0.05
SAVE_FILES=1
POSTERIOR_THRESH=0.9
PROGRAM_PATH=/apps/eb/ancestryinfer/20230111-foss-2022a-Python-2.7.18/bin

###run ancestryHMM
combine_all_individuals_hmm_v5.pl HMM.parental.files.list_allchrs HMM.hybrid.files.list_allchrs $GENOME_PROP $AIMS_COUNTS $GEN_ADMIX $FOCAL_CHROMS $READ_LENGTH $ERROR_RATE _allchrs

###post-processing scripts
convert_rchmm_to_ancestry_tsv_v3.pl current.samples.list current.samples.read.list $SAVE_FILES $FOCAL_CHROMS
perl ${PROGRAM_PATH}/transpose_tsv.pl ancestry-probs-par1_transposed_allchrs.tsv
perl ${PROGRAM_PATH}/transpose_tsv.pl ancestry-probs-par2_transposed_allchrs.tsv
parsetsv_to_genotypes_v2.pl ancestry-probs-par1_allchrs.tsv ancestry-probs-par2_allchrs.tsv $POSTERIOR_THRESH ancestry-probs_allchrs.tsv_rec.txt
Rscript ${PROGRAM_PATH}/identify_intervals_ancestryinfer.R ancestry-probs_allchrs.tsv_rec.txt ${PROGRAM_PATH}
