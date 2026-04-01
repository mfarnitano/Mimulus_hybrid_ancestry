#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p batch --mem=8000 
#SBATCH --time=05:00:00

WORKING_DIR=/scratch/mcf96392/ancestry_IM767_v2
SCRIPTS_PATH=~/ancestry_pipeline_scripts/		###folder with python scripts

cd ${WORKING_DIR}
mkdir ${WORKING_DIR}/ANCESTRY_OUTPUTS

cp ancestry-probs* ANCESTRY_OUTPUTS/
cp all.indivs.hmm.combined_allchrs* ANCESTRY_OUTPUTS/

python ${SCRIPTS_PATH}/summarize_ancestry_bysample.py ${WORKING_DIR}/ANCESTRY_OUTPUTS/ancestry-probs_allchrs.tsv_rec.txt ${WORKING_DIR}/ANCESTRY_OUTPUTS/ancestry_sample_summary.txt
python ${SCRIPTS_PATH}/summarize_ancestry_bysite.py ${WORKING_DIR}/ANCESTRY_OUTPUTS/ancestry-probs_allchrs.tsv_rec.txt ${WORKING_DIR}/ANCESTRY_OUTPUTS/ancestry_site_summary.txt
