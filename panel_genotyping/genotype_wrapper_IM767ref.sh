#!/bin/bash
#SBATCH --job-name=genotype_wrapper		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=5		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=3:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/home/%u/logs/%x.%j.out		# Standard output and error log, make sure this folder exists
#SBATCH --error=/home/%u/logs/%x.%j.err		# Standard output and error log, make sure this folder exists

#prep and submission script for genotype_each.sh

###SETUP

PROJECT=Tn5_sequencing
BATCH_NAME=panel_IM767ref

SCRIPTS_DIR=~/Git_repos/${PROJECT}/gutnas_panel
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

GENOME=~/Genomes/JGI_finals/Mguttatusvar_IM767_887_v2.0.fa
FASTQ_TABLE=${WORKING_DIR}/gutnas_panel_fastq_table_bigger.txt #two columns, directory and prefix

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p genotyping_temp
mkdir -p gvcfs


###prep genome
printf "\n...preparing genome\n" | tee >(cat >&2)
if [ ! -f ${GENOME}.fai ]; then
	ml SAMtools/1.16.1-GCC-11.3.0
	samtools faidx $GENOME
	module purge
fi
if [ ! -f ${GENOME}.amb ]; then
	ml BWA/0.7.17-GCCcore-11.3.0
	bwa index $GENOME
	module purge
fi
if [ ! -f ${GENOME%.fasta}.dict ]; then
	ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
	gatk CreateSequenceDictionary -R $GENOME
	module purge
fi

if [ -f ${WORKING_DIR}/gvcf_map.txt ]; then
	rm ${WORKING_DIR}/gvcf_map.txt
fi

REF_CODE=IM767_v2

###loop
while read -r line; do
	INPUT_PREFIX=$(echo -n $line | tr -s ' ' | cut -d' ' -f2)
	SEQBATCH="panel"
	printf "\n...submitting sample: $INPUT_PREFIX\n" | tee >(cat >&2)
	sbatch ${SCRIPTS_DIR}/genotype_each.sh $PROJECT $BATCH_NAME $GENOME $INPUT_PREFIX $SEQBATCH $REF_CODE
	printf "%s\t%s.%s.g.vcf.gz\n" $INPUT_PREFIX $INPUT_PREFIX $REF_CODE >> ${WORKING_DIR}/gvcf_map.txt
done < $FASTQ_TABLE

printf "\n...submitted all samples\n" | tee >(cat >&2)
