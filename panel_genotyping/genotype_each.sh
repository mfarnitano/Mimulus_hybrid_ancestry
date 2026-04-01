#!/bin/bash
#SBATCH --job-name=genotype_each		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=40gb			                                # Total memory for job
#SBATCH --time=168:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/home/%u/logs/%x.%j.out		# Standard output and error log, make sure this folder exists
#SBATCH --error=/home/%u/logs/%x.%j.err		# Standard output and error log, make sure this folder exists

###Creates gvcf for single sample from filtered bam, requires input variables
###See submission script genotype_wrapper.sh

###SETUP
PROJECT=$1
BATCH_NAME=$2
GENOME=$3
INPUT_PREFIX=$4
SEQ_BATCH=$5
REF_CODE=$6

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=10
MEM=40

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
cd $WORKING_DIR
mkdir -p ${WORKING_DIR}/qualimap
mkdir -p ${WORKING_DIR}/genotyping_temp
mkdir -p ${WORKING_DIR}/gvcfs

###MODULES #updated for cluster change 09/23
# ml SAMtools/1.16.1-GCC-11.3.0
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
# ml Qualimap/2.2.1-foss-2021b-R-4.1.2
module list

printf "\nstarting analysis for sample ${INPUT_PREFIX}\n" | tee >(cat >&2)

###individual genotyping with GATK
printf "\n...adding read groups\n" | tee >(cat >&2)
gatk AddOrReplaceReadGroups \
	-I ${WORKING_DIR}/bams/${INPUT_PREFIX}.${REF_CODE}.fds.bam \
	-O $WORKING_DIR/genotyping_temp/${INPUT_PREFIX}.${REF_CODE}.rg.bam \
	-RGID $INPUT_PREFIX \
	-LB ${SEQ_BATCH} \
	-PL illumina \
	-PU $INPUT_PREFIX \
	-SM $INPUT_PREFIX

printf "\n...indexing bam file\n" | tee >(cat >&2)
module purge
ml SAMtools/1.16.1-GCC-11.3.0

samtools index $WORKING_DIR/genotyping_temp/${INPUT_PREFIX}.${REF_CODE}.rg.bam

module purge
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17

printf "\n...genotyping with gatk HaplotypeCaller\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" HaplotypeCaller  \
	-R $GENOME \
	-I $WORKING_DIR/genotyping_temp/${INPUT_PREFIX}.${REF_CODE}.rg.bam \
	-O $WORKING_DIR/gvcfs/${INPUT_PREFIX}.${REF_CODE}.g.vcf.gz \
	-ERC GVCF

printf "...Finished sample ${INPUT_PREFIX}\n" | tee >(cat >&2)
