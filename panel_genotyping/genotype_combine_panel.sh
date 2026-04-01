#!/bin/bash
#SBATCH --job-name=genotype_combine                     # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=32		                            # Number of cores per task
#SBATCH --mem=120gb			                                # Total memory for job
#SBATCH --time=168:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/home/%u/logs/%x.%j.out		# Standard output and error log, make sure this folder exists
#SBATCH --error=/home/%u/logs/%x.%j.err		# Standard output and error log, make sure this folder exists
#SBATCH --array=1-14%14

###SETUP
CHR_LIST=~/Genomes/JGI_finals/Mguttatusvar_IM767_887_v2.0.chrlist.txt
CHR=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${CHR_LIST}) #pulls one line from input list for each array task

BATCH_NAME=panel_IM767v2
WORKING_DIR=/scratch/mcf96392/panel_IM767ref
GENOME=~/Genomes/JGI_finals/Mguttatusvar_IM767_887_v2.0.fa

NTHREADS=32
MEM=110


###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p VCFs

###MODULES #fixed for 09/23 cluster update
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
module list

###create chr_list
printf "\n...creating chromosome list\n" | tee >(cat >&2)


##merge gvcfs into database
printf "\n...merging gvcfs into database, selecting only major (chromosome) linkage groups\n" | tee >(cat >&2)
cd ${WORKING_DIR}/gvcfs
gatk --java-options "-Xmx${MEM}G" GenomicsDBImport \
	--sample-name-map ${WORKING_DIR}/gvcf_map.txt \
	--genomicsdb-workspace-path ${WORKING_DIR}/${BATCH_NAME}_${CHR}_vcf_db \
	-L $CHR
cd $WORKING_DIR

###joint genotyping
printf "\n...Genotyping VCFs...outputing all-sites VCF with min call confidence = 30" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" GenotypeGVCFs \
	-R $GENOME \
	-L $CHR \
	-V gendb://${WORKING_DIR}/${BATCH_NAME}_${CHR}_vcf_db \
	-stand-call-conf 30 \
	-all-sites \
	-O ${WORKING_DIR}/VCFs/${BATCH_NAME}.${CHR}.allsites.vcf.gz

printf "\n...Finished genotyping\n" | tee >(cat >&2)
