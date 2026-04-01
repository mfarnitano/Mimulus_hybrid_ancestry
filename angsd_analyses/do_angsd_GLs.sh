#!/bin/bash
#SBATCH --job-name=angsd_GLs                     # Job name
#SBATCH --partition=highmem_p	                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=9		                            # Number of cores per task
#SBATCH --mem=480gb			                                # Total memory for job
#SBATCH --time=96:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/home/%u/logs/%x.%j.out		# Standard output and error log, make sure this folder exists
#SBATCH --error=/home/%u/logs/%x.%j.err		# Standard output and error log, make sure this folder exists

###SETUP
WORKING_DIR=$1
GROUPID=$2
BAMLIST=$3 #single column, full path to bams
SNPLIST=$4 #two columns, chr and pos
SNPLIST_TAG=$5
CHR=$6 #allChrs if blank

printf "Script called with inputs %s %s %s %s %s %s \n" $1 $2 $3 $4 $5 $6 | tee >(cat >&2)


if [ -z $CHR ]; then
	CHR=allChrs
fi
printf "CHR is %s\n" $CHR | tee >(cat >&2)

NTHREADS=8
MEM=480

###SETUP DIRS

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

###MODULES
ml angsd/0.940-GCC-12.3.0


#make hq SNPs list
if [ ! -f ${SNPLIST_TAG}.${CHR}.pos.txt ]; then
	if [ $CHR == "allChrs" ]; then
		cut -f1,2 $SNPLIST > ${SNPLIST_TAG}.${CHR}.pos.txt
		rm -f ${SNPLIST_TAG}.${CHR}.chrs.txt
		for i in Chr_{01..14}; do printf '%s\n' $i >> ${SNPLIST_TAG}.${CHR}.chrs.txt; done
	else
		grep $CHR $SNPLIST | cut -f1,2 > ${SNPLIST_TAG}.${CHR}.pos.txt
		echo $CHR > ${SNPLIST_TAG}.${CHR}.chrs.txt
	fi
	angsd sites index ${SNPLIST_TAG}.${CHR}.pos.txt
fi

printf "\n...starting genotyping for %s\n" ${SNPLIST_TAG}.${CHR} | tee >(cat >&2)
###run angsd to get genotype likelihoods
if [ ! -f ${WORKING_DIR}/genolike.${GROUPID}.${SNPLIST_TAG}.${CHR}.beagle.gz ]; then
	angsd -GL 2 -out genolike.${GROUPID}.${SNPLIST_TAG}.${CHR} -nThreads $NTHREADS \
		-rf ${SNPLIST_TAG}.${CHR}.chrs.txt -sites ${SNPLIST_TAG}.${CHR}.pos.txt \
		-doGlf 2 -doCounts 1 -doMajorMinor 1 -doMaf 1 -bam ${BAMLIST}

fi

printf "\n...finished genotyping for %s\n" ${SNPLIST_TAG}.$CHR | tee >(cat >&2)
