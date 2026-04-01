#!/bin/bash
#SBATCH --job-name=fastq_align_each		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=128gb			                                # Total memory for job
#SBATCH --time=120:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/home/%u/logs/%x.%j.out		# Standard output and error log, make sure this folder exists
#SBATCH --error=/home/%u/logs/%x.%j.err		# Standard output and error log, make sure this folder exists

###Creates filtered bam file for single sample from raw fastq, requires input variables
###See submission script fastq_align_wrapper.sh

###SETUP
PROJECT=$1
BATCH_NAME=$2
GENOME=$3
INPUT_DIR=$4
INPUT_PREFIX=$5
REF_CODE=$6

SCRIPTS_DIR=~/Git_repos/${PROJECT}/gutnas_panel
LOG_DIR=/home/mcf96392/logs/gutnas_panel
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=10
MEM=128

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p ${WORKING_DIR}/fastqc
mkdir -p ${WORKING_DIR}/fastq_align_temp
mkdir -p ${WORKING_DIR}/bams

###MODULES #updated for cluster change 09/23
# ml FastQC/0.11.9-Java-11
# ml Trimmomatic/0.39-Java-13
# ml BWA/0.7.17-GCCcore-11.3.0
# ml SAMtools/1.16.1-GCC-11.3.0
# ml picard/2.27.5-Java-15
# ml Qualimap/2.2.1-foss-2021b-R-4.1.2
# module list

printf "\nstarting analysis for sample ${INPUT_PREFIX}\n" | tee >(cat >&2)

###fastqc
printf "...running fastqc" | tee >(cat >&2)
ml FastQC/0.11.9-Java-11
fastqc ${INPUT_DIR}/${INPUT_PREFIX}_read_1.fastq.gz ${INPUT_DIR}/${INPUT_PREFIX}_read_2.fastq.gz
mv ${INPUT_DIR}/${INPUT_PREFIX}_*_fastqc* ${WORKING_DIR}/fastqc/
module purge

###find overrepresented sequences
printf "\n...finding overrepresented sequences\n" | tee >(cat >&2)
for i in ${WORKING_DIR}/fastqc/${INPUT_PREFIX}_*.zip; do
	fdir=${i%.zip}
  unzip -o -q $i -d ${WORKING_DIR}/fastqc
  grep -A 2 ">>Overrepresented sequences" ${fdir}/fastqc_data.txt |
    tail -n1 |
    sed 's/ /_/' >> ${WORKING_DIR}/fastqc/${INPUT_PREFIX}_tempfile.txt
done

paste ${WORKING_DIR}/fastqc/${INPUT_PREFIX}_tempfile.txt ${WORKING_DIR}/fastqc/${INPUT_PREFIX}_tempfile.txt |
	grep -v "pass\|fail\|warn\|No_Hit"
  cut -f4,5 |
  sort -u |
  sed -e 's/^/>/' |
  sed 's/\t/\n/' > ${WORKING_DIR}/fastqc/${INPUT_PREFIX}_overrepresented.fa

 ###Trim
# if test "$(wc -l < ${WORKING_DIR}/fastqc/${INPUT_PREFIX}_overrepresented.fa)" -gt 0; then
printf "\n...trimming with Trimmomatic\n" | tee >(cat >&2)
ml Trimmomatic/0.39-Java-13
java -jar ${EBROOTTRIMMOMATIC}/trimmomatic-0.39.jar \
	PE -threads $NTHREADS \
	${INPUT_DIR}/${INPUT_PREFIX}_read_1.fastq.gz ${INPUT_DIR}/${INPUT_PREFIX}_read_2.fastq.gz \
	-baseout ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.fq.gz \
	ILLUMINACLIP:${WORKING_DIR}/fastqc/${INPUT_PREFIX}_overrepresented.fa:1:30:15 \
	SLIDINGWINDOW:5:20 \
	MINLEN:30
module purge

###map to reference and sort
printf "\n...mapping to reference with bwa mem and sorting with samtools\n" | tee >(cat >&2)
ml BWA/0.7.17-GCCcore-11.3.0
ml SAMtools/1.16.1-GCC-11.3.0
bwa mem -t $NTHREADS $GENOME \
	${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}_1P.fq.gz ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}_2P.fq.gz |
	samtools sort --threads $NTHREADS -O bam -o ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.${REF_CODE}.s.bam -

###index bam file
printf "\n...indexing bam file\n" | tee >(cat >&2)
samtools index ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.${REF_CODE}.s.bam
module purge

###mark and remove duplicates
printf "\n...marking and removing PCR duplicates with picard\n" | tee >(cat >&2)
ml picard/2.27.5-Java-15
java -jar ${EBROOTPICARD}/picard.jar \
	MarkDuplicates \
	I=${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.${REF_CODE}.s.bam \
	O=${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.${REF_CODE}.ds.bam \
	M=${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.dedupe-metrics.txt \
	REMOVE_DUPLICATES=true
module purge

###filter bam file and index
printf "\n...filtering bamfile for MAPQ>=29, both reads mapped and properly paired, passes platform QC\n" | tee >(cat >&2)
ml SAMtools/1.16.1-GCC-11.3.0
samtools view --threads $NTHREADS -q 29 -f 2 -F 524 -b ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.${REF_CODE}.ds.bam > ${WORKING_DIR}/bams/${INPUT_PREFIX}.${REF_CODE}.fds.bam

printf "\n...indexing filtered bam file\n" | tee >(cat >&2)
samtools index ${WORKING_DIR}/bams/${INPUT_PREFIX}.${REF_CODE}.fds.bam
module purge

printf "\n...generating coverage summary stats with qualimap\n" | tee >(cat >&2)
ml Qualimap/2.2.1-foss-2021b-R-4.1.2
qualimap bamqc -bam ${WORKING_DIR}/bams/${INPUT_PREFIX}.${REF_CODE}.fds.bam -c -outdir ${WORKING_DIR}/qualimap/${INPUT_PREFIX}
module purge

printf "\n...completed sample ${INPUT_PREFIX}" | tee >(cat >&2)
