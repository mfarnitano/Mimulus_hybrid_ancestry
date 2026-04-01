#!/bin/bash
#SBATCH --job-name=genotype_filter	                # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=32		                            # Number of cores per task
#SBATCH --mem=60gb			                                # Total memory for job
#SBATCH --time=120:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/home/%u/logs/%x.%j.out		# Standard output and error log, make sure this folder exists
#SBATCH --error=/home/%u/logs/%x.%j.err		# Standard output and error log, make sure this folder exists
#SBATCH --array=1-14%14

###SETUP
CHR_LIST=~/Genomes/JGI_finals/Mguttatusvar_IM767_887_v2.0.chrlist.txt
CHR=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${CHR_LIST}) #pulls one line from input list for each array task
printf "\nCHR is %s\n\n" $CHR | tee >(cat >&2)

BATCH_NAME=panel_IM767v2
WORKING_DIR=/scratch/mcf96392/panel_IM767ref
GENOME=~/Genomes/JGI_finals/Mguttatusvar_IM767_887_v2.0.fa

SCRIPTS_DIR=~/Git_repos/Tn5_sequencing/gutnas_panel

GENE_BED=~/Genomes/JGI_finals/Mguttatusvar_IM767_887_v2.1.gene.gff3.1-based.bed
SAMPLELIST=${WORKING_DIR}/VCFs/gutnas_panel_key33_samples.args
POPFILE=${WORKING_DIR}/VCFs/gutnas_panel_key33_popfile.txt
##made with the following:
##printf "#chrom\tchromStart\tchromEnd\n" > Mguttatusvar_IM767_887_v2.1.gene.gff3.1-based.bed
##awk -v OFS='\t' '$3=="gene" {print $1,$4,$5}' Mguttatusvar_IM767_887_v2.1.gene.gff3 >> Mguttatusvar_IM767_887_v2.1.gene.gff3.1-based.bed

NTHREADS=32
MEM=50		##slightly less than requested

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p VCFs

###MODULES
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
ml HTSlib/1.15.1-GCC-11.3.0
ml BCFtools/1.15.1-GCC-11.3.0

module list

#inputs
RAW_VCF=${WORKING_DIR}/VCFs/${BATCH_NAME}.${CHR}.allsites.vcf.gz
PREFIX=${WORKING_DIR}/VCFs/${BATCH_NAME}.${CHR}

printf "\n...extracting biallelic SNPs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${RAW_VCF} \
	-select-type SNP --restrict-alleles-to BIALLELIC \
	-sn $SAMPLELIST -O ${PREFIX}.SNPs.vcf.gz

printf "\n...extracting invariant sites\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${RAW_VCF} \
	-select-type NO_VARIATION \
	-sn $SAMPLELIST -O ${PREFIX}.INVTs.vcf.gz

printf "\n...filtering SNPs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" VariantFiltration -V ${PREFIX}.SNPs.vcf.gz -O ${PREFIX}.SNPs.f.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 40.0" --filter-name "QUAL40" \
	-filter "SOR > 3.0" --filter-name "SOR4" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	-filter "ReadPosRankSum < -12.5" --filter-name "ReadPosRankSum-12.5" \
	-filter "ReadPosRankSum > 12.5" --filter-name "ReadPosRankSum12.5" \
	--verbosity ERROR

printf "\n...filtering INVTs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" VariantFiltration -V ${PREFIX}.INVTs.vcf.gz -O ${PREFIX}.INVTs.f.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "SOR > 3.0" --filter-name "SOR4" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	--verbosity ERROR

printf "\n...selecting passing SNPs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants \
	-V ${PREFIX}.SNPs.f.vcf.gz --exclude-filtered -O ${PREFIX}.SNPs.fp.vcf.gz

printf "\n...selecting passing INVTs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants \
	-V ${PREFIX}.INVTs.f.vcf.gz --exclude-filtered -O ${PREFIX}.INVTs.fp.vcf.gz

# printf "\n...prepare bed file of repeat-masked regions\n" | tee >(cat >&2)
# printf '#chrom\tchromStart\tchromEnd\n' > ${REPEATMASK}.bed
# cut -f1,4,5 ${REPEATMASK} >> ${REPEATMASK}.bed

module purge
ml VCFtools/0.1.16-GCC-13.3.0

printf "\n...filter individual SNP genotypes by depth and GQ, and exclude repeat-masked regions\n" | tee >(cat >&2)
vcftools --gzvcf ${PREFIX}.SNPs.fp.vcf.gz -c --minGQ 15 --minDP 6 --maxDP 100 \
	--bed ${GENE_BED} \
	--recode --recode-INFO-all | bgzip -c > ${PREFIX}.SNPs.fpi.genic.vcf.gz

printf "\n...filter individual INVT genotypes by depth\n" | tee >(cat >&2)
vcftools --gzvcf ${PREFIX}.INVTs.fp.vcf.gz -c --minDP 6 --maxDP 100 \
	--bed ${GENE_BED} \
	--recode --recode-INFO-all | bgzip -c > ${PREFIX}.INVTs.fpi.genic.vcf.gz

module purge
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
ml HTSlib/1.15.1-GCC-11.3.0
ml BCFtools/1.15.1-GCC-11.3.0

tabix -p vcf ${PREFIX}.SNPs.fpi.genic.vcf.gz
tabix -p vcf ${PREFIX}.INVTs.fspi.genic.vcf.gz

printf "\n...merge SNP and INVT sites\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" MergeVcfs \
	-I ${PREFIX}.SNPs.fpi.genic.vcf.gz -I ${PREFIX}.INVTs.fpi.genic.vcf.gz \
	-O ${PREFIX}.merged.fpi.genic.vcf.gz

printf "\n...sort merged VCF\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SortVcf \
	-I ${PREFIX}.merged.fpi.genic.vcf.gz \
	-SD ${GENOME%.fa*}.dict -O ${PREFIX}.merged.fspi.genic.vcf.gz

printf "\n...filter SNPs by mincalled\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${PREFIX}.SNPs.fpi.genic.vcf.gz \
	--max-nocall-number 6 --exclude-filtered -O ${PREFIX}.SNPs.fpi.genic.27called.vcf.gz

printf "\n...filter merged VCF by mincalled\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${PREFIX}.merged.fspi.genic.vcf.gz \
	--max-nocall-number 6 --exclude-filtered -O ${PREFIX}.merged.fspi.genic.27called.vcf.gz

printf "\n...completing all filtering steps\n" | tee >(cat >&2)

printf "\n...making genotype tables\n" | tee >(cat >&2)

printf 'CHROM\tPOS\tREF\tALT\t' > ${PREFIX}.SNPs.fpi.genic.table
printf 'CHROM\tPOS\tREF\tALT\t' > ${PREFIX}.merged.fspi.genic.table

bcftools query -l ${PREFIX}.SNPs.fpi.genic.vcf.gz | tr '\n' '\t' >> ${PREFIX}.SNPs.fpi.genic.table
bcftools query -l ${PREFIX}.merged.fspi.genic.vcf.gz | tr '\n' '\t' >> ${PREFIX}.merged.fspi.genic.table

printf '\n' >> ${PREFIX}.SNPs.fpi.genic.table
printf '\n' >> ${PREFIX}.merged.fspi.genic.table

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ${PREFIX}.SNPs.fpi.genic.vcf.gz >> ${PREFIX}.SNPs.fpi.genic.table
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ${PREFIX}.merged.fspi.genic.vcf.gz >> ${PREFIX}.merged.fspi.genic.table

cat ${PREFIX}.SNPs.fpi.genic.table | sed 's#0/0#0#g' | sed 's#0|0#0#g' | sed 's#0/1#1#g' | sed 's#1/0#1#g' | sed 's#0|1#1#g' | sed 's#1|0#1#g' | sed 's#1/1#2#g' | sed 's#1|1#2#g' | sed 's#./.#NA#g' > ${PREFIX}.SNPs.fpi.genic.gt.table
cat ${PREFIX}.merged.fspi.genic.table | sed 's#0/0#0#g' | sed 's#0|0#0#g' | sed 's#0/1#1#g' | sed 's#1/0#1#g' | sed 's#0|1#1#g' | sed 's#1|0#1#g' | sed 's#1/1#2#g' | sed 's#1|1#2#g' | sed 's#./.#NA#g' > ${PREFIX}.merged.fspi.genic.gt.table

###other steps
python ${SCRIPTS_DIR}/genocounts_groups.py ${POPFILE} ${PREFIX}.merged.fspi.genic.gt.table > ${PREFIX}.merged.fspi.genic.counts.txt
awk '$6<3 && $9<2 && $12<3 && $15<3 {print $1,$2,$3,$4}' ${PREFIX}.merged.fspi.genic.counts.txt > ${PREFIX}.merged.fspi.genic.highqual.sites	###at least 75% called in each subgroup

python ${SCRIPTS_DIR}/genocounts_groups.py ${POPFILE} ${PREFIX}.SNPs.fpi.genic.gt.table > ${PREFIX}.SNPs.fpi.genic.counts.txt
awk '$6<3 && $9<2 && $12<3 && $15<3 {print}' ${PREFIX}.SNPs.fpi.genic.counts.txt > ${PREFIX}.SNPs.fpi.genic.highqual.sites.counts.txt ###at least 75% called in each subgroup
awk '{gutcount=$5+$11; guttotal=$7-$6+$13-$12; gutref=guttotal*2-gutcount; gutprop=gutcount/(2*guttotal); nascount=$8; nastotal=$10-$9; nasprop=nascount/(2*nastotal); nasref=nastotal*2-nascount; if (nasprop-gutprop>=0.8) {print $1,$2,$3,$4,gutref,gutcount,nasref,nascount}}' ${PREFIX}.SNPs.fpi.genic.highqual.sites.counts.txt > ${PREFIX}.SNPs.fpi.genic.AIMs_counts.txt

cut -f1,2,3,4 ${PREFIX}.SNPs.fpi.genic.AIMs_counts.txt > /scratch/mcf96392/AIMs/AIMs_July2025_panel33_${CHR}.AIMs.txt
cut -f1,2,5,6,7,8 ${PREFIX}.SNPs.fpi.genic.AIMs_counts.txt > /scratch/mcf96392/AIMs/AIMs_July2025_panel33_${CHR}.AIMs_counts.txt
cut -f1,2,3,4 ${PREFIX}.SNPs.fpi.genic.highqual.sites.counts.txt > /scratch/mcf96392/AIMs/SNPs_July2025_panel33_${CHR}.txt

# some old number of AIMs: 191921
# some old number of all SNPs: 3493514
