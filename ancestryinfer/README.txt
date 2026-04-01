
#Instructions for running the ancestryinfer pipeline (Matt's version)

### 1. Prepare files:

A. Create modified genome reference files with no "_" characters in the chromosome names
(Do this for both parental genomes)

sed 's/Chr_/Chr/g' /path/to/genome.fa > /path/to/genome.mod.fa

B. Index both modified genomes

C. Prepare/obtain AIMs list: two files, "*.AIMs.fixed.txt" and "*counts.realcMs.fixed.txt"

D. Move/copy all fastq files into your working directory for this project,
and ensure they all have the following exact filename structure:

SAMPLENAME_read_1.fastq.gz
SAMPLENAME_read_2.fastq.gz

E. Create read list: a file with one row per sample, two columns: read_1 filenames and read_2 filenames (no folder paths)

example:
SAMPLE1_read_1.fastq.gz	SAMPLE1_read_2.fastq.gz
NAME_read_1.fastq.gz	NAME_read_2.fastq.gz
OTHERNAME_read_1.fastq.gz	OTHERNAME_read_2.fastq.gz

### 2. Run pre_run_general.sh script

-Quick script, makes some intermediate files that are useful later
-Moves to working directory a modified perl script for the 'samtools' stage of the pipeline (run_samtools_to_hmm_v10_mod.pl)

### 3. Run map_batch_general.sh script

-array script: specify size of array (one per sample) in header or on command line (sbatch --array=1-200%30 map_batch_general.sh)
-this script takes a while (multiple days usually to get through lots of samples)
-to check if it worked, you should have nonempty ...par1.sam and ...par2.sam files for each sample

### 4. Run samtools_batch_general.sh script

-array script: specify size of array (one per sample) in header or on command line (sbatch --array=1-200%30 map_batch_general.sh)
-produces a bunch of intermediate files for each sample
-the important final output for each should be ...sam.hmm.combined.pass.formatted
-takes a while to get through lots of samples (but quicker than mapping step)

### 5. Run hmm_batch_general.sh script

-runs the actual ancestry_HMM
-shouldn't take more than 24 hours even with many samples
-key outputs:
all.indivs.hmm.combined_allchrs.filtered_focalchroms	[the full formatted input file for the hmm]
ancestry-probs_allchrs.tsv_rec.txt	[called genotypes for all individuals/sites]
ancestry-probs_allchrs.tsv_rec.txt_ancestrytransitions_allchrs	[estimated locations of all genotype transitions for each individual]
-also note that the log file for this script has the hmm output statistics (log-likelihood, estimated generations since admixture)

### 6. Run pull_final_files.sh

-moves key output files into a subfolder ANCESTRY_OUTPUTS for easier access
-runs a couple summarizing scripts:
	-summarize_ancestry_bysample.py gives genomewide summary for each individual,
	-summarize_ancestry_bysite.py gives allele/genotype frequencies across all individuals for each site

### Note: to easily transpose (switch rows and columns) output tables, use this script:
curl https://www1.cuni.cz/~obo/textutils/transpose > transpose.pl
