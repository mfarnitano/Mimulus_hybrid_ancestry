#!/usr/bin/python

#script to get allele frequencies from ancestry call table
#takes in transposed call table, with columns for SNPs and rows for samples

### input file format
#header: #sample[\tchrom_position]...
#body:	$samplename[\tgenotype]... #genotypes are 0,1,2 or NA

#usage: python ./summarize_ancestry_bysample.py [input] [output]
#if no output filename given, prints to stdout

import sys
if len(sys.argv) > 2:
    output=open(sys.argv[2],'w')
else:
    output=sys.stdout
header=[]
counts=[]
totals=[]
hets=[]
with open(sys.argv[1],'r') as input:
    start=True
    output.write("sample\tcounted_sites\tmissing_sites\thet_sites\tp2_hom_sites\tp2_alleles\tproportion_p2\theterozygosity\n")
    for line in input:
        data=line.strip().split()
        if start:
            start=False
            continue
        else:
            samplename=data[0]
            genotypes=data[1:]
            counted_sites=0
            missing_sites=0
            het_sites=0
            p2_hom_sites=0
            p2_alleles=0
            for g in genotypes:
                try:
                    geno=int(g)
                    counted_sites+=1
                    p2_alleles+=geno
                    if geno==1:
                        het_sites+=1
                    elif geno==2:
                        p2_hom_sites+=1
                except ValueError:
                    missing_sites+=1
            try:
                proportion_p2=p2_alleles/(2*counted_sites)
                heterozygosity=het_sites/counted_sites
            except ZeroDivisionError:
                proportion_p2="NA"
                heterozygosity="NA"
            output.write(samplename+'\t'+str(counted_sites)+'\t'+str(missing_sites)+'\t'+str(het_sites)+'\t'+str(p2_hom_sites)+'\t'+str(p2_alleles)+'\t'+str(proportion_p2)+'\t'+str(heterozygosity)+'\n')

if output is not sys.stdout:
    output.close()
#DONE
