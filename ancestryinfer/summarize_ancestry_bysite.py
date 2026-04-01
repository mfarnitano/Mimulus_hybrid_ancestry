#!/usr/bin/python

#script to get allele frequencies from ancestry call table
#takes in transposed call table, with columns for SNPs and rows for samples

### input file format
#header: #sample[\tchrom_position]...
#body:	$samplename[\tgenotype]... #genotypes are 0,1,2 or NA
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
    for line in input:
        if start:
            header=["Site"]+line.strip().split()[1:]
            nfields=len(header)-1
            counts=[0]*nfields #count alt alleles
            totals=[0]*nfields #count total alleles
            hets=[0]*nfields #count het genotypes
            start=False
        else:
            fields=line.strip().split()
            genotypes=fields[1:]
            for i,g in enumerate(genotypes):
                if g!="NA":
                    counts[i]=counts[i]+int(g)
                    totals[i]=totals[i]+2
                    if g=="1":
                        hets[i]=hets[i]+1
counts=[str(x) for x in counts]
totals=[str(x) for x in totals]
hets=[str(x) for x in hets]

output.write('\t'.join(header)+'\n')
output.write("P2_alleles\t"+'\t'.join(counts)+'\n')
output.write("Total_alleles\t"+'\t'.join(totals)+'\n')
output.write("Het_genotypes\t"+'\t'.join(hets)+'\n')
#done
