#!/usr/bin/python


### takes hmm outputs (or similar genotype call table) and summarizes a single genotype per window
### from ancestryinfer pipeline, use input 'ancestry-probs_allchrs.tsv_rec.txt_transposed'

###input: rows are sites, columns are samples (with a header row), each value is 0,1,2, or NA
###sites are listed as chr:pos in the first column
###output is one window per line, 0,1,2, or NA
###current strategy: if a proportion ("cutoff") of non-NA sites agree, make a call, otherwise set to NA
###optional minimum number of agreeing sites: "minsupport", otherwise NA

import sys
import math

calltable=sys.argv[1]

if len(sys.argv) > 2:
    output=open(sys.argv[2],'w')
else:
    output=sys.stdout

windowsize=50000
cutoff=0.9
minsupport=5

#starting null data
nsamples=None
current_genotypes={}
current_chr=None
current_window=1
current_windowend=current_window+windowsize-1
nsites=0

with open(calltable,'r') as calls:

    for line in calls:
        linelist=line.strip().split()
        nsamples=len(linelist)-1
        chrpos=linelist[0]
        #read and write header line
        if chrpos=="chrom_position" or chrpos=="id":    #if header row starts differently, add an option here
            current_genotypes={"0":[0]*nsamples,"1":[0]*nsamples,"2":[0]*nsamples,"NA":[0]*nsamples}
            writestring="chrom:windowstart:windowend\tnsites\t"
            writecols='\t'.join(linelist[1:])
            output.write(writestring+writecols+"\n")
            continue
        gts=linelist[1:]
        chr=chrpos.split(':')[0]
        pos=int(chrpos.split(':')[1])

        if current_chr is None:
            #start first window of first chrom
            current_chr=chr
            current_windowend=math.ceil(pos/windowsize)*windowsize
            current_window=current_windowend-windowsize+1
            for i,g in enumerate(gts):
                current_genotypes[g][i]+=1
            nsites=1

        elif (chr != current_chr) or (pos > current_windowend):
            #write existing window
            write_gts=[]
            for i in range(nsamples):
                counted=current_genotypes["0"][i]+current_genotypes["1"][i]+current_genotypes["2"][i]
                if counted==0:
                    write_gts.append("NA")
                    continue
                aa=current_genotypes["0"][i]/counted
                ab=current_genotypes["1"][i]/counted
                bb=current_genotypes["2"][i]/counted
                if aa>=cutoff and current_genotypes["0"][i]>=minsupport:
                    write_gts.append("0")
                elif ab>=cutoff and current_genotypes["1"][i]>=minsupport:
                    write_gts.append("1")
                elif bb>=cutoff and current_genotypes["2"][i]>=minsupport:
                    write_gts.append("2")
                else:
                    write_gts.append("NA")
            writestring=current_chr+':'+str(current_window)+':'+str(current_windowend)+'\t'+str(nsites)+'\t'
            writecols='\t'.join(write_gts)
            output.write(writestring+writecols+'\n')
            #start new window
            current_chr=chr
            current_windowend=math.ceil(pos/windowsize)*windowsize
            current_window=current_windowend-windowsize+1
            current_genotypes={"0":[0]*nsamples,"1":[0]*nsamples,"2":[0]*nsamples,"NA":[0]*nsamples}
            for i,g in enumerate(gts):
                current_genotypes[g][i]+=1
            nsites=1
        else:
            for i,g in enumerate(gts):
                current_genotypes[g][i]+=1
            nsites+=1

    ###make sure to print last window
    write_gts=[]
    for i in range(nsamples):
        counted=current_genotypes["0"][i]+current_genotypes["1"][i]+current_genotypes["2"][i]
        if counted==0:
            write_gts.append("NA")
            continue
        aa=current_genotypes["0"][i]/counted
        ab=current_genotypes["1"][i]/counted
        bb=current_genotypes["2"][i]/counted
        if aa>=cutoff and current_genotypes["0"][i]>=minsupport:
            write_gts.append("0")
        elif ab>=cutoff and current_genotypes["1"][i]>=minsupport:
            write_gts.append("1")
        elif bb>=cutoff and current_genotypes["2"][i]>=minsupport:
            write_gts.append("2")
        else:
            write_gts.append("NA")
    writestring=current_chr+':'+str(current_window)+':'+str(current_windowend)+'\t'+str(nsites)+'\t'
    writecols='\t'.join(write_gts)
    output.write(writestring+writecols+'\n')
##DONE
