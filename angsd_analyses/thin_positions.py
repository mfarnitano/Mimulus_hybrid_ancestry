#!/usr/bin/python


###script to randomly thin a list of sites to 1 per (bin length)

###USAGE: python ./thin_positions.py [bin_length] [input] > [output]
###input should be two columns, chr and pos

import sys
import math
import random

binlength=int(sys.argv[1])

if len(sys.argv) > 3:
    output=open(sys.argv[3],'w')
else:
    output=sys.stdout
with open(sys.argv[2],'r') as input:
    current_chr=""
    current_bin=0
    bin_contents=[]
    first=True

    for line in input:
        cols=line.strip().split()
        try:
            chr,pos=cols[0],int(cols[1])
        except ValueError:
            #header or nonviable row to skip
            continue
        if first:
            current_chr=chr
            current_bin=math.floor(pos/binlength)*binlength
            first=False
        bin=math.floor(pos/binlength)*binlength
        if current_chr == chr and bin == current_bin:
            #continue same bin
            bin_contents.append((chr,pos))
        else:
            #end bin
            kept=random.choice(bin_contents)
            output.write(kept[0]+"\t"+str(kept[1])+"\n")

            #start new bin
            current_chr=chr
            current_bin=math.floor(pos/binlength)*binlength
            bin_contents=[(chr,pos)]
    #last bin
    kept=random.choice(bin_contents)
    output.write(kept[0]+"\t"+str(kept[1])+"\n")

if len(sys.argv) > 3:
    output.close()
