#!/usr/bin/env python3
"""
Convert multiple mpileup from freebayes or methylextract result files into a sync file for popoolation2

 SYNC file format
 2L  79  G   0:0:0:15:0:0    0:0:0:38:0:0
 2L  80  A   12:0:0:0:0:0    38:0:0:0:0:0
 2L  81  A   14:0:0:0:0:0    43:0:0:0:0:0
 2L  82  A   14:0:0:0:0:0    42:0:0:0:0:0
 column 1: reference contig
 column 2: position in the reference contig
 column 3: reference character
 column >3: allele frequencies for all populations in the form A-count:T-count:C-count:G-count:N-count:deletion-count
"""

import argparse
import re

parser = argparse.ArgumentParser(description=None)
parser.add_argument('-o', help='Sync output file name')
parser.add_argument('infiles', nargs='+', help='Freebayes output files (min 2)')

args = parser.parse_args()

alldata = {}
allfiles = []
filetype = ''

# read in all files
for fname in args.infiles:
    filename = fname.split('/')[-1].split('.')[0]
    allfiles.append(filename)
    for line in open(fname):
        # Check for the file format
        if line[0] == '#':
            if line[0:18] == "##source=freeBayes":
                filetype = 'fb'
            elif line[0:22] == "##source=MethylExtract":
                filetype = 'me'
            continue
        line = line.strip('\n').split('\t')
        chrom = line[0]
        pos = line[1]
        if len(line) == 6:
            RB = line[2]
            AB = line[3]
            varnames = line[4].split(':')
            varvalues = line[5].split(':')
        else:
            RB = line[3]
            AB = line[4]
            varnames = line[8].split(':')
            varvalues = line[9].split(':')
        if chrom not in alldata:
            alldata[chrom] = {}
        if pos not in alldata[chrom]:
            alldata[chrom][pos] = {}
        data = {}
        for k, v in zip(varnames, varvalues):
            data[k] = v
        if filetype == "fb":
            alldata[chrom][pos][filename] = [RB, AB, data['RO'], data['AO']]
        else:
            af = re.findall("AF=([^;]+)", line[7])
            af = af[0].split(',')
            alldata[chrom][pos][filename] = [RB, AB, data['DP4'].split(','), af]


# Print out the sync file
fo = open(args.o, "w")
for chrom in alldata:
    for pos in alldata[chrom]:
        snpstrings = []
        b2p = {"A": 0, "T": 1, "C": 2, "G": 3, "N": 4}
        for ff in allfiles:
            refstring = ["0", "0", "0", "0", "0", "0"]
            if ff in alldata[chrom][pos]:
                mydata = alldata[chrom][pos][ff]
                refbase = mydata[0]
                altbases = mydata[1].split(',')
                if refbase != "." and len(refbase) == 1:
                    if filetype == 'fb':
                        refstring[b2p[refbase]] = mydata[2]
                        altcounts = mydata[3].split(',')
                        i = 0
                        for abase in altbases:
                            if abase == "N":
                                continue
                            if abase != "." and len(abase) == 1:
                                refstring[b2p[abase]] = altcounts[i]
                            i += 1
                    else:
                        refstring[b2p[refbase]] = str(int(float(mydata[2][0]))+int(float(mydata[2][1])))
                        altcounts = int(int(float(mydata[2][0]))+int(float(mydata[2][1]))+float(mydata[2][2]))+int(float(mydata[2][3]))
                        for abase, afactor in zip(altbases, mydata[3]):
                            if len(abase) == 1:
                                refstring[b2p[abase]] = str(int(round(altcounts*float(afactor))))
            rstring = ":".join(refstring)
            snpstrings.append(rstring)
        fo.write("\t".join([chrom, pos, refbase, "\t".join(snpstrings)]) + "\n")
fo.close()
