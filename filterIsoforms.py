#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import os
import numpy as np

path = sys.argv[1]
infile = sys.argv[2]
minimum_ratio = float(sys.argv[3])
minimum_reads = int(sys.argv[4])
maximum_5_overhang = float(sys.argv[5])
maximum_3_overhang = float(sys.argv[6])
minimum_5_overhang = float(sys.argv[7])
minimum_3_overhang = float(sys.argv[8])

out = open(path + '/Isoform_Consensi_filtered.fasta', 'w')

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            readDict[line[1:]] = ''
            lastHead = line[1:]
        else:
            readDict[lastHead] += line
    return readDict

isoforms = read_fasta(infile)

count = {}
for isoform in isoforms:
    parts = isoform.split('_')
    chromosome = parts[0]
    direction = parts[2]
    start = int(parts[3])
    end = int(parts[4])
    number = int(parts[-1])

    for base in np.arange(start, end):
        if not count.get(chromosome + '_' + direction):
            count[chromosome + '_' + direction] = {}
        if not count[chromosome + '_' + direction].get(base):
            count[chromosome + '_' + direction][base] = number
        else:
            count[chromosome + '_' + direction][base] += number

for isoform in isoforms:
    parts = isoform.split('_')
    chromosome = parts[0]
    direction = parts[2]
    start = int(parts[3])
    end = int(parts[4])
    number = int(parts[-1])
    coverage_list = []
    overhang5 = float(parts[5])
    overhang3 = float(parts[6])
    for base in np.arange(start, end):
        coverage_list.append(count[chromosome + '_' + direction][base])

    max_coverage = max(coverage_list)
    if minimum_3_overhang <= overhang3 <= maximum_3_overhang \
       and minimum_5_overhang <= overhang5 <= maximum_5_overhang:
        if number >= minimum_reads \
           and number/max(coverage_list) >= minimum_ratio:
            out.write('>' + isoform + '_'
                      + str(round(number/max(coverage_list), 2))
                      + '\n' + isoforms[isoform] + '\n')
