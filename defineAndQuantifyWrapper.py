#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--content_file', type=str)
parser.add_argument('-p', '--path', type=str)
#parser.add_argument('-a','--genome_annotation',type=str)
#parser.add_argument('-g','--gmap_genome',type=str)
#parser.add_argument('-l','--gene_list',type=str)
#parser.add_argument('-i','--illumina_content_file',type=str, default='-')
#parser.add_argument('-r','--refine',type=str,default=None,choices=['g','gi','i','n'])
# if set to 'i' or 'gi' and no illumina_content_file is provided, the i will be ignored

parser.add_argument('-u', '--upstream_buffer', type=str)
parser.add_argument('-d', '--downstream_buffer', type=str)
parser.add_argument('-s', '--subsample_consensus', type=str)
parser.add_argument('-g', '--genome_annotation', type=str)
parser.add_argument('-r', '--minimum_ratio', type=str)
parser.add_argument('-R', '--minimum_reads', type=str)
parser.add_argument('-i', '--minimum_5_overhang', type=str)
parser.add_argument('-t', '--minimum_3_overhang', type=str)
parser.add_argument('-I', '--maximum_5_overhang', type=str)
parser.add_argument('-T', '--maximum_3_overhang', type=str)

args = parser.parse_args()

content_file = args.content_file         # file containing paths to psl and fasta/q files you want analyzed
path = args.path + '/'	         #path where you want your output files to go
upstream_buffer = args.upstream_buffer
downstream_buffer = args.downstream_buffer
subsample_consensus = args.subsample_consensus
genome_annotation = args.genome_annotation
minimum_ratio = args.minimum_ratio
minimum_reads = args.minimum_reads
minimum_5_overhang = args.minimum_5_overhang
minimum_3_overhang = args.minimum_3_overhang
maximum_5_overhang = args.maximum_5_overhang
maximum_3_overhang = args.maximum_3_overhang

os.system('python3 spliceSites.py %s %s %s %s %s ' %(content_file, path, '0.05', genome_annotation, 'g'))
os.system('python3 defineAndQuantifyIsoforms.py %s %s %s %s' %(content_file, path,downstream_buffer,upstream_buffer)) # This script sort raw reads into isoform bins. The two number variables determine the window around TSS and TES in which read ends can fall and still be matched to the site.

for line in open(content_file):
    subpath = line.strip().split('\t')[2]
    os.system('python3 createConsensi.py %s %s ' %(subpath,subsample_consensus))
    os.system('python3 filterIsoforms.py %s %s %s %s %s %s %s %s ' %(subpath,subpath+'/Isoform_Consensi.fasta',minimum_ratio,minimum_reads,maximum_5_overhang,maximum_3_overhang,minimum_5_overhang,minimum_3_overhang))
