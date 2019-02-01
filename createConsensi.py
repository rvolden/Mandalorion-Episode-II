#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import os
import argparse
import numpy as np

def argParser():
    parser = argparse.ArgumentParser(description = 'Makes consensus sequences \
                                                    from R2C2 reads.',
                                     add_help = True,
                                     prefix_chars = '-')
    parser.add_argument('--path', '-p', type=str, action='store', default=os.getcwd(),
                        help='Directory where all the files are/where they will end up.\
                              Defaults to your current directory.')
    parser.add_argument('--subsample', '-s', type=int, action='store')
    parser.add_argument('--config', '-c', type=str, action='store', default='',
                        help='If you want to use a config file to specify paths to\
                              programs, specify them here. Use for poa, racon, water,\
                              blat, and minimap2 if they are not in your path.')
    parser.add_argument('--matrix', '-m', type=str, action='store',
                        default='NUC.4.4.mat',
                        help='Score matrix to use for poa.\
                              Defaults to NUC.4.4.mat.')
    return vars(parser.parse_args())

def configReader(configIn):
    '''Parses the config file.'''
    progs = {}
    for line in open(configIn):
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t')
        progs[line[0]] = line[1]
    # should have minimap, poa, racon, water, consensus
    # check for extra programs that shouldn't be there
    possible = set(['poa', 'minimap2', 'water', 'consensus', 'racon', 'blat'])
    inConfig = set()
    for key in progs.keys():
        inConfig.add(key)
        if key not in possible:
            raise Exception('Check config file')
    # check for missing programs
    # if missing, default to path
    for missing in possible-inConfig:
        if missing == 'consensus':
            path = 'consensus.py'
        else:
            path = missing
        progs[missing] = path
        sys.stderr.write('Using ' + str(missing)
                         + ' from your path, not the config file.\n')
    return progs

args = argParser()
path = args['path']
temp_folder = path + '/parsed_reads'
subsample = args['subsample']
score_matrix = args['matrix']

if args['config']:
    progs = configReader(args['config'])
    minimap2 = progs['minimap2']
    poa = progs['poa']
    racon = progs['racon']
    water = progs['water']
    consensus = progs['consensus']
else:
    minimap2, poa, racon, water = 'minimap2', 'poa', 'racon', 'water'
    consensus = 'consensus.py'

consensus = 'python3 ' + consensus

def read_fastq_file(seq_file):
    '''
    Takes a FASTQ file and returns a list of tuples
    In each tuple:
        name : str, read ID
        seed : int, first occurrence of the splint
        seq : str, sequence
        qual : str, quality line
        average_quals : float, average quality of that line
        seq_length : int, length of the sequence
    '''
    read_list, lineNum = [], 0
    lastPlus = False
    for line in open(seq_file):
        line = line.rstrip()
        if not line:
            continue
        # make an entry as a list and append the header to that list
        if lineNum % 4 == 0 and line[0] == '@':
            splitLine = line[1:].split('_')
            root, seed = splitLine[0], int(splitLine[1])
            read_list.append([])
            read_list[-1].append(root)
            read_list[-1].append(seed)

        # sequence
        if lineNum % 4 == 1:
            read_list[-1].append(line)

        # quality header
        if lineNum % 4 == 2:
            lastPlus = True

        # quality
        if lineNum % 4 == 3 and lastPlus:
            read_list[-1].append(line)
            avgQ = sum([ord(x)-33 for x in line])/len(line)
            read_list[-1].append(avgQ)
            read_list[-1].append(len(read_list[-1][2]))
            read_list[-1] = tuple(read_list[-1])

        lineNum += 1
    return read_list

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

def determine_consensus(name, fasta, fastq):
    '''Aligns and returns the consensus'''
    print('determine')
    corrected_consensus = ''
    out_F = fasta
    fastq_reads = read_fastq_file(fastq)
    out_Fq = temp_folder + '/subsampled.fastq'
    out = open(out_Fq, 'w')

    indexes = np.random.choice(np.arange(0, len(fastq_reads)),
                               min(len(fastq_reads), subsample), replace=False)
    print('XXXXXXXXXXXXXXXXXX', indexes)
    subsample_fastq_reads = []
    for index in indexes:
        subsample_fastq_reads.append(fastq_reads[index])

    for read in subsample_fastq_reads:
        out.write('@' + read[0] + '_' + str(read[1]) + '\n'
                  + read[2] + '\n+\n' + read[3] + '\n')
    out.close()

    poa_cons = temp_folder + '/consensus.fasta'
    final = temp_folder + '/corrected_consensus.fasta'
    overlap = temp_folder +'/overlaps.sam'
    pairwise = temp_folder + '/prelim_consensus.fasta'

    max_coverage = 0
    reads = read_fasta(out_F)
    repeats = str(len(reads))
    for read in reads:
        coverage = int(read.split('_')[3])
        if coverage > max_coverage:
             best = read
             max_coverage = coverage

    out_cons_file = open(poa_cons, 'w')
    out_cons_file.write('>' + best + '\n'
                        + reads[best].replace('-', '') + '\n')
    out_cons_file.close()


    final = poa_cons
    for i in np.arange(1, 2):
        try:
            if i == 1:
                input_cons = poa_cons
                output_cons = poa_cons.replace('.fasta',
                                               '_' + str(i) + '.fasta')
            else:
                input_cons = poa_cons.replace('.fasta',
                                              '_'+str(i-1) + '.fasta')
                output_cons = poa_cons.replace('.fasta',
                                               '_' + str(i) + '.fasta')

            print(i, minimap2, input_cons, out_Fq, overlap)
            os.system('%s --secondary=no -ax map-ont\
                      %s %s > %s 2> ./minimap2_messages.txt'
                      % (minimap2, input_cons, out_Fq, overlap))
            print('minimap2 done')
            os.system('%s -q 5 -t 1 \
                      %s %s %s >%s 2> ./racon_messages.txt'
                      %(racon,out_Fq, overlap, input_cons, output_cons))
            print('racon done')
            final = output_cons
        except:
            pass

    print(final)
    reads = read_fasta(final)
#    if len(reads)==0:
#        reads = read_fasta(poa_cons)
    for read in reads:
        corrected_consensus = reads[read]
    return corrected_consensus, repeats

combined_consensus_file = open(path + '/Isoform_Consensi.fasta', 'w')
combined_consensus_file.close()

for line in open(path + '/isoform_list'):
    fasta = line.split('\t')[0]
    fastq = line.split('\t')[1]
    name = line.split('\t')[2].strip()
    print(line)
    print(fasta, fastq, name)
    corrected_consensus, repeats = determine_consensus(name, fasta, fastq)
    combined_consensus_file = open(path + '/Isoform_Consensi.fasta', 'a')

    combined_consensus_file.write('>' + name + '_' + repeats + '\n'
                                  +corrected_consensus + '\n')
    combined_consensus_file.close()
