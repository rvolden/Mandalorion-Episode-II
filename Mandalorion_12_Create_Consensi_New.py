#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import os
import numpy as np

path=sys.argv[1]
temp_folder=path+'/parsed_reads'
subsample=int(sys.argv[2])
#os.system('mkdir '+temp_folder)


poa = 'poa'
score_matrix = '/home/vollmers/scripts/ONT/NUC.4.4.mat'
  # Change this if you want to subsample more or less reads per isoform

racon='/home/vollmers/Downloads/racon/bin/racon'
minimap2='/home/vollmers/Downloads/minimap2-2.5_x64-linux/minimap2'
consensus='/home/vollmers/scripts/ONT/R2C2/consensus.py'


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
    read_list = []
    length = 0
    for line in open(seq_file):
        length += 1
    lineNum = 0
    seq_file_open = open(seq_file, 'r')
    while lineNum < length:
        name_root = seq_file_open.readline().strip()[1:].split('_')
        name, seed = name_root[0], int(name_root[1])
        seq = seq_file_open.readline().strip()
        plus = seq_file_open.readline().strip()
        qual = seq_file_open.readline().strip()
        quals = []
        for character in qual:
            number = ord(character) - 33
            quals.append(number)
        average_quals = np.average(quals)
        seq_length = len(seq)
        read_list.append((name, seed, seq, qual, average_quals, seq_length))

        lineNum += 4
    return read_list


def read_fasta(infile):
    reads={}
    sequence=''
    name=''
    for line in open(infile):
        a=line.strip()
        if a[0]=='>':
            if sequence!='':

                reads[name]=sequence
            name=a[1:].split()[0]
            sequence=''
        else:
            sequence+=a
    reads[name]=sequence
    return reads

def reverse_complement(sequence):
  Seq=''
  complement = {'A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
  for item in sequence[::-1]:
    Seq=Seq+complement[item]
  return Seq


def determine_consensus(name, fasta,fastq):
    '''Aligns and returns the consensus'''
    print('determine')
    corrected_consensus = ''

    out_F = fasta

    fastq_reads=read_fastq_file(fastq)

    out_Fq=temp_folder + '/subsampled.fastq'

    out=open(out_Fq,'w')


    indexes=np.random.choice(np.arange(0,len(fastq_reads),1),min(len(fastq_reads),subsample))
    print('XXXXXXXXXXXXXXXXXX',indexes)
    subsample_fastq_reads=[]
    for index in indexes:
        subsample_fastq_reads.append(fastq_reads[index])

    for read in subsample_fastq_reads:

        out.write('@'+read[0]+'_'+str(read[1])+'\n'+read[2]+'\n+\n'+read[3]+'\n')

    out.close()

    poa_cons = temp_folder + '/consensus.fasta'
    final = temp_folder + '/corrected_consensus.fasta'
    overlap = temp_folder +'/overlaps.sam'
    pairwise = temp_folder + '/prelim_consensus.fasta'


#    PIR = temp_folder + '/' + name + 'alignment.fasta'
#    os.system('%s -read_fasta %s -hb -pir %s \
#               %s >./poa_messages.txt 2>&1' \
#              %(poa, out_F, PIR, score_matrix))
#    reads = read_fasta(PIR)
#    print('poa done')


#    if repeats == '2':
#        Qual_Fasta = open(pairwise, 'w')
#        for read in reads:
#            if 'CONSENS' not in read:
#                Qual_Fasta.write('>' + read + '\n' + reads[read] + '\n')
#        Qual_Fasta.close()
#        os.system('%s %s %s %s >> %s' \
#                  %(consensus, pairwise, out_Fq, name, poa_cons))

#    else:
#        for read in reads:
#          if 'CONSENS0' in read:
#            out_cons_file = open(poa_cons, 'w')
#            out_cons_file.write('>' + name + '\n' \
#                                + reads[read].replace('-', '') + '\n')
#            out_cons_file.close()
#    print('cosensus done')

    max_coverage=0
    reads=read_fasta(out_F)
    repeats=str(len(reads))
    for read in reads:
        coverage=int(read.split('_')[3])
        if coverage>max_coverage:
             best=read
             max_coverage=coverage

    out_cons_file = open(poa_cons, 'w')
    out_cons_file.write('>' + best + '\n' + reads[best].replace('-', '') + '\n')
    out_cons_file.close()


    final=poa_cons
    for i in np.arange(1,2,1):
        try:
            if i==1:
                input_cons=poa_cons
                output_cons=poa_cons.replace('.fasta','_'+str(i)+'.fasta')
            else:
                input_cons=poa_cons.replace('.fasta','_'+str(i-1)+'.fasta')
                output_cons=poa_cons.replace('.fasta','_'+str(i)+'.fasta')

            print(i,minimap2,input_cons,out_Fq,overlap)
            os.system('%s --secondary=no -ax map-ont \
                      %s %s > %s 2> ./minimap2_messages.txt' \
                      % (minimap2, input_cons, out_Fq, overlap))
            print('minimap2 done')
            os.system('%s --sam --bq 5 -t 1 \
                      %s %s %s %s >> ./racon_messages.txt 2>&1' \
                      %(racon,out_Fq, overlap, input_cons, output_cons))
            print('racon done')
            final=output_cons
        except:
            pass

    print(final)
    reads = read_fasta(final)
#    if len(reads)==0:
#        reads = read_fasta(poa_cons)

    for read in reads:
        corrected_consensus = reads[read]

    return corrected_consensus, repeats



combined_consensus_file=open(path+'/Isoform_Consensi.fasta','w')
combined_consensus_file.close()


determine_consensus



for line in open(path+'/isoform_list'):
    fasta=line.split('\t')[0]
    fastq=line.split('\t')[1]
    name=line.split('\t')[2].strip()
    print(line)
    print(fasta,fastq,name)
    corrected_consensus,repeats=determine_consensus(name,fasta,fastq)
    combined_consensus_file=open(path+'/Isoform_Consensi.fasta','a')

    combined_consensus_file.write('>'+name+'_'+repeats+'\n'+corrected_consensus+'\n')
    combined_consensus_file.close()
