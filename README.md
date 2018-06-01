# Mandalorion-Episode-II #
Attack of the Isoforms

Takes R2C2/C3POa data and defines/quantifies isoforms.
Consensi can be made from the isoform clusters.

## Dependencies ##

- [minimap2 2.7-r654](https://github.com/lh3/minimap2)
- [racon](https://github.com/isovic/racon)


## Inputs ##

-c, --content_file

tab-delimited file containing information for all samples you want to analyze.
Each sample should have in one line (separated by tabs):
/path/to/alignments.psl /path/to/reads.fasta /path/to/output/ /path/to/subreads.fastq 

-f, --config_file

text file containing paths to minimap2 and racon in the following format:

minimap2[tab]/path/to/minimap2
racon[tab]/path/to/racon

-p, --path

This is where the bed file containing splice sites will be saved.

-m, --score_matrix

path/to/NUC.4.4.mat 

-u, --upstream_buffer and -d, --downstream_buffer

These values define how lenient TSS and polyA sites are defined. We use an -u of 5 and -d of 30 meaning that read end aligning within 5nt upstream and 30nt downstream of a TSS or polyA site are combined.

-s, --subsample_consensus

This number defines how many randomly sampled subreads are used to create a Isoform consensus. We use 200 here. 

-g, --genome_annotation

Genome annotation of choice in gtf format. We have only used Gencode annotations and the gtf format is not very well defined so be careful when using other annotations. 

-r, --minimum_ratio

Minimum ratio of reads aligned to a locus that have to be assigned to an isoform for the isoform to be reported. We use 0.05

-R, --minimum_reads

Minimum number of reads that have to be assigned to an isoform for the isoform to be reported. We use 3

-i, --minimum_5_overhang and -I, --maximum_5_overhang

Only isoforms with median number of un-aligned bases on the 5' end between these numbers are reported. We use 0 and 100 but this numbers may be different based on your cDNA adapters 

-t, --minimum_3_overhang and -T, --maximum_3_overhang

Only isoforms with median number of un-aligned bases on the 3' end between these numbers are reported. We use 0 and 60 but this numbers may be different based on your cDNA adapters


* takes 5'-3' oriented full length cDNA reads
* subreads
* alignment (reads to genome)
* gtf for annotations

Example command:
```
python3 defineAndQuantifyWrapper.py -c content_file -p /path/to/data/ -u 5 -d 30 -s 200 -r 0.05 -R 3 -i 0 -t 0 -I 100 -T 60 -g gencode.v26.annotation.gtf -m NUC.4.4.mat -f example_config
```
