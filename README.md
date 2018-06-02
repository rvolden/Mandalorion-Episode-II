# Mandalorion Episode II #
Attack of the Isoforms

Takes R2C2/C3POa data and defines high confidence isoforms. 

We plan to ultimately support straight up ONT 1D and 1D2 reads or PacBio reads but for now reads will have to be in the format C3POa produces

C3POa produces fasta file containing consensus reads as well as a fastq file containing subreads. 
You will need both of those files as well as alignments of the consensus reads to a genome. We recommend you use minimap2 and then convert the sam file it produces to psl format using the sam2psl tool from jvarkit.

- [sam2psl](http://lindenb.github.io/jvarkit/SamToPsl.html)

Mandalorion takes read alignments and performs the following steps.

1.) It parses a genome annotation file (gtf) to extract all annotated splice sites.

2.) It identifies un-annotated splice sites based on read alignments of all the samples in the content file

3.) It analyzes every individual sample and groups read alignments based on splice site usage.

4.) It identifies all TSS and polyA sites used by these grouped read alignments (often multiple per group, doesn't rely on annotation).

5.) It defines isoforms based on splice site, TSS and polyA site usage

6.) It generates a consensus sequence for each isoform

7.) It filters all consensus sequences based on read number and un-aligned bases at their ends

The main output file is the Isoform_Consensi_filtered.fasta file that is generated for every sample in the content file

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


Example command:
```
python3 defineAndQuantifyWrapper.py -c content_file -p /path/to/data/ -u 5 -d 30 -s 200 -r 0.05 -R 3 -i 0 -t 0 -I 100 -T 60 -g gencode.v26.annotation.gtf -m NUC.4.4.mat -f example_config
```
