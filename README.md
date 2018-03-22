# Mandalorion-Episode-II #
Attack of the Isoforms

Takes R2C2/C3POa data and defines/quantifies isoforms.
Consensi can be made from the isoform clusters.

## Inputs ##
* takes 5-3 trimmed reads
* subreads
* alignment (reads to genome)
* gtf for annotations

Example command:
```
python3 Mandalorion_define_and_quantify_isoforms.py -c /home/vollmers/data/R2C2_RAD4_CD27CD38_Combined/content_file -p /home/vollmers/data/R2C2_RAD4_CD27CD38_Combined/ -u 5 -d 30 -s 200 -r 0.05 -R 3 -i 0 -t 0 -I 100 -T 60 -g /home/vollmers/Genomes/gencode.v26.annotation.gtf
```
