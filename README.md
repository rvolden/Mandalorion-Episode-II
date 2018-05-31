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
python3 defineAndQuantifyWrapper.py -c content_file -p /path/to/data/ -u 5 -d 30 -s 200 -r 0.05 -R 3 -i 0 -t 0 -I 100 -T 60 -g gencode.v26.annotation.gtf -m NUC.4.4.mat -f example_config
```
