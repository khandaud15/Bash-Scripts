#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=4
#PBS -N rseqc
#PBS -j oe
# PBS -o /home/glasser/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR


geneBody_coverage.py -i 03_sorted/CEC5_S7_L002_R1_001.bam,03_sorted/CEC6_S8_L002_R1_001.bam -r /home/mk9uc/sequencing/data/ExtData/UCSC/mm10/Mus_musculus.GRCm38.83.bed -o RseQC_geneBody_coverage
