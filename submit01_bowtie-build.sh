#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -N btbuild
#PBS -j oe
# PBS -o /home/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR

bowtie-build2 /home/sequencing/data/ExtData/UCSC/hg38/hg38.fa /home/sequencing/data/ExtData/UCSC/hg38/hg38 
# /local/bin/bowtie2-build /home/sequencing/data/ExtData/UCSC/mm9/mm9.fa /home/sequencing/data/ExtData/UCSC/mm9/mm9
# /local/bin/bowtie2-build /home/sequencing/data/ExtData/UCSC/mm10/mm10.fa /home/sequencing/data/ExtData/UCSC/mm10/mm10

