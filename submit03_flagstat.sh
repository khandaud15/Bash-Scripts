#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -M my@email.adr
#PBS -N flagstat
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

#cd $PBS_O_WORKDIR
module load samtools
# Input folder
DIRIN= /sfs/qumulo/qproject/cphg-millerlab/CAD_QTL/coronary_QTL/epigenome/atac/bam/


# Output file
FOUT=flagstat.txt

for file in `find $DIRIN -type f -name "*.bam" | sort`; do
		echo $file >> $FOUT;
        samtools flagstat $file >> $FOUT;
done
