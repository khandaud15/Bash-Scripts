#!/bin/bash

module load samtools
DIRIN=/sfs/qumulo/qproject/CPHG/MILLER/CAD_QTL/coronary_QTL/expression/pass1


# Output file
FOUT=flagstat.txt

for file in `find $DIRIN -type f -name "*.bam" | sort`; do
                echo $file >> $FOUT;
        samtools flagstat $file >> $FOUT;
done
