#!/bin/bash

dir1="/data/xwang/Testis/FASTQ1"
dir2="~/testis"

files=`find $dir1 -name '*LaneALL_R1.fastq'`

for name1 in $files; do
  name2=`basename $name1`
  name3=${name2/_R1.fastq/}
  echo $name3
  qsub -v arg=$name3 $dir2/bowtie.sh
# qsub -v arg=$name3 $dir2/rsem.sh
# qsub -v arg=$name3 $dir2/tophat.sh
done
  
# files=`find $dir2 -name '*.sh.e*'`

# for name1 in $files; do
#   name2=`basename $name1`
#   echo $name2 >> 1
#   grep "at least" $name2 >> 1
# done 

