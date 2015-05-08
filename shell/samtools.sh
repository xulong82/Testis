#!/bin/bash

module load samtools

dir1="/data/xwang/Testis/RSEM"
dir2="/home/xwang/Dropbox/GitHub/Testis/shiny/gviz/bam"

files=`find $dir1 -name '*.transcript.sorted.bam'`

cd $dir1

for name1 in $files; do
  name2=`basename $name1`
  name3=${name2/.transcript.sorted.bam/}
  echo $name3
# samtools view -b -q 1 $name2 \
  samtools view -b $name2 \
	ENSMUST00000084214 \
	ENSMUST00000084215 \
	ENSMUST00000105830 \
	ENSMUST00000105831 > $dir2/$name3.bam
done
  
