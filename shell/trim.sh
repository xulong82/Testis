#!/bin/sh
#
# RNA-seq data analysis pipeline 
# Author: XuLong Wang (xulong.wang@jax.org)

# Read trimming with Trimmomatic

echo $0

module load java/1.7.0

dir1="/data/xwang/Testis/FASTQ"
dir2="/data/xwang/Testis/FASTQ1"

files=`find $dir1 -name '*_R1_ALL.fastq'`

for name1 in $files; do
  name2=`basename $name1`
  name3=${name2/_R[12]_ALL.fastq/}
  echo $name3
  java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 20 -phred33 \
            "$dir1"/"$name3"_R1_ALL.fastq \
            "$dir1"/"$name3"_R2_ALL.fastq \
            "$dir2"/"$name3"_R1.fastq \
            "$dir2"/"$name3"_unpaired_R1.fastq \
            "$dir2"/"$name3"_R2.fastq \
            "$dir2"/"$name3"_unpaired_R2.fastq \
            LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:60
done

