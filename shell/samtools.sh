#!/bin/bash
#PBS -l mem=128gb,nodes=1:ppn=1,walltime=10:00:00

module load bowtie
module load samtools
module load java

# dir1="/data/xwang/Testis/RSEM"
# dir2="/home/xwang/Dropbox/GitHub/Testis/shiny/gviz/bam"
# 
# files=`find $dir1 -name '*.transcript.sorted.bam'`
# 
# cd $dir1
# 
# for name1 in $files; do
#   name2=`basename $name1`
#   name3=${name2/.transcript.sorted.bam/}
#   echo $name3
# # samtools view -b -q 1 $name2 \
#   samtools view -b $name2 \
#         ENSMUST00000080449 \
# 	ENSMUST00000084214 \
# 	ENSMUST00000084215 \
# 	ENSMUST00000105830 \
# 	ENSMUST00000105831 > $dir2/$name3.bam
# done

# cd /data/xwang/Testis/bam
# 
# sams=`find ./ -name '*.sam'`
# bams=`find ./ -name '*.sorted.bam'`
# mm10="/data/xwang/REF/Mm10Genome"
# 
# for idx in $sams; do
#   idx2=`basename $idx`
#   file=${idx2/.sam/}
#   echo $file
#   echo "sam to bam, sorting, indexing"
#   samtools view -bSF 4 -o ${file}.bam ${file}.sam
#   samtools sort ${file}.bam ${file}.sorted
#   samtools index ${file}.sorted.bam
#   samtools view -b ${file}.sorted.bam "chr12:76403176-76406936" > ${file}.hspa2.bam
# done

  
cd /data/xwang/Testis/bam

mm10="/data/xwang/REF/Mm10Genome"
hspa2=`find ./ -name '*.hspa2.bam'`

for idx in $hspa2; do
  idx2=`basename $idx`
  file=${idx2/.bam/}
  echo $file
  echo "sam to bam, sorting, indexing"
  samtools sort ${file}.bam ${file}.sorted
  samtools index ${file}.sorted.bam
done

