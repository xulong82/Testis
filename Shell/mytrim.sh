#!/bin/sh
#
# RNA-seq data analysis pipeline 
# Author: XuLong Wang (xulong.wang@jax.org)

# Read trimming with Trimmomatic

echo $0
begin=`date +%h`

module load java/1.7.0

dir_sou="/hpcdata/xwang/AD/howell_2014/retina"
dir_des="/hpcdata/xwang/AD/howell_2014/retina_trim"
#
files=`find $dir_sou -name '*.fastq'`

for filename1 in $files; do
  filename2=`basename $filename1`
  filename3=${filename2/_R[12].fastq/}
  echo $filename3
  java -jar /home/xwang/bin/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 20 -phred33 \
            "$dir_sou"/"$filename3"_R1.fastq \
            "$dir_sou"/"$filename3"_R2.fastq \
            "$dir_des"/"$filename3"_R1.fastq \
            "$dir_des"/"$filename3"_unpaired_R1.fastq \
            "$dir_des"/"$filename3"_R2.fastq \
            "$dir_des"/"$filename3"_unpaired_R2.fastq \
	    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:60
done

end=`date +%h`
echo $((end-begin))

