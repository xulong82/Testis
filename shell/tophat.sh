#!/bin/sh
#PBS -l mem=64gb,nodes=1:ppn=20,walltime=10:00:00

file=${arg}

module load samtools/0.1.19
module load bowtie/0.12.9
module load bowtie2/2.2.0
module load tophat/2.0.7

cd /data/xwang/Testis/FASTQ1

gtf="/data/xwang/C3H/C3H.gtf"
index="/data/xwang/C3H/C3H"

tophat --num-threads 20 \
       --GTF $gtf \
       --library-type fr-unstranded \
       --no-coverage-search \
       --output-dir ../Tophat/"$file" \
       $index "$file"_R1.fastq "$file"_R2.fastq

# cd /hpcdata/xwang/DCR
# for myfile in VpIpp2; do
# cufflinks --num-threads 32 \
#           --output-dir ./cufflinks/"$myfile" \
#           ./tophat/"$myfile"/accepted_hits.bam
# done

#!/bin/sh
# RNA-seq data analysis pipeline with Tophat and Cufflinks
# Author: XuLong Wang (xulong.wang@jax.org)

# Step 2: Assemble transcripts for each sample with Cufflinks
# Usage: sh mycufflinks.sh > mylog

#echo $0;
#
cd /hpcdata/xwang

echo "begin: `date`";

#
#for myfile in VnInp1 VnInp2 VnIpp1 VnIpp2 VpIpp1
#do
#echo $myfile
#cufflinks --num-threads 32 \
#	  --output-dir ./cufflinks/"$myfile" \
#	  ./tophat/"$myfile"/accepted_hits.bam
#done
#
cd ./DCR/mm10/cufflinks
mygtf="/hpcdata/xwang/Mus_musculus_iGenome/UCSC/mm10/Annotation/Genes/genes.gtf"
mygenome="/hpcdata/xwang/Mus_musculus_iGenome/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
#cd ./DCR/cufflinks
#mygtf="/hpcdata/xwang/GTF/Mus_musculus.GRCm38.73.gtf"
#mygenome="/hpcdata/xwang/Genome/GRCm38.73/GRCm38.73.fa"

cuffmerge --num-threads 32 \
	  --ref-gtf $mygtf \
	  --ref-sequence $mygenome \
	  ./assemblies.txt

echo "end: `date`";
#!/bin/sh
# AUTHOR: XULONG WANG (XULONG.WANG@JAX.ORG)

# Step 4: Cuffdiff
# Usage: sh mycuffdiff.sh > mylog

cd /hpcdata/xwang/DCR/mm10

echo "begin: `date`";

cuffdiff --num-threads 32 \
	 --output-dir cuffdiff \
	 --multi-read-correct \
	 ./merged_asm/merged.gtf \
	 ./tophat/VnInp1/accepted_hits.bam,./tophat/VnIpp1/accepted_hits.bam,./tophat/VpIpp1/accepted_hits.bam \
	 ./tophat/VnInp2/accepted_hits.bam,./tophat/VnIpp2/accepted_hits.bam,./tophat/VpIpp2/accepted_hits.bam

#	 --labels WT,Het \
#	 --frag-bias-correct ./Genome/mm9.fa \
echo "end: `date`";

