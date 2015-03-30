#!/bin/sh
 
module load bowtie/1.0.0
module load rsem

dir1="/data/fysun/R8-d14_polysome_profiling-RNA-seq"
dir2="/data/fysun/RSEM"
ref="/data/fysun/REF/GRCm38.transcripts"

# files=`find $dir1 -name '*LaneALL_R1_ALL.fastq'`

# filename3="M3PLM_GES14_04918_GTGAAA_LaneALL"
filename3="W2NONP_GES14_04908_ACTTGA_LaneALL"

# for filename1 in $files; do
# filename2=`basename $filename1`
# filename3=${filename2/_R1_ALL.fastq/}
  echo $filename3
  rsem-calculate-expression -p 20 \
 			    --bowtie-phred33-quals \
 			    --forward-prob 0.5 \
 			    --paired-end \
                            "$dir1"/"$filename3"_R1_ALL.fastq \
                            "$dir1"/"$filename3"_R2_ALL.fastq \
 			    "$ref" \
 			    "$dir2"/"$filename3"
# done

#  rsem-calculate-expression -p 10 \
# 			    --bowtie-phred33-quals \
# 			    --forward-prob 0.5 \
# 			    --paired-end \
#                            "$dir1"/W3PLM_GES14_04912_CTTGTA_LaneALL_R1_ALL.fastq \
#			    "$dir1"/W3PLM_GES14_04912_CTTGTA_LaneALL_R2_ALL.fastq \
# 			    "$ref" \
# 			    "$dir2"/"$filename3"
# Options
# Make the genome bam (Caution: This take 10 hours per sample)
#        		    --output-genome-bam 
# Generating a Wiggle file
#rsem-bam2wig sorted_bam_input wig_output.wig wiggle_name

# Prepare reference sequences with RSEM
#rsem-prepare-reference --gtf ./GTF/mm10knownGene.gtf \
#		       --transcript-to-gene-map ./GTF/mm10knownIsoforms.txt \
#                       ./Genome/mm10.fa \
#                       ./RSEM/mm10rsem
#
#rsem-prepare-reference --gtf ./GTF/Mus_musculus.GRCm38.73.gtf \
#                       ./Genome/GRCm38.73 \
#                       ./RSEM/NCBIM37.59
#
#rsem-prepare-reference --gtf ./GTF/Mus_musculus.NCBIM37.59.gtf \
#                       ./Genome/mm9/regular \
#                       ./RSEM/NCBIM37
#
