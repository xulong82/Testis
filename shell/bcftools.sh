#!/bin/bash

# module load samtools
# module load bcftools
# module load tabix 

cd ${HOME}/Dropbox/GitHub/Testis/shiny/gviz

# samtools faidx seqs.fa  # index the fa file 

for name1 in M1IN_GES14_04904_TGACCA M2IN_GES14_04905_ACAGTG M3IN_GES14_04906_GCCAAT; do
  name2=${name1}_LaneALL.bam
  echo $name2

# samtools mpileup -uf seqs.fa bam/$name2 | bcftools call -mv -Ov > vcf/$name1.vcf
# bcftools filter -i'(%QUAL>20 && DP>50)' vcf/$name1.vcf | bgzip > vcf/${name1}_filter.vcf.gz 
# tabix vcf/${name1}_filter.vcf.gz

#  ~/Applications/bcftools consensus -f seqs.fa vcf/${name1}_filter.vcf.gz > consensus/$name1.fa
done

bcftools merge M1IN_GES14_04904_TGACCA_filter.vcf.gz \
	M2IN_GES14_04905_ACAGTG_filter.vcf.gz \
	M3IN_GES14_04906_GCCAAT_filter.vcf.gz | bgzip > MIN_merged.vcf.gz

tabix MIN_merged.vcf.gz

~/Applications/bcftools consensus -f seqs.fa vcf/MIN_merged.vcf.gz > MIN.fa

