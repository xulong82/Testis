#!/bin/sh
#
# RNA-seq data analysis pipeline
# Author: XuLong Wang (xulong.wang@jax.org)

module load bowtie/1.0.0
module load rsem

cd /data/fysun/REF
# Build the reference index
bowtie-build GRCm38.transcripts.fa GRCm38.transcripts
rsem-prepare-reference GRCm38.transcripts.fa GRCm38.transcripts
