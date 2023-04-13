#!/bin/bash
#BSUB -P RSEM_COUNT
#BSUB -oo RSEM_COUNT.out -eo RSEM_COUNT.err
#BSUB -n 10
#BSUB -M 30G
#BSUB -J RSEM_COUNT

####################################################
############ GTF and STAR genome Index #############

GENOME=/home/nmishra/REF/Mus_musculus/REF_GRCm39_vM27/STAR-index/2.7.9a
#GTF=/research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/Trial/mm39-tRNAs.gtf
GTF=/research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/Trial/mm39-tRNAs_exclude_Pseudogene_High_Confidense.gtf
NUMBER_THREAD=10

####################################################
######### Trimgalore adapter trimming ##############
### Conda environ for latest cutadapt and trmgalore.
source /hpcf/apps/conda3/install/5.1.0/bin/activate cutadaptenv

mkdir -p GENOME_data/rsem


# $1 == STAR alignments bam file
# $2 == name of STAR index
# $3 == sample name
mkdir -p ./RSEM_COUNTS/
RSEM_OUT=./RSEM_COUNTS

#rsem-calculate-expression \
#-p 8 \
#--no-bam-output  \
#--alignments  \
#--paired-end  \
#--strandedness reverse \
#2113859_JBQS011_Aligned.toTranscriptome.out.bam \
#/research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/GENOME_data/rsem/tRNA/vM27_tRNA \
#trial.RSEM

#rsem-calculate-expression \
#--paired-end \
#--star \
#--star-path ~/.conda/envs/cutadaptenv/bin/STAR \
#--paired-end \
#-p 8 \
#2113866_JBQS018_val_1.fq.gz \
#2113866_JBQS018_val_2.fq.gz \
#/research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/GENOME_data/rsem/tRNA/vM27_tRNA \
#paired_end_quals


rsem-calculate-expression --paired-end --estimate-rspd --append-names --star --star-gzipped-read-file --output-genome-bam 2113866_JBQS018_val_1.fq.gz 2113866_JBQS018_val_2.fq.gz /research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/GENOME_data/rsem/tRNA/vM27_tRNA rsem_output_CTRL1 > CTRL1.log
#rsem-calculate-expression --paired-end \
#                               --paired-end \
#                               -p 8 \
#                               2113863_JBQS015_Aligned.toTranscriptome.out.bam \
#                               GENOME_data/rsem/vM27


#rsem-calculate-expression --paired-end 2113866_JBQS018_val_1.fq.gz 2113866_JBQS018_val_2.fq.gz /research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/GENOME_data/rsem/tRNA/vM27_tRNA JBQS018_val.RSEM

#rsem-calculate-expression --star --star-path ~/.conda/envs/cutadaptenv/bin/STAR --star-bzipped-read-file -p 8 --strandedness reverse --paired-end 2113866_JBQS018_val_1.fq.gz 2113866_JBQS018_val_2.fq.gz /research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/GENOME_data/rsem/tRNA/vM27_tRNA JBQS018_val.RSEM
