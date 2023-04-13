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

mkdir -p GENOME_data/rsem/tRNA


#rsem-prepare-reference --gtf /home/nmishra/REF/Mus_musculus/REF_GRCm39_vM27/gencode.vM27.annotation.gtf \
#--star \
#--star-path ~/.conda/envs/cutadaptenv/bin/STAR \
#-p 8 \
#/home/nmishra/REF/Mus_musculus/REF_GRCm39_vM27/RSEM/GRCm39.primary_assembly.genome.fa
#GENOME_data/rsem/vM27


#rsem-prepare-reference --gtf /home/nmishra/REF/Mus_musculus/REF_GRCm39_vM27/gencode.vM27.annotation.gtf \
#                       --star \
#                       /home/nmishra/REF/Mus_musculus/REF_GRCm39_vM27/RSEM/GRCm39.primary_assembly.genome.fa \
#                       GENOME_data/rsem/vM27


rsem-prepare-reference --gtf $GTF \
		       -p 8 \
		       --star \
		       /home/nmishra/REF/Mus_musculus/REF_GRCm39_vM27/RSEM/GRCm39.primary_assembly.genome.fa \
		       GENOME_data/rsem/tRNA/vM27_tRNA
