#!/bin/bash
#BSUB -P STAR_COUNT
#BSUB -oo STAR_COUNT.out -eo STAR_COUNT.err
#BSUB -n 15
#BSUB -M 30G
#BSUB -J STAR_COUNT

####################################################
############ GTF and STAR genome Index #############
OUTPUT_DIR=/research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/Trial
GENOME=/home/nmishra/REF/Mus_musculus/REF_GRCm39_vM27/STAR-index/2.7.9a
#GTF=/research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/Trial/mm39-tRNAs.gtf
GTF=/research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/Trial/mm39-tRNAs_exclude_Pseudogene_High_Confidense.gtf
NUMBER_THREAD=10

####################################################
############### Parsing FASTQ file #################
PREFIX=$1
FQ1=`ls|grep $PREFIX | grep "_R1_"| grep "fastq.gz"`
FQ2=`ls|grep $PREFIX | grep "_R2_"| grep "fastq.gz"`

####################################################
######### Trimgalore adapter trimming ##############
### Conda environ for latest cutadapt and trmgalore.
source /hpcf/apps/conda3/install/5.1.0/bin/activate cutadaptenv
mkdir -p ./FASTQC/TRIM #Output file for FASTQ

#/home/nmishra/TrimGalore-0.6.6/trim_galore \
trim_galore \
-j $NUMBER_THREAD \
--clip_R1 4 \
--clip_R2 4 \
--three_prime_clip_R1 4 \
--three_prime_clip_R2 4 \
-a TGGAATTCTCGGGTGCCAAGG \
-a2 TGGAATTCTCGGGTGCCAAGG \
--basename $PREFIX \
--length 23 \
--paired \
--fastqc_args "--outdir ./FASTQC/TRIM" $FQ1 $FQ2

TRIM1=${PREFIX}_val_1.fq.gz
TRIM2=${PREFIX}_val_2.fq.gz
####################################################
############## STAR GeneCounts #####################
STAR \
--genomeDir $GENOME \
--readFilesCommand zcat \
--runThreadN $NUMBER_THREAD \
--sjdbGTFfile $GTF \
--alignEndsType EndToEnd \
--outFilterMismatchNmax 50 \
--outSAMattributes NH HI AS nM NM ch \
--outFilterMultimapScoreRange 10 \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMprimaryFlag AllBestScore \
--outFileNamePrefix ${PREFIX}_ \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 999 \
--outFilterMismatchNoverReadLmax 0.33 \
--outSAMunmapped Within \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--outFilterMatchNmin 5 \
--alignSJDBoverhangMin 1000 \
--outWigType wiggle \
--outWigStrand Stranded \
--outWigNorm RPM \
--twopassMode Basic \
--readFilesIn $TRIM1 $TRIM2

####################################################
############## RSEM GeneCounts #####################
# $1 == STAR alignments bam file
# $2 == name of STAR index
# $3 == sample name
mkdir -p ./RSEM_COUNTS/
RSEM_OUT=./RSEM_COUNTS

rsem-calculate-expression \
-p 8 \
--no-bam-output  \
--alignments  \
--paired-end  \
--strandedness reverse \
${PREFIX}_Aligned.toTranscriptome.out.bam \
/research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/GENOME_data/rsem/tRNA/vM27_tRNA \
${PREFIX}.RSEM

