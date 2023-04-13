#!/bin/bash
#BSUB -P STAR_COUNT
#BSUB -oo STAR_COUNT.out -eo STAR_COUNT.err
#BSUB -n 12
#BSUB -M 25G
#BSUB -J STAR_COUNT

####################################################
######### Load library from module list ############
module load star/2.7.1a
#module load python/2.7.8
#module load trimgalore
module load fastqc/0.11.8
module load pigz
module load subread/1.5.1
## Latest cutadapt and trmgalore, installed in the environment.
conda activate cutadaptenv

####################################################
############ GTF and STAR genome Index #############
GENOME=/research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/GeneSeq/tRNA/STAR/2.7/mm39-mature-tRNA
GTF=/research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/GeneSeq/tRNA/gencode.vM27.tRNAs.gtf
NUMBER_THREAD=6

####################################################
############### Parsing FASTQ file #################
PREFIX=$1
FQ1=`ls|grep $PREFIX | grep "_R1_"| grep "fastq.gz"`
FQ2=`ls|grep $PREFIX | grep "_R2_"| grep "fastq.gz"`

####################################################
######### Trimgalore adapter trimming ##############
### Conda environ for latest cutadapt and trmgalore.
conda activate cutadaptenv 
mkdir -p ./FASTQC/TRIM #Output file for FASTQ

/home/nmishra/TrimGalore-0.6.6/trim_galore \
--paired \
--cores 2 \
--length 17 \
--retain_unpaired \
-a TGGAATTCTCGGGTGCCAAGG \
--trim-n 4 \
--output_dir ./ --fastqc_args "--outdir ./FASTQC/TRIM" \
--basename $PREFIX /
$FQ1 $FQ2

TRIM1=${PREFIX}_val_1.fq.gz
TRIM2=${PREFIX}_val_2.fq.gz
conda deactivate cutadaptenv

####################################################
############## STAR GeneCounts #####################
STAR \
--genomeDir $GENOME \
--readFilesCommand zcat \
--runThreadN $NUMBER_THREAD \
--sjdbGTFfile $GTF \
--alignEndsType EndToEnd \
--outFilterMismatchNmax 1 \
--outFilterMultimapScoreRange 0 \
--quantMode TranscriptomeSAM GeneCounts \
--outFileNamePrefix {$PREFIX}_ \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 10 \
--outSAMunmapped Within \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--outFilterMatchNmin 16 \
--alignSJDBoverhangMin 1000 \
--alignIntronMax 1 \
--outWigType wiggle \
--outWigStrand Stranded \
--outWigNorm RPM \
--readFilesIn $TRIM1 $TRIM2

####################################################
################ Subraed featureCounts #############
samtools index ${$PREFIX}_Aligned.sortedByCoord.out.bam
featureCounts -T $NUMBER_THREAD -s 2 -a $GTF -p -d 19 --countReadPairs -M  -O --fracOverlap 0.3 -Q 10 -o ${PREFIX}_featureCounts.counts ${$PREFIX}_Aligned.sortedByCoord.out.bam

####################################################
############### HTSeq counts #######################
module load python/2.7.15-rhel7
module load htseq/0.6.1p1-deprecated

htseq-count \
-f bam \
-r pos \
-s no \
-a 10 \
-t exon \
-i gene_id \
-m intersection-nonempty \
-n $NUMBER_THREAD
${$PREFIX}_Aligned.sortedByCoord.out.bam \
$GTF > {$PREFIX}_HTSeq.counts
