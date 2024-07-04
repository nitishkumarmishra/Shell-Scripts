#!/bin/sh
#SBATCH --time=50:00:00          # Run time in hh:mm:ss
#SBATCH --mem=50g       # Maximum memory required per CPU (in megabytes)
#SBATCH --nodes=1
#SBATCH --partition=batch
#SBATCH --job-name=HTSeq1
#SBATCH --error=/work/unmc_gudalab/nitish123/Radha/RNAseq/Batch1/portal.us.medgenome.com/P2001043_03302020/job.%J.err
#SBATCH --output=/work/unmc_gudalab/nitish123/Radha/RNAseq/Batch1/portal.us.medgenome.com/P2001043_03302020/job.%J.out

module load compiler/gcc
module load star/2.7
module load htseq/py37/0.11

STAR \
--readFilesIn UR-15111_R1.fastq.gz,UR-15112_R1.fastq.gz,UR-15113_R1.fastq.gz,UR-15114_R1.fastq.gz,UR-15115_R1.fastq.gz UR-15111_R2.fastq.gz,UR-15112_R2.fastq.gz,UR-15113_R2.fastq.gz,UR-15114_R2.fastq.gz,UR-15115_R2.fastq.gz \
--alignIntronMax 1000000 \
--alignIntronMin 20 \
--alignMatesGapMax 1000000 \
--alignSJDBoverhangMin 1 \
--alignSJoverhangMin 8 \
--chimJunctionOverhangMin 15 \
--chimMainSegmentMultNmax 1 \
--chimOutType WithinBAM \
--chimSegmentMin 15 \
--genomeDir /work/unmc_gudalab/nitish123/STAR_hg38/STAR_hg38_V33/ \
--genomeLoad NoSharedMemory \
--limitSjdbInsertNsj 1200000 \
--outFilterIntronMotifs None \
--outFilterMatchNminOverLread 0.33 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--outFilterMultimapNmax 20 \
--outFilterScoreMinOverLread 0.33 \
--outFilterType BySJout \
--outSAMattributes NH HI AS nM NM ch \
--outSAMstrandField intronMotif \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat \
--twopassMode Basic
