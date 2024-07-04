#!/bin/sh
#SBATCH --time=75:99:00    # Run time in hh:mm:ss
#SBATCH --mem=50g      # Memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=$1
#SBATCH --error=$1.%J.err
#SBATCH --output=$1.%J.out
module load samtools/0.1 bwa/0.7
## Input file
R1='_end1.fq'
R2='_end2.fq'
####### bwa alignment ########
#R1="raw_data/SRR098401_1.fastq"
#R2="raw_data/SRR098401_2.fastq"
INDEX="genome/hg19"
FA="genome/hg19.fa"
samtools faidx $FA
/home/unmc_gudalab/nitish123/software/mrsFAST/mrsfast/mrsfast --index $FA
/home/unmc_gudalab/nitish123/software/mrsFAST/mrsfast/mrsfast --search $FA --pe --seq1 $$1R1 --seq2 $$1R2 --threads $SLURM_NTASKS_PER_NODE -o $1.mrsfast --best --disable-nohits
# The following three commands helped to remove errors identified by samtools / PicardTools, but may not be necessary for your work
grep -vP "MD:Z:$" $1.mrsfast.sam > $1.mrsfast_fixed.sam
samtools view -bS $1.mrsfast_fixed.sam -o $1.mrsfast.bam
java -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1.114/picard-tools-1.114/SortSam.jar INPUT=$1.mrsfast.bam OUTPUT=$1.mrsfast.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT

