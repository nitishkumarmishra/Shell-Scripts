#!/bin/sh
#*********** This program only run on Tusker #############
#SBATCH --time=2:00:00    # Run time in hh:mm:ss
#SBATCH --mem=10g      # Memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15
#SBATCH --job-name=PACA46
#SBATCH --error=PACA46.%J.err
#SBATCH --output=PACA46.%J.out

########## run vcf-filter ###########
perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' PACA46/L001_filtered_snps_PASS.vcf > fpfilter/PACA46_PASS.var

bam-readcount -q1 -b15 -w1 -l fpfilter/PACA46_PASS.var -f /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Reference-Platypus/hg19_UCSC.fa /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Reference-Platypus/TumorVsNormal/PACA46/Tumor.bam > fpfilter/PACA46_PASS.readcount

perl /home/unmc_gudalab/nitish123/variant-filter/variant-filter-master/fpfilter.pl --var-file PACA46/L001_filtered_snps_PASS.vcf --readcount-file fpfilter/PACA46_PASS.readcount --output-file fpfilter/PACA46_PASS.fpfilter
