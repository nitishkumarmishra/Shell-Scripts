#!/bin/sh
#*********** This program only run on Tusker #############
#SBATCH --time=1:00:00    # Run time in hh:mm:ss
#SBATCH --mem=2g      # Memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15
#SBATCH --job-name=PACA86
#SBATCH --error=PACA86.%J.err
#SBATCH --output=PACA86.%J.out

########## run vcf-filter ###########
#perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' PACA86/L001_filtered_snps_PASS.vcf > fpfilter/PACA86_PASS.var

#bam-readcount -q1 -b15 -w1 -l fpfilter/PACA86_PASS.var -f /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Reference-Platypus/hg19_UCSC.fa /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Reference-Platypus/TumorVsNormal/PACA86/Tumor.bam > fpfilter/PACA86_PASS.readcount

#perl /home/unmc_gudalab/nitish123/variant-filter/variant-filter-master/fpfilter.pl --var-file PACA86/L001_filtered_snps_PASS.vcf --readcount-file fpfilter/PACA86_PASS.readcount --output-file fpfilter/PACA86_PASS.fpfilter
grep "^#" PACA86/L001_filtered_snps_PASS.vcf >header
grep -v "^#" PACA86/L001_filtered_snps_PASS.vcf >fpfilter/PACA86_vcf
grep -v "^#" fpfilter/PACA86_PASS.fpfilter >fpfilter/PACA86_fpfilter
paste -d "	" fpfilter/PACA86_vcf fpfilter/PACA86_fpfilter |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$19"\t"$8"\t"$9"\t"$10"\t"$11}' >PACA86_fpfilter_vcf
grep "PASS" PACA86_fpfilter_vcf >PACA86_fpfilter_vcf1
cat header PACA86_fpfilter_vcf1 >PACA86_fpfilter.vcf
rm header fpfilter/PACA86_vcf fpfilter/PACA86_fpfilter PACA86_fpfilter_vcf PACA86_fpfilter_vcf1
