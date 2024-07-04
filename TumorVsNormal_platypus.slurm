#!/bin/sh
#SBATCH --time=24:99:00    # Run time in hh:mm:ss
#SBATCH --mem=50g      # Memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=PACA1130
#SBATCH --error=PACA1130.%J.err
#SBATCH --output=PACA1130.%J.out

######## Load Module ######################
module load samtools/0.1 bwa/0.7

######## Samtool Index and SamToBam #######

bwa mem -t $SLURM_NTASKS_PER_NODE -M /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/reads-src/SRP014665/SRS353180/SRX173375/SRR532052_processed/PACA-1130-N_end1.fq /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/reads-src/SRP014665/SRS353180/SRX173375/SRR532052_processed/PACA-1130-N_end2.fq > Normal.sam
bwa mem -t $SLURM_NTASKS_PER_NODE -M /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/reads-src/SRP014665/SRS353185/SRX173403/SRR532054_processed/PACA-1130-T_end1.fq /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/reads-src/SRP014665/SRS353185/SRX173403/SRR532054_processed/PACA-1130-T_end2.fq > Tumor.sam

java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1.114/picard-tools-1.114/SortSam.jar SO=coordinate INPUT=Normal.sam OUTPUT=Normal.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1.114/picard-tools-1.114/SortSam.jar SO=coordinate INPUT=Tumor.sam OUTPUT=Tumor.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

#############################################
rm  Normal.sam Tumor.sam
samtools index Normal.bam
samtools index Tumor.bam

######## run Platypus #######################
python /shared/unmc1/unmc_gudalab/nitish123/Tools/Platypus_0.7.4/Platypus.py callVariants --bamFiles=Normal.bam,Tumor.bam --refFile=/shared/unmc1/unmc_gudalab/nitish123/DbGAP/Reference-Platypus/hg19_UCSC.fa --minMapQual=30 --minBaseQual=30 --minGoodQualBases=25 --assemble=1 --nCPU=$SLURM_NTASKS_PER_NODE --output=L001.vcf

############ Short VCF file #################
grep "^#" L001.vcf >>L002.vcf
grep -w "chr1" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr2" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr3" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr4" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr5" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr6" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr7" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrX" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr8" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr9" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr10" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr11" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr12" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr13" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr14" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr15" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr16" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr17" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr18" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr20" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrY" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr19" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr22" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr21" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr6_ssto_hap7" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr6_mcf_hap5" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr6_cox_hap2" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr6_mann_hap4" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr6_apd_hap1" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr6_qbl_hap6" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr6_dbb_hap3" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr17_ctg5_hap1" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr4_ctg9_hap1" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr1_gl000192_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000225" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr4_gl000194_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr4_gl000193_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000222" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000212" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr7_gl000195_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000224" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000219" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr17_gl000205_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000215" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000216" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000217" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr9_gl000199_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000211" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000213" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000220" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000218" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr19_gl000209_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000221" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000214" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000228" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000227" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr1_gl000191_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr19_gl000208_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr9_gl000198_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr17_gl000204_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000233" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000237" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000230" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000243" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000241" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000240" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000232" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000234" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr11_gl000202_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000238" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr17_gl000203_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr8_gl000197_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000235" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000239" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chr21_gl000210_random" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000231" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000229" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrM" L001.vcf |sort -nk2 >>L002.vcf
grep -w "chrUn_gl000226" L001.vcf |sort -nk2 >>L002.vcf

################## Call SNP #######################
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T SelectVariants -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -V L002.vcf -selectType SNP -o L001_raw_snps.vcf

java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T VariantFiltration -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -V L001_raw_snps.vcf --filterExpression "QD < 7.0" --filterName "LowQD" --filterExpression "MQ < 40.0" --filterName "LowMQ" --filterExpression "MMLQ <25"  --filterName "LowMMLQ" --filterExpression "HapScore >5" --filterName "Low_HapScore" --filterExpression "SbPval > 0.5" --filterName  "Low_SbPval"  --filterExpression "BRF >=0.50" --filterName  "lowBRF" -o L001_filtered_snps.vcf

#java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T VariantFiltration -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -V L001_raw_snps.vcf --filterExpression "QD < 2.0" --filterName "QD" --filterExpression "MQ < 40.0" --filterName "MQ" --filterExpression "HapScore >4" --filterName "HapScore" --filterExpression "SbPval > 0.75" --filterName  "alleleBias"  --filterExpression "TR < 4" --filterName "HARD_TO_VALIDATE" -o L001_filtered_snps.vcf

################## Call Indel ######################
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T SelectVariants -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -V L002.vcf -selectType INDEL -o L001_raw_indels.vcf

java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T VariantFiltration -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -V L001_raw_indels.vcf --filterExpression "QD < 5.0" --filterName "LowQD" --filterExpression "MQ < 30.0" --filterName "LowMQ" --filterExpression "MMLQ <25"  --filterName "LowMMLQ"  --filterExpression "HapScore >4" --filterName "Low_HapScore" --filterExpression "SbPval > 0.5" --filterName  "Low_SbPval"  --filterExpression "BRF >=0.50" --filterName  "lowBRF" -o L001_filtered_indels.vcf

#java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T VariantFiltration -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -V L001_raw_indels.vcf --filterExpression "QD < 2.0" --filterName "QD" --filterExpression "MQ < 30.0" --filterName "MQ" --filterExpression "HapScore >4" --filterName "HapScore" --filterExpression "SbPval > 0.85" --filterName  "alleleBias" --filterExpression "TR < 2" --filterName "HARD_TO_VALIDATE" -o L001_filtered_indels.vcf

########## Annovar:: convert2annovar.pl ###########
/shared/unmc1/unmc_gudalab/nitish123/Tools/annovar/convert2annovar.pl --format vcf4old --includeinfo L001_filtered_snps.vcf > L001_filtered_snps.annovar1
/shared/unmc1/unmc_gudalab/nitish123/Tools/annovar/convert2annovar.pl --format vcf4old --includeinfo L001_filtered_indels.vcf > L001_filtered_indels.annovar1

################ Select only PASS #################
grep "PASS"  L001_filtered_snps.annovar1 > L001_filtered_snps.annovar
grep "PASS"  L001_filtered_indels.annovar1 > L001_filtered_indels.annovar
rm L001_filtered_snps.annovar1
rm L001_filtered_indels.annovar1

################ Annovar.pl ################
/shared/unmc1/unmc_gudalab/nitish123/Tools/annovar/table_annovar.pl  L001_filtered_snps.annovar  /shared/unmc1/unmc_gudalab/nitish123/Tools/annovar/humandb/ -buildver hg19 -out L001_final_snps -remove -protocol refGene,knownGene,ensGene,snp138,cosmic68,cosmic68wgs,cytoBand,genomicSuperDups,1000g2012apr_all,ljb23_all -operation g,g,g,f,f,f,r,r,f,f -nastring . -csvout

/shared/unmc1/unmc_gudalab/nitish123/Tools/annovar/table_annovar.pl  L001_filtered_indels.annovar /shared/unmc1/unmc_gudalab/nitish123/Tools/annovar/humandb/ -buildver hg19 -out  L001_final_indels -remove -protocol refGene,knownGene,ensGene,snp138,cosmic68,cosmic68wgs,cytoBand,genomicSuperDups,1000g2012apr_all,ljb23_all -operation g,g,g,f,f,f,r,r,f,f -nastring . -csvout
