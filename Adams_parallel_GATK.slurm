#!/bin/sh

### to run program  give comand sbatch pipeline.sh PACA-1130-N      #####

#SBATCH --time=75:99:00    # Run time in hh:mm:ss
#SBATCH --mem=50g      # Memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=50-T-PACA
#SBATCH --error=PACA-50-T.%J.err
#SBATCH --output=PACA-50-T.%J.out
module load samtools/0.1 bwa/0.7
## Input file
N1='_end1.fq'
N2='_end2.fq'

###### bwa alignment ########
bwa mem -t $SLURM_NTASKS_PER_NODE -M /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa $1$N1 $1$N2 > L001.sam
#bwa mem -t $SLURM_NTASKS_PER_NODE -M /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa $T1fq $T2fq > L002.sam

####### Piccard based mark duplicate and sort on base of co-ordinare ##########
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1.114/picard-tools-1.114/AddOrReplaceReadGroups.jar I=L001.sam O=L001.bam SO=coordinate ID=$1 PL=ILLUMINA SM=L001 PU=L001 LB=$1
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1.114/picard-tools-1.114/MarkDuplicates.jar INPUT=L001.bam OUTPUT=L001.marked.bam METRICS_FILE=L001.matrix CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT

######### GATK based local realignment ########
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -o L001.bam.list -nt $SLURM_NTASKS_PER_NODE -I L001.marked.bam
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -I L001.marked.bam -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -targetIntervals L001.bam.list -o L001.marked.realigned.bam
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1.114/picard-tools-1.114/FixMateInformation.jar INPUT=L001.marked.realigned.bam OUTPUT=L001_bam.marked.realigned.fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T BaseRecalibrator  -I L001_bam.marked.realigned.fixed.bam  -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -knownSites /shared/unmc1/unmc_gudalab/nitish123/dbSNP/dbSNP-138/dbsnp_138.hg19_UCSC.vcf -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -nct $SLURM_NTASKS_PER_NODE -o L001.recal_data.grp
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T PrintReads -I L001_bam.marked.realigned.fixed.bam -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -BQSR L001.recal_data.grp -nct $SLURM_NTASKS_PER_NODE -o L001_bam.marked.realigned.fixed.final.bam
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T HaplotypeCaller -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -I L001_bam.marked.realigned.fixed.final.bam -D /shared/unmc1/unmc_gudalab/nitish123/dbSNP/dbSNP-138/dbsnp_138.hg19_UCSC.vcf -o L001_snps.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -A Coverage -A AlleleBalance -A HaplotypeScore -A InbreedingCoeff

##### Call SNP ########
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T SelectVariants -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -V L001_snps.vcf -selectType SNP -o L001_raw_snps.vcf
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T VariantFiltration -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -V L001_raw_snps.vcf --filterExpression "DP < 10 " --filterName "LowCoverage" --filterExpression "QD < 2.0" --filterName "LowQD"  --filterExpression "FS > 60.0" --filterName "PoorFS" --filterExpression "MQ < 40.0" --filterName "LowMQ"  --filterExpression "MQRankSum < -12.5" --filterName "LowMQSum"  --filterExpression "ReadPosRankSum < -8.0" --filterName "LowRank" -o L001_filtered_snps.vcf

###### Call Indel #######
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T SelectVariants -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -V L001_snps.vcf -selectType INDEL -o L001_raw_indels.vcf

java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T VariantFiltration -R /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/hg19-Pancreatic/hg19_UCSC.fa -V L001_raw_indels.vcf --filterExpression "DP < 10" --filterName "LowCoverage" --filterExpression "QD < 2.0" --filterName "LowQD"  --filterExpression "FS > 200.0" --filterName "PoorFS" --filterExpression "ReadPosRankSum < -20.0" --filterName "LowRank" -o L001_filtered_indels.vcf
######### Annovar annotation::  convert2annovar.pl ##########
/shared/unmc1/unmc_gudalab/nitish123/Tools/annovar/convert2annovar.pl --format vcf4old --includeinfo L001_filtered_snps.vcf > L001_filtered_snps.annovar1
/shared/unmc1/unmc_gudalab/nitish123/Tools/annovar/convert2annovar.pl --format vcf4old --includeinfo L001_filtered_indels.vcf > L001_filtered_indels.annovar1


grep "PASS"  L001_filtered_snps.annovar1 > L001_filtered_snps.annovar
grep "PASS"  L001_filtered_indels.annovar1 > L001_filtered_indels.annovar
#rm L001_filtered_snps.annovar1
#rm L001_filtered_indels.annovar1

########## table_annovar.pl ###########
/shared/unmc1/unmc_gudalab/nitish123/Tools/annovar/table_annovar.pl L001_filtered_snps.annovar  /shared/unmc1/unmc_gudalab/nitish123/Tools/annovar/humandb/  -buildver hg19 -out L001_final_snps -remove -protocol refGene,knownGene,ensGene,cytoBand,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp138,ljb23_all,cosmic68 -operation g,g,g,r,r,f,f,f,f,f -nastring . -csvout

/shared/unmc1/unmc_gudalab/nitish123/Tools/annovar/table_annovar.pl L001_filtered_indels.annovar /shared/unmc1/unmc_gudalab/nitish123/Tools/annovar/humandb/  -buildver hg19 -out  L001_final_indels -remove -protocol refGene,knownGene,ensGene,cytoBand,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp138,ljb23_all,cosmic68 -operation g,g,g,r,r,f,f,f,f,f -nastring . -csvout
