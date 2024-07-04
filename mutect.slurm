#!/bin/sh
### sbatch pipeline.sh Path-For-Normal(upto / of last) PACA-104N path-For-Tumor PACA-104T   #####

#SBATCH --time=24:99:00    # Run time in hh:mm:ss
#SBATCH --mem=50g      # Memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15
#SBATCH --job-name=mutect50-PACA
#SBATCH --error=mutect-PACA-50.%J.err
#SBATCH --output=mutect-PACA-50.%J.out

### to run program  give comand sbatch pipeline.sh PACA-1130-N      #####
module load samtools/0.1 bwa/0.7
#bwa index human_g1k_v37.fasta

## Input file
N1='_end1.fq'
N2='_end2.fq'

###### bwa alignment ########
bwa mem -t $SLURM_NTASKS_PER_NODE -M human_g1k_v37.fasta $1$2$N1 $1$2$N2 > Normal.sam
bwa mem -t $SLURM_NTASKS_PER_NODE -M human_g1k_v37.fasta $3$4$N2 $3$4$N2 > Tumor.sam

####### Piccard based mark duplicate and sort on base of co-ordinare Normal ##########
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1.114/picard-tools-1.114/AddOrReplaceReadGroups.jar I=Normal.sam O=Normal.bam SO=coordinate ID=$2 PL=ILLUMINA SM=Normal PU=Normal LB=$2
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1.114/picard-tools-1.114/MarkDuplicates.jar INPUT=Normal.bam OUTPUT=Normal.marked.bam METRICS_FILE=Normal.matrix CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT


####### Piccard based mark duplicate and sort on base of co-ordinare Tumor ##########
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1.114/picard-tools-1.114/AddOrReplaceReadGroups.jar I=Tumor.sam O=Tumor.bam SO=coordinate ID=$4 PL=ILLUMINA SM=Tumor PU=Tumor LB=$4
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1.114/picard-tools-1.114/MarkDuplicates.jar INPUT=Tumor.bam OUTPUT=Tumor.marked.bam METRICS_FILE=Tumor.matrix CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT

######### GATK based local realignment Normal #######################################
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R human_g1k_v37.fasta -o Normal.bam.list -nt $SLURM_NTASKS_PER_NODE -I Normal.marked.bam
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -I Normal.marked.bam -R human_g1k_v37.fasta -targetIntervals Normal.bam.list -o Normal.marked.realigned.bam
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1.114/picard-tools-1.114/FixMateInformation.jar INPUT=Normal.marked.realigned.bam OUTPUT=Normal_bam.marked.realigned.fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T BaseRecalibrator  -I Normal_bam.marked.realigned.fixed.bam  -R human_g1k_v37.fasta -knownSites dbsnp_138.b37.vcf -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -nct $SLURM_NTASKS_PER_NODE -o Normal.recal_data.grp
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T PrintReads -I Normal_bam.marked.realigned.fixed.bam -R human_g1k_v37.fasta -BQSR Normal.recal_data.grp -nct $SLURM_NTASKS_PER_NODE -o Normal_bam.marked.realigned.fixed.final.bam

######### GATK based local realignment Tumor ########
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R human_g1k_v37.fasta -o Tumor.bam.list -nt $SLURM_NTASKS_PER_NODE -I Tumor.marked.bam
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -I Tumor.marked.bam -R human_g1k_v37.fasta -targetIntervals Tumor.bam.list -o Tumor.marked.realigned.bam
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1.114/picard-tools-1.114/FixMateInformation.jar INPUT=Tumor.marked.realigned.bam OUTPUT=Tumor_bam.marked.realigned.fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T BaseRecalibrator  -I Tumor_bam.marked.realigned.fixed.bam  -R human_g1k_v37.fasta -knownSites dbsnp_138.b37.vcf -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -nct $SLURM_NTASKS_PER_NODE -o Tumor.recal_data.grp
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T PrintReads -I Tumor_bam.marked.realigned.fixed.bam -R human_g1k_v37.fasta -BQSR Tumor.recal_data.grp -nct $SLURM_NTASKS_PER_NODE -o Tumor_bam.marked.realigned.fixed.final.bam

#Samtools index file.bam file.bai
#####################################################################################
java -Xmx20g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/muTect-1.1.4/muTect-1.1.4.jar -T MuTect -R human_g1k_v37.fasta --dbsnp  --cosmic  --input_file:normal Normal_bam.marked.realigned.fixed.final.bam --input_file:tumor Tumor_bam.marked.realigned.fixed.final.bam --out PACA-50.out --coverage_file PACA-50_coverage_file --vcf PACA-50.vcf --log_to_file PACA-50.log --performanceLog PACA-50_performance_log -nt  $SLURM_NTASKS_PER_NODE

########### Remove sam file to save space ###########################################
rm Normal.sam Tumor.sam

