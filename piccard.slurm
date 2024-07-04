#!/bin/sh
#SBATCH --time=24:99:00    # Run time in hh:mm:ss    
#SBATCH --mem=60g      # Memory
#SBATCH --mem-per-cpu=20480       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=piccard-PACA-4520-T_SOLiD
#SBATCH --error=/shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/PACA-4520-T_SOLiD-piccard%J.err
#SBATCH --output=/shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/PACA-4520-T_SOLiD-piccard%J.out
#module load samtools/0.1
java -Xmx40g -jar /shared/unmc1/unmc_gudalab/nitish123/Tools/picard-tools-1/picard-tools-1.111/SamToFastq.jar I=reads-src/SRP014665/SRS352730/SRX172895/SRR531328_processed/PACA-4520-T_SOLiD.bam_qsort.bam F=reads-src/SRP014665/SRS352730/SRX172895/SRR531328_processed/PACA-4520-T_SOLiD_end1.fq F2=reads-src/SRP014665/SRS352730/SRX172895/SRR531328_processed/PACA-4520-T_SOLiD_end2.fq
