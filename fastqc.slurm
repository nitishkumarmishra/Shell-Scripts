#!/bin/sh
#SBATCH --time=5:99:00    # Run time in hh:mm:ss    
#SBATCH --mem=20g      # Memory
#SBATCH --mem-per-cpu=10240       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=101-T_SOLiD.bam_qsort_end1
#SBATCH --error=PACA-101-T_SOLiD.bam_qsort_end1%J.err
#SBATCH --output=PACA-101-T_SOLiD.bam_qsort_end1%J.out
module load fastqc/0.10
fastqc reads-src/SRP014665/SRS352735/SRX172901/SRR531341_processed/PACA-101-T_SOLiD.bam_qsort_end1.fq -o /shared/unmc1/unmc_gudalab/nitish123/DbGAP/Pancreatic/FastQC-Report-SOLiD
