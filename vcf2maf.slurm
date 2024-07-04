#!/bin/sh
#SBATCH --time=5:99:00    # Run time in hh:mm:ss
#SBATCH --mem=30g      # Memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=PACA104
#SBATCH --error=PACA104.%J.err
#SBATCH --output=PACA104.%J.out

###########################################
perl /shared/unmc1/unmc_gudalab/nitish123/Tools/vcf2maf-master/vcf2maf.pl --input-vcf PACA104_fpfilter.vcf --output-maf PACA104_fpfilter.maf --tumor-id PACA104-TUMOR --normal-id PACA-104-NORMAL


