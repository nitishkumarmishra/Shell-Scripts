#!/bin/sh
#SBATCH --time=15:00:00          # Run time in hh:mm:ss
#SBATCH --mem=40g       # Maximum memory required per CPU (in megabytes)
#SBATCH --nodes=1
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=fastqc
#SBATCH --error=/work/unmc_gudalab/nitish123/Radha/job.%J.err
#SBATCH --output=/work/unmc_gudalab/nitish123/Radha/job.%J.out

module load compiler/gcc
module load fastqc/0.11

for i in */*.fastq
do
        fastqc -o fastqc_results/ $i
done

