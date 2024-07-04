#!/bin/sh
#SBATCH --time=24:99:00    # Run time in hh:mm:ss    
#SBATCH --mem=40g      # Memory
#SBATCH --mem-per-cpu=2048       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=phs000245-700064-xenograft_fastq
#SBATCH --error=/shared/unmc1/unmc_gudalab/nitish123/TripleBreast/reads-src/SRP003602/SRS115049/SRX097143/SRR343359_processed/phs000245-700064-xenograft_fastq.%J.err
#SBATCH --output=/shared/unmc1/unmc_gudalab/nitish123/TripleBreast/reads-src/SRP003602/SRS115049/SRX097143/SRR343359_processed/phs000245-700064-xenograft_fastq.%J.out
module load samtools/0.1
/shared/unmc1/unmc_gudalab/nitish123/Tools/bedtools2-2.19.1/bin/bedtools bamtofastq -i /shared/unmc1/unmc_gudalab/nitish123/TripleBreast/reads-src/SRP003602/SRS115049/SRX097143/SRR343359_processed/phs000245-700064-xenograft_from_primary_qsort.bam -fq /shared/unmc1/unmc_gudalab/nitish123/TripleBreast/reads-src/SRP003602/SRS115049/SRX097143/SRR343359_processed/phs000245-700064-xenograft_from_primary_qsort_end1.fq -fq2 /shared/unmc1/unmc_gudalab/nitish123/TripleBreast/reads-src/SRP003602/SRS115049/SRX097143/SRR343359_processed/phs000245-700064-xenograft_from_primary_qsort_end2.fq

