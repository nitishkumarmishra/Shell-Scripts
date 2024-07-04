#!/bin/sh
#SBATCH --time=60:00:00          # Run time in hh:mm:ss
#SBATCH --mem=50g       # Maximum memory required per CPU (in megabytes)
#SBATCH --nodes=1
#SBATCH --partition=batch
#SBATCH --job-name=HTSeq_2
#SBATCH --output=slurmout/arrayJob_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/arrayJob_%A_%a.err # File to which STDERR will be written
module load compiler/gcc
module load star/2.7
module load htseq/py37/0.11


start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=samples.txt
REF=/work/unmc_gudalab/nitish123/STAR_hg38/STAR_hg38_V33/
GTF=/work/unmc_gudalab/nitish123/STAR_hg38/gencode.v33.annotation.gtf

outpath='02-STAR_alignment'
[[ -d ${outpath} ]] || mkdir ${outpath}

do
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}
echo "SAMPLE: ${sample}"
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

call= "STAR --runThreadN 8 \
--readFilesIn Dataset2/${sample}1.fastq.gz Dataset2/${sample}.fastq.gz \ 
--alignIntronMax 1000000 \
--alignIntronMin 20 \
--alignMatesGapMax 1000000 \
--alignSJDBoverhangMin 1 \
--alignSJoverhangMin 8 \
--alignSoftClipAtReferenceEnds Yes \
--chimJunctionOverhangMin 15 \
--chimMainSegmentMultNmax 1 \
--chimOutType WithinBAM \
--chimSegmentMin 15 \
--genomeDir $REF \
--sjdbGTFfile $GTF \
--genomeLoad NoSharedMemory \
--limitSjdbInsertNsj 1200000 \
--outFilterIntronMotifs None \
--outFilterMatchNminOverLread 0.33 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--outFilterMultimapNmax 20 \
--outFilterScoreMinOverLread 0.33 \
--outFilterType BySJout \
--outSAMattributes NH HI AS nM NM ch \
--outSAMstrandField intronMotif \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--outFileNamePrefix ${outpath}/${sample}/${sample}_ \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat \
--twopassMode Basic" 

echo $call
eval $call
done
end=`date +%s`
runtime=$((end-start))
echo $runtime
