#!/bin/bash
#module load perl
module load flash
module load fastqc
module load java
module load bowtie/1.2.2
#module load bowtie/1.1.2 
module load bedtools
#module load viennarna/2.4.14 

vienna_perl="/rgs01/project_space/cab/automapper/common/tchang1/Softwares/viennaRNA/ViennaRNA-2.4.17/lib/site_perl/5.10.1"
experimental_perl="/rgs01/project_space/cab/automapper/common/tchang1/Softwares/sRNAtools/program/perl/experimental-0.022/modules/lib/site_perl/5.10.1"
vienna_bin="/rgs01/project_space/cab/automapper/common/tchang1/Softwares/viennaRNA/ViennaRNA-2.4.17/bin"
vienna_path="/rgs01/project_space/cab/automapper/common/tchang1/Softwares/viennaRNA/ViennaRNA-2.4.17/lib"

export PATH=$PATH:$vienna_bin
export PERL5LIB=$PERL5LIB:$vienna_perl:$experimental_perl
export LD_LIBRARY_PATH=$vienna_path:$LD_LIBRARY_PATH
#export PERL5LIB=$PERL5LIB:$experimental_perl


for fq in fastq/blancgrp_212001_miRNA-seq/*_R1_001.fastq.gz; do
#for fq in fastq/blancgrp_212001_miRNA-seq/2113858_JBQS010_S87_L002_R1_001.fastq.gz; do
	read1=$fq
	read2=${fq%%_R1_001.fastq.gz}_R2_001.fastq.gz

	sample=$( basename $read1 )
	sample=${sample%%_R1_001.fastq.gz}

	#echo $read1
	#echo $read2
	
	### s1 ... skipped
	mkdir -p miRNA_merged
	#cmd="flash $read1 $read2 -o $sample -d miRNA_merged --max-overlap 100"
	
	### s2
	#trim_folder="miRNA_trimmed" ### automatic
	trim_folder="miRNA_trimmed_v15" ### specify adpater
	illumina_smRNA_truseq_a='TGGAATTCTCGGGTGCCAAAG'
	illumina_smRNA_adapter_v1='TCGTATGCCGTCTTCTGCTTG'
        illumina_smRNA_adapter_v15='ATCTCGTATGCCGTCTTCTGCTTG'
	mkdir -p $trim_folder/fastqc
	#cmd="trim_galore --paired $read1 $read2 --length 18 -e 0.1 -q 30 --stringency 3 -o $trim_folder --fastqc --fastqc_args \"--nogroup --outdir $trim_folder/fastq\"" 
	cmd="trim_galore --paired $read1 $read2 -a $illumina_smRNA_adapter_v15 --length 18 -e 0.1 -q 30 --stringency 3 -o $trim_folder --fastqc --fastqc_args \"--nogroup --outdir $trim_folder/fastq\""	

	#bsub -M 4000 -J mrg -e mrg.err -o mrg.out "$cmd"


	### s3
	read1="$trim_folder/${sample}_R1_001_val_1.fq.gz"
	read2="$trim_folder/${sample}_R2_001_val_2.fq.gz"

	if [ -s $read1 ] && [ -s $read2 ]; then
		cmd="/rgs01/project_space/cab/automapper/common/tchang1/Softwares/bbmap/reformat.sh in=$read1 in2=$read2 out=miRNA_trimmed/${sample}_interleaved.fq.gz"
		#bsub -M 4000 -J mrg -e mrg.2.err -o mrg.2.out "$cmd"
	fi

	if [ -s miRNA_trimmed/${sample}_interleaved.fq.gz ];then
		cmd="perl /rgs01/project_space/cab/automapper/common/tchang1/Softwares/sRNAtools/fq2collapedFa.pl -i miRNA_trimmed/${sample}_interleaved.fq.gz -o miRNA_trimmed/${sample}_interleaved.fasta.gz"
		#bsub -M 24000 -J mrg -e mrg.3.err -o mrg.3.out "$cmd"
		#$cmd

	fi

	if [ -s miRNA_trimmed/${sample}_interleaved.fasta.gz ]; then
		#cp /rgs01/project_space/cab/automapper/common/tchang1/Softwares/sRNAtools/DBCONFIG.txt .
		#cmdtest="echo -e $PERL5LIB"
		#species="hsa"
		species="mmu"	
		cmd="perl /rgs01/project_space/cab/automapper/common/tchang1/Softwares/sRNAtools/run_animal.pl -infile miRNA_trimmed/${sample}_interleaved.fasta.gz -species $species -ncpu 2 -outdir miRNA_sRNA_${species}/$sample"
		#$cmd
		bsub -q priority -M 24000 -J srna_${sample} -eo srna.${sample}.err -oo srna.${sample}.out "$cmd"
	fi

	

done 

