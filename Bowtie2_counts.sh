#!/bin/bash
#BSUB -P Bowtie2_COUNT
#BSUB -oo Bowtie2_COUNT.out -eo Bowtie2_COUNT.err
#BSUB -n 12
#BSUB -M 30G
#BSUB -J Bowtie2_COUNT

### Conda environ for latest bowtie2 ########.
source /hpcf/apps/conda3/install/5.1.0/bin/activate cutadaptenv
bowtie2 --quiet --align-paired-reads --min-score G,1,8 --local -D 20 -R 3 -N 1 -L 10 -i S,1,0.5 -p 10 -x /research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/GeneSeq/tRNA/Bowtie2.4.4/mm39-mature-tRNA/mature_tRNA -1 2113858_JBQS010_val_1.fq.gz -2 2113858_JBQS010_val_2.fq.gz


