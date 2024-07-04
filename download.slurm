#!/bin/sh
#SBATCH --ntasks=25 # Number of task
#SBATCH --time=80:00:00          # Run time in hh:mm:ss
#SBATCH --mem=80g       # Maximum memory required per CPU (in megabytes)
#SBATCH --nodes=1
#SBATCH --partition=batch
#SBATCH --job-name=WXS_TARGET
#SBATCH --error=/work/unmc_gudalab/nitish123/Osteosarcoma/WXS/WXS.err
#SBATCH --output=/work/unmc_gudalab/nitish123/Osteosarcoma/WXS/WXS.out
gdc-client download -m gdc_manifest_20200328_212052.txt -t gdc-user-token.2020-03-30T20_24_39.490Z.txt
