#!/bin/bash

#SBATCH --job-name=star_mapping
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --output=star_mapping.%J.out
#SBATCH --error=star_mapping.%J.err
#SBATCH --account=def-ben
#SBATCH --mail-user=cauretc@oregonstate.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load star/2.7.9a

~/project/cauretc/scripts/KO_rnaseq/star_mapping.py $@ --threads 6

#sbatch ~/project/cauretc/scripts/KO_rnaseq/star_mapping_sub.sh --fastq /home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data_3rd_run/scanw/*trim*.fastq.gz --genomeSTARDir /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/genomeDir --output_directory /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/star_mapping --prefix_out
