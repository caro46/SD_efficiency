#!/bin/bash

#SBATCH --job-name=star_index_genome
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --time=01:00:00
#SBATCH --mem=100G
#SBATCH --output=star_index_genome.%J.out
#SBATCH --error=star_index_genome.%J.err
#SBATCH --account=def-ben
#SBATCH --mail-user=cauretc@oregonstate.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load star/2.7.9a

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir $1 \
--genomeFastaFiles $2 \
--sjdbGTFfile $3

#sbatch ~/project/cauretc/scripts/KO_rnaseq/genome_index.sh /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/genomeDir /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla_10_1_ref_genome/XENLA_10.1_genome.fa /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla_10_1_ref_genome/XENLA_10.1_GCF_XBmodels.gtf 
