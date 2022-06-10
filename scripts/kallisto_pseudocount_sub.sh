#!/bin/bash

#SBATCH --job-name=star_mapping
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=1G
#SBATCH --output=star_mapping.%J.out
#SBATCH --error=star_mapping.%J.err
#SBATCH --account=def-ben
#SBATCH --mail-user=cauretc@oregonstate.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load star/2.7.9a

#for indexing
python3 ~/project/cauretc/scripts/KO_rnaseq/kallisto_pseudocount.py $@ --threads 1
#for quantification
#python3 ~/project/cauretc/scripts/KO_rnaseq/kallisto_pseudocount.py $@ --threads 6

#sbatch ~/project/cauretc/scripts/KO_rnaseq/kallisto_pseudocount_sub.sh --ref_transcriptome /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla10_1_seq/XENLA_10.1_GCF.transcripts.fa.gz --step index
#sbatch ~/project/cauretc/scripts/KO_rnaseq/kallisto_pseudocount_sub.sh --ref_transcriptome /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla10_1_seq/XENLA_10.1_GCF.transcripts.fa.gz --step quant --fastq /home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data_3rd_run/scanw/*trim*.fastq.gz --output_directory /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/kallisto_dir --prefix_out
