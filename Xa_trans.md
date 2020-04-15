## Download the transcriptome assembly
Publication: [Pownall et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5997840/)
```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111639/suppl/GSE111639_X.andrei_transcriptome_assembly.fasta.gz
```

## Make a database and run a blast search on cedar
```
sbatch ~/additional_scripts/makeblast_blast.sh ../GSE111639_X.andrei_transcriptome_assembly.fasta.gz X.andrei_transcriptome X.andrei_transcriptome_blastable dmw_Xa_Xlae_exons.fa dmw_Xa_Xlae_exons_Xa_transcriptome_blast.out
```
`makeblast_blast.sh`
```bash
#!/bin/bash

#SBATCH --account=def-ben  # The account to use
#SBATCH --time=00:02:00       # The duration in HH:MM:SS format
#SBATCH --cpus-per-task=1     # The number of cores
#SBATCH --mem=512M            # Total memory for this task

module load gcc/7.3.0 blast+/2.9.0

# Create the nucleotide database based on `ref.fa`.
gunzip -c $1 | makeblastdb -in - -title $2 -dbtype nucl -out $3
#blasting a query against the new database
blastn -evalue 1e-2 -query $4 -db $3 -out $5
```
