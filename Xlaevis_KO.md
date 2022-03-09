# Intro
With BJE we are going to analyze the RNAseq data obtained from KO individuals of the genes on W of *X. laevis*

# Preliminary

## Xenbase data
Update on *X. laevis* assembly: W is assembled on 2L! dm-w, ccdc69.w and scan.w are all on 2L.
Blasted the ccdc69.w sequence on genome + Xenopus mRNA: no transcript, same for scan.w. (Dec. 7, 2020)

## Transcripts of W genes

(Dec. 7-8, 2020)

One of the first things to do is to compare the scanw and ccdc69w transcripts sequences to those predicted by Mawaribuchi et al.
BJE produced trinity assemblies and used the dmw batch to get dmw and dmrt1 sequences so no need to produce the database: 
```bash
#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --account=def-ben
#SBATCH --output=blastn.%J.out
#SBATCH --error=blastn.%J.err
# Run: sbatch running_blast.sh ~/projects/rrg-ben/ben/Final_Sequences/XL_ccdc69w_ex1_and_ex2.fasta ~/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/trinity_assembly_all_batches_dmw/dmw_trinity_assembly_all_batches.Trinity.fasta_blastable ccdc_ex1_2_and_flanking_to_dmw_assemb.out
module load nixpkgs/16.09  gcc/7.3.0 blast+/2.9.0
blastn -query $1 -db $2 -out $3
```
Only one match for ccdc69w `TRINITY_DN767_c0_g2_i1` which actually correspond to ccdc69.L... scanw matches also correspond to autosomal genes.

`TRINITY_DN54529_c0_g1_i1` from the ccdc assembly seemed be promissing: tiny alignment but no snp and seemed to be a region specific to the W gene. 

Extracted the sequence

```awk -v seq="TRINITY_DN54529_c0_g1_i1" -v RS='>' '$1 == seq {print RS $0}' /home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/trinity_assembly_all_batches_ccdc/ccdc_trinity_assembly_all_batches_NEW.Trinity.fasta```

and blasted on xenbase: best hit for ccdc69.S except the tiny region from the transcript that is specific to the W... When blasting scanw on the ccdc assembly no perfect match except very small ones that are not specific to a region when reblasted on xenbase...

## Sex-linked genes in Xenbase transcriptome

(Dec. 10, 2020)

Piprek et al. 2018 identified important genes linked with sexual differentiation. If we want to use the xenbase transcriptome, we need to make sure that those are there. Sometimes 'a' or 'b' is specified for the paralogs, sometimes not despite the presence of both in the laevis genome (ex. amh) - I thus check at the gene level.

First using grep to identify which ones might need to be checked by blast.
```
grep -E 'gata4|sox9|dmrt1|amh|fgf9|ptgds|fshr|cyp17a1|xdm-w|fst|foxl2|cyp19a1' XENLA_10.1_GCF.transcripts.fa >Piprek_2018_genes_grep_Xenla10_1_trans.txt
```
There: gata4.L/S, sox9.L/S, dmrt1.L/S, amh.L/S, fgf9.S (no fgf9.L), ptgds.S (no ptgds.L), fshr.L/S, cyp17a1.L (no cyp17a1.S), xdm-w (=dm-w), fst.L/S, foxl2.L/S, cyp19a1.L (no cyp19a1.S)

fgf9.L is missing from transcriptome despite being annotated on the genome and expressed in the RNA-seq data when clicking on the gene on xenbase. Same for ptgds.L but not for cyp17a1.S and cyp19a1.S for which the paralog might actually be missing.

I also previously blasted the W genes sequence (dmw, scanw and ccdc69w) onto the downloaded transcriptome and only dmw was found (Dec. 8, 2020)

## Genome versus transcriptome

# Mapping against genome - STAR

`star/2.7.9a` installed and available as a module on graham.

## Indexing genome

```
wget https://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz
wget https://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1/XENLA_10.1_GCF_XBmodels.gtf

#submitted March8, 2022
sbatch ~/project/cauretc/scripts/KO_rnaseq/genome_index.sh /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/genomeDir /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla_10_1_ref_genome/XENLA_10.1_genome.fa /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla_10_1_ref_genome/XENLA_10.1_GCF_XBmodels.gtf
```

`genome_index.sh` can be found in the `scripts` directory.

Important note: the fasta file for the genome needs to be unzipped.

For the annotation file I used `XENLA_10.1_GCF_XBmodels.gtf` since already in the good format. Based on notes from Xenbase: NCBI-Xenbase gene models (gtf) - in theory to me the gene models should be better.

The indexing took about 35 min (with 6 threads). 

(March 9, 2022):

Star can directly count the reads per gene using by default only uniquely mapped genes - I am giving it a try since the results are expected to be the same as the default behavior of htseq. I am also outputting TranscriptomeSAM.

# Raw count - HTSEQ
On ComputeCanada: 'HTSeq' framework, version 0.9.1
To load prior to be able to use commands such as `htseq-count`:
```
module load nixpkgs/16.09  gcc/5.4.0 mefit/1.0
```

# Analysis of counts - EdgeR
## Installation
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("edgeR")
```
`export R_LIBS_USER=/home/cauretc/project/cauretc/R_libs` did not work to set the libs path but got installed in the good directory still (where I was)
