# Intro
With BJE we are going to analyze the RNAseq data obtained from KO individuals of the genes on W of *X. laevis*.
Custom scripts to process and analyze data are in the `scripts` subfolder.

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

## Mapping reads

(March 9, 2022):

Star can directly count the reads per gene using by default only uniquely mapped genes - I am giving it a try since the results are expected to be the same as the default behavior of htseq. I am also outputting TranscriptomeSAM (to potentially used with salmon or stringtie). Some files did not have the same naming convention (ex. Run1 - ccdc_12_S3_L001_R1_trim_001.fastq.gz Run2 - ccdc11_S9_L001_R1_trim_001.fastq.gz) so Andrew helped with a python trick (try to split based on `_`, if errors with the expected number of splits, split the first string where the numbers start)

```
sbatch ~/project/cauretc/scripts/KO_rnaseq/star_mapping_sub.sh --fastq /home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data_3rd_run/scanw/*trim*.fastq.gz --genomeSTARDir /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/genomeDir --output_directory /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/star_mapping --prefix_out
```
Submitted for 10h but only needed 03h and 05min for the whole scanw dataset.

T27 failed - not sure why. Resubmitted by itself on March10.

```
sbatch ~/project/cauretc/scripts/KO_rnaseq/star_mapping_sub.sh --fastq /home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmw/*trim*.fastq.gz --genomeSTARDir /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/genomeDir --output_directory /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/star_mapping --prefix_out Run1_
```
Run 1 submitted for 6h, 2nd run for 4. 10h for ccdc Run 1 and 4h for ccdc Run 2.

I realized run1 had samples split on 2 lanes. On March 16, used `concat_PE_per_ind.py` to concatenate reads per individuals, then rerun star (dmw 3h30). Despite an error of slurm `sbatch: error: Batch job submission failed: Socket timed out on send/recv operation`, multiple jobs were submitted for ccdc, for sanity I resubmitted star on March 16 (time:).
```
python3 ~/project/cauretc/scripts/KO_rnaseq/concat_PE_per_ind.py --files /home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmw/*trim*.fastq.gz --output_directory /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/merged_fastq_files
```

### Checking for reads on the W specific genes:

#### Location on laevis 10.1

- dmw : Chr2L:182693575-182720555

- scanw (based on blast exons 1>6 including introns from v9): Chr2L:182715762-182725391

- ccdc69w (based on blast exons 1-2 + introns from v9): Chr2L:182743547-182743978 

#### Presence

(March 9-10, 2022)
```
module load StdEnv/2020  gcc/9.3.0 samtools/1.13
samtools index SCANW_T14_L002Aligned.sortedByCoord.out.bam
samtools tview -p Chr2L:182693575 SCANW_T14_L002Aligned.sortedByCoord.out.bam
samtools tview -p Chr2L:182715762 SCANW_T14_L002Aligned.sortedByCoord.out.bam
samtools tview -p Chr2L:182725391 SCANW_T14_L002Aligned.sortedByCoord.out.bam
samtools tview -p Chr2L:182743547 SCANW_T14_L002Aligned.sortedByCoord.out.bam

```
OK for dmw nothing for the others. BJE highlighted that because of low expression some individuals might not show those genes in the reads (same for T9).

Starting looking at the W specific region (a bit before cdk4 which starts around Chr2L:182818673 based on blast of the S copy):

- cdk6.L (ex around 182816699): a lot of good quality reads 

Going back to dmw and moving towards ccdc:

- 182777321: some high quality reads (no annotation on xenbase, same at 182777001)

**CCL:** I checked for scanw and ccdc (at least around the exons limits, more for 2 ind) for 3 wt ind: T14, T9 and T15 (others wt: 6, KO: 19, 27, 28, 30, 31). No reads except dmw. On the opposite KO T28 does not have reads for dmw while KO T6 has a few. I think it is enough for manual checking/first look at he mapping.

# Raw count - HTSEQ - not used
On ComputeCanada: 'HTSeq' framework, version 0.9.1
To load prior to be able to use commands such as `htseq-count`:
```
module load nixpkgs/16.09  gcc/5.4.0 mefit/1.0
```

Expected to give similar results as star. 

# Star output - count

Star output is in the form of 
```
N_unmapped      2151158 2151158 2151158
N_multimapping  1474468 1474468 1474468
N_noFeature     6838700 26316149        26227181
N_ambiguous     176505  36256   37036
XBmRNA83513     6       4       2
XBmRNA83514     2473    1271    1202
```
(Ex. from `SCANW_T9_L002ReadsPerGene.out.tab`)

Made `star_count_merged.py` to convert into a format ready for edgr. Ex. of the output format:

```
Genes, SCANW_T14_L002, SCANW_T15_L002, ...
XBmRNA83513, ..., ...,
XBmRNA83514, ..., ...,
```
Can be run like this:
```
 python3 ~/project/cauretc/scripts/KO_rnaseq/star_count_merged.py --files dmw*ReadsPerGene.out.tab --output_directory /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/star_mapping --prefix_out dmw_merged_all_ind_ReadsPerGene
```
# Raw count - stringtie

Online [manual](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#deseq). "Our focus in developing StringTie was on building a system that can assemble and quantitate transcripts regardless of whether gene annotation is available." (FAQ)

Stringtie requires a `XS` tag, need to specify for `star` `--outSAMstrandField intronMotif` (specified in the manual for compatibility with `Cufflinks/Cuffdiff`) or maybe simply `--outSAMattributes NH HI NM MD AS nM jM jI XS` (default: `NH HI AS nM`). Probably better to also do a `--twopassMode Basic`.

For each individual we will run `stringtie` with `-G <ref_ann.gff> -o [<path/>]<out.gtf> <read_alignments.bam>`, then use `stringtie --merge -G <ref_ann.gff> list_ind_file.txt -o output_merged.gtf`, follow by `stringtie -e -G output_merged.gtf -o [<path/>]<out_exp_est.gtf> <read_alignments.bam>` for each individual. Then, `prepDE.py -i sample_lst.txt`. Format of `sample_lst.txt`:
```
sample1 <PATH_TO_sample1.gtf>
sample2 <PATH_TO_sample2.gtf>
```

# Analysis of counts - EdgeR
## Installation
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("edgeR")
```
`export R_LIBS_USER=/home/cauretc/project/cauretc/R_libs` did not work to set the libs path but got installed in the good directory still (where I was)

## Analysis

Important: dmw was split over 2 runs, which needs to be accounted for. ccdc was also split over 2 runs but only males were in the 2nd runs. The most important here is to compare wt females versus KO females.
`edger_DE_KO.R` - need to be updated for it (March 25)