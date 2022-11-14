# Intro
With BJE we are going to analyze the RNAseq data obtained from KO individuals of the genes on W of *X. laevis*.
Custom scripts to process and analyze data are in the `scripts` subfolder.

## Individuals info
```
SCANW_T14_L002  wt_f
SCANW_T15_L002  wt_f
SCANW_T19_L002  scanwKO
SCANW_T27_L002  scanwKO
SCANW_T28_L002  scanwKO
SCANW_T30_L002  scanwKO
SCANW_T31_L002  scanwKO
SCANW_T6_L002   wt_f
SCANW_T9_L002   wt_f
dmw_12_Run2_L001    wt_f
dmw_14_Run1_L001L002    dmwKO
dmw_15_Run2_L001    wt_f
dmw_16_Run1_L001L002    dmwKO
dmw_17_Run1_L001L002    wt_f
dmw_20_Run1_L001L002    wt_f
dmw_22_Run2_L001    wt_f
dmw_26_Run1_L001L002    dmwKO
dmw_28_Run1_L001L002    dmwKO
dmw_29_Run1_L001L002    dmwKO
dmw_35_Run1_L001L002    dmwKO
dmw_9_Run2_L001 wt_f
ccdc_11_Run2_L001   wt_m
ccdc_12_Run1_L001L002   wt_m
ccdc_13_Run2_L001   wt_m
ccdc_14_Run1_L001L002 ccdc_KO
ccdc_25_Run1_L001L002   wt_m
ccdc_2_Run2_L001    wt_m
ccdc_30_Run1_L001L002   ccdc_KO
ccdc_32_Run1_L001L002   ccdc_KO
ccdc_34_Run1_L001L002   wt_f
ccdc_35_Run1_L001L002   ccdc_KO
ccdc_36_Run1_L001L002   ccdc_KO
ccdc_3_Run1_L001L002    wt_f
ccdc_42_Run1_L001L002   ccdc_KO
ccdc_9_Run1_L001L002    wt_m
```
Names of individuals = headers of outputs from `star_count_merged.py`. Also contains lane and run information, important for dmw analysis of females wt vs females KO.

Made a file `/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/ind_KO_status.txt` with the same info.

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
```
sbatch ~/project/cauretc/scripts/KO_rnaseq/star_mapping_sub.sh --fastq /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/merged_fastq_files_Run1/ccdc*trim*.fastq.gz --genomeSTARDir /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/genomeDir --output_directory /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/star_mapping --prefix_out Run1_
```
(June 8th, 2022):

Deleting all the `star_mapping/*Aligned.toTranscriptome.out.bam` produced by `TranscriptomeSAM` since we are not going to use them. 
Sometime ago I added an option `count` in the `star_mapping.py`: for ex. `--count stringtie` will run with the needed option to obtain files compatible with stringtie (also do a twopass mode). `--count star` will do the commands from March 19.
```
#submitted for 10h - done in 08:37:12 
sbatch ~/project/cauretc/scripts/KO_rnaseq/star_mapping_sub.sh --fastq /home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data_3rd_run/scanw/*trim*.fastq.gz --genomeSTARDir /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/genomeDir --output_directory /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/star_mapping_for_stringtie --prefix_out --count stringtie

#submitted for 15h - done in
sbatch ~/project/cauretc/scripts/KO_rnaseq/star_mapping_sub.sh --fastq /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/merged_fastq_files_Run1/ccdc*trim*.fastq.gz --genomeSTARDir /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/genomeDir --output_directory /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/star_mapping_for_stringtie --prefix_out Run1_ --count stringtie
```
### Paths of the data
Since I had to merge only 1 Run, the files are in separate directories. To make it easier/quicker, below are the path of the reads data (trimmed only or trimmed and merged when needed), along with some information (number ind/sex - ! females refer to genotypic fem = ZW but not phenotypic !) 
```
#Run 1, ccdc69w, 11 ind. (3 wt males, 2 wt females, 6 KO females)
/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/merged_fastq_files_Run1/ccdc*trim*.fastq.gz
#Run 1, dmw, 8 ind. (2 wt females, 6 KO females)
/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/merged_fastq_files_Run1/dmw*trim*.fastq.gz
#Run 2, ccdc69w, 3 ind. (3 wt males)
/home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data_2nd_run/ccdc/*trim*.fastq.gz
#Run 2, dmw, 4 ind. (4 wt females)
/home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data_2nd_run/dmw/*trim*.fastq.gz
#Run 3, scanw, 9 ind. (4 wt females, 5 KO females)
/home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data_3rd_run/scanw/*trim*.fastq.gz
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
On graham: v. 2.1.5.

# Pseudo-alignment - kallisto

[Kallisto page](https://pachterlab.github.io/kallisto/starting). Command ex.:
```
kallisto index -i transcripts.idx transcripts.fasta.gz
kallisto quant -b 100 -i index -o output read1.fastq.gz read2.fastq.gz
```
(June 10th, 22)

```
sbatch ~/project/cauretc/scripts/KO_rnaseq/kallisto_pseudocount_sub.sh --ref_transcriptome /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla10_1_seq/XENLA_10.1_GCF.transcripts.fa.gz --step index
sbatch ~/project/cauretc/scripts/KO_rnaseq/kallisto_pseudocount_sub.sh --ref_transcriptome /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla10_1_seq/XENLA_10.1_GCF.transcripts.fa.gz --step quant --fastq /home/cauretc/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data_3rd_run/scanw/*trim*.fastq.gz --output_directory /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/kallisto_dir --prefix_out
```
`kallisto_pseudocount_sub.sh` needs to be updated for the resources needed depending on the step run.

For the indexing, I initially requested 1G of RAM, too low, 10G worked. It took 5min36, 1 thread. For quantification, BJE previously used 32G for his analysis, doing the same, except requesting 6 tasks (done in 01:28:12 for scanw, June 13th: dmw run1:, run2:, ccdc: run1:, run2:).

The `kallisto quant` produces a `.h5`, `.tsv` and `.json` files within each of the directory specify with the `-o`. For each samples they are name the same within each directory. BJE previously used a trinity utility script to merge all the individuals and obtain a better format.
```
#get the location of trinity on graham
module show trinity/2.14.0
#path to the needed script
#/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/trinity/2.14.0/trinityrnaseq-v2.14.0/util/abundance_estimates_to_matrix.pl
#get a file with the absolute path for each of the abundance file for scanw (within 'kallisto_dir')
ls -d $PWD/SCANW*/*.tsv >scanw_abundance_listing_target_files.txt
```
(July 19th, 2022)

Obtained all the `[...]_listing_target_files.txt` using the ex. command above.
To obtain the matrices, did not use sbatch since <2min of run. Requires to load EdgeR.
```
module load trinity/2.14.0
module load r/4.1.2
export R_LIBS_USER=/home/cauretc/project/cauretc/R_libs
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/trinity/2.14.0/trinityrnaseq-v2.14.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix dmw  --gene_trans_map none --name_sample_by_basedir --quant_files dmw_abundance_listing_target_files.txt
```
Ex. of out - ccdc69w.isoform.counts.matrix
```
ccdc_11_Run2_L001_kallisto_bout_out     ccdc_12_Run1_L001L002_kallisto_bout_out ccdc_13_Run2_L001_kallisto_bout_out     ccdc_14_Run1_L001L002_kallisto_bout_out ccdc_25_Run1_L001L002_kallisto_bout_out ccdc_2_Run2_L001_kallisto_bout_out      ccdc_30_Run1_L001L002_kallisto_bout_out ccdc_32_Run1_L001L002_kallisto_bout_out ccdc_34_Run1_L001L002_kallisto_bout_out ccdc_35_Run1_L001L002_kallisto_bout_out ccdc_36_Run1_L001L002_kallisto_bout_out ccdc_3_Run1_L001L002_kallisto_bout_out  ccdc_42_Run1_L001L002_kallisto_bout_out ccdc_9_Run1_L001L002_kallisto_bout_out
lcl|NC_054387.1_mrna_XM_018235619.2_76046       5.02266e-07     109     3.15872e-07     62      1.75048e-06     1.74321e-06     64.4464 127.407 125     0.0495367       126     7.51099e-08     3.07987e-07     0
```

# Analysis of counts - EdgeR
## Installation
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("edgeR")
```
`export R_LIBS_USER=/home/cauretc/project/cauretc/R_libs` did not work to set the libs path but got installed in the good directory still (where I was)

Version: 
```
[1] edgeR_3.36.0 limma_3.50.3
```

## Analysis

Important: dmw was split over 2 runs, which needs to be accounted for. ccdc was also split over 2 runs but only males were in the 2nd runs. The most important here is to compare wt females versus KO females.
`edger_DE_KO.R` - need to be updated for it (March 25). Updated April 13, 2022

### Using count from STAR
Nothing differentially expressed between wild type and KO females for dmw (required at least 5 reads for a gene to be considered). For scanw and ccdc69w, none of the identified DE genes are common and none of them are in Piprek et al. 2018. 

Updates July 29th: requiring on average 1 read/ind: nothing sig for dmw. 19/5 down/up for scanw, 3/7 for ccdc69w. Nothing common with Piprek.

### kallisto

During the analysis only the first part of the header was kept. I cannot directly used the annotation gff file as previously but all the information needed for annotation is contain in the header of the transcripts file. I put the info from header in a file as follow:
```
zgrep ">" /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla10_1_seq/XENLA_10.1_GCF.transcripts.fa.gz | sed -e 's/>//' > /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla10_1_seq/XENLA_10.1_GCF.transcripts_headers_only.txt
sed 's/[ ]\[/;/g ; s/\]//g' /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla10_1_seq/XENLA_10.1_GCF.transcripts_headers_only.txt >/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla10_1_seq/XENLA_10.1_GCF.transcripts_headers_only1.txt
mv /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla10_1_seq/XENLA_10.1_GCF.transcripts_headers_only1.txt /home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla10_1_seq/XENLA_10.1_GCF.transcripts_headers_only.txt
```
The `edger` analysis for kallisto output is within the `edger_DE_KO.R` after the `STAR` analysis.
Some overlap found for ccdc and scanw (mrtfa.S and paics.2.L) but not for dmw. paics is especially interesting ([Curzon et al. 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8268354/)).

Some other especially interesting genes:

- ccdc: LOC108704262 (product=oocyte zinc finger protein XlCOF29-like), tex261.L (testis expressed 261 L homeolog), zeb2.L (product=zinc finger E-box binding homeobox 2 L homeolog), paqr3.L (product=progestin and adipoQ receptor family member 3 L homeolog)

# DESeq2

```
 module load r/4.1.2 #same as previously for edgeR
```
```R
.libPaths("/home/cauretc/project/cauretc/R_libs")
BiocManager::install("DESeq2")
```
Needed to update some dependencies.

Version:

``` 
 [1] DESeq2_1.34.0               SummarizedExperiment_1.24.0
 [3] Biobase_2.54.0              MatrixGenerics_1.6.0       
 [5] matrixStats_0.62.0          GenomicRanges_1.46.1       
 [7] GenomeInfoDb_1.30.1         IRanges_2.28.0             
 [9] S4Vectors_0.32.4            BiocGenerics_0.40.0        
```
Run for ouputs from STAR on Aug. 3, 22. For analysis based on kallisto res., better to load data using `tximport`.

See `kallisto_DE_analysis.R` for analysis (made / Run Aug. 10) using `tximport` to import at the transcript and gene level count. For ccdc69w: 860/32, scanw: 128/14, dmw: 508/16 DE transcript/gene level with DeSeq2.

# Gene ontology

A lot of ways to do it. Some examples:

- overrepresentation of gene ontology categories using WebGestalt (Wang et al. 2013a) in Chandler et al. 2014

- enriched GO terms in ranked gene lists (GOrilla - Eden et al. 2009) in Wright et al. 2017

Xenbase has Data reports which includes "GO terms associated with Xenbase genepages". This contain the gene identifiers and the GO identifiers (ex. GO:0003723) associated with it. To retrieve GO term names associated with each GO IDs, we can use [YeastMine](https://yeastmine.yeastgenome.org/yeastmine/bag.do) as highlighted in the [FAQs of the geneontology webstite](http://geneontology.org/docs/faq/). I gave a try to directly use the go enrichment analysis on the [geneontology website](http://geneontology.org/) with as inputs the list of DE genes for scanw and dmw KO (DESeq2_wtVSKO_all_genes_05sign_gene_level.csv) obtained from kallisto + deseq2 analysis as well as the list I sent to BJE before (DESeq2_wtVSKO_all_genes_05sign_annotation_for_BJE.csv) - only used the genes for scanw - obtained from outputs of star; nothing showed up ("No statistically significant results"). [Young et al. 2010](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-2-r14#MOESM1) has good references (though it probably starts to be a bit old) and explanations in the intro + their application GOseq seems nice.  