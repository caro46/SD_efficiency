# Intro
With BJE we are going to analyze the RNAseq data obtained from KO individuals of the genes on W of *X. laevis*

# Preliminary

## Xenbase data
Update on *X. laevis* assembly: W is assembled on 2L! dm-w, ccdc69.w and scan.w are all on 2L.
Blasted the ccdc69.w sequence on genome + Xenopus mRNA: no transcript, same for scan.w.

## Transcripts of W genes
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
