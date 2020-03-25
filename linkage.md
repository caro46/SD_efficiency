# R packages on cedar

```bash
module load  nixpkgs  gcc/7.3.0 r/3.6.1
#when doing module spider r/3.6.1 suggests a specific nixpkgs version but only one can be loaded ....
export R_LIBS=~/R/x86_64-pc-linux-gnu-library/3.6/

R
>install.packages("ggplot2")
>install.packages("dplyr")
>install.packages("qtl")
>install.packages("onemap")
>install.packages("vcfR")
# did it by steps to avoid cedar failing
```
# Linkage
Our crosses are usually outcrosses so need to use `Onemap` followed by `R/qtl`.

```R
library(onemap)
library(vcfR)

setwd("/home/cauretc/scratch/clivii_radseq_unknown_sex")
vcfR.object <- read.vcfR("mpileup_MQ20_varsOnly_cliviiFamily_chr1-2.vcf.gz")
vcf_example_out <- onemap_read_vcfR(vcfR.object = vcfR.object, parent1="all_Xc_BJE4530", parent2="all_Xc_BJE4531",cross = "outcross")
write_onemap_raw(vcf_example_out, file.name = "Xclivii_family_unknown_sex_chr1-2.raw", cross="outcross")
onemap_example_out_geno <- read_onemap(inputfile = "Xclivii_family_unknown_sex_chr1-2.raw")

onemap_example_out_geno <- read_onemap(inputfile = "Xclivii_family_all_info_chr1-2.raw")
```
```
cat Xclivii_family_unknown_sex_chr1-2.raw Xclivii_family_dmw_info_only_pheno.raw >Xclivii_family_all_info_chr1-2.raw

```
```
bins <- find_bins(onemap_example_out_geno, exact = FALSE)
bins_no_redudancy <- create_data_bins(onemap_example_out_geno, bins)

segreg_test <- test_segregation(bins_no_redudancy)
print(segreg_test)

dist <- select_segreg(segreg_test, distorted = TRUE, numbers = TRUE) #to show the markers numbers with segregation distortion
dist

no_dist <- select_segreg(segreg_test, distorted = FALSE, numbers = TRUE) #to show the markers numbers without segregation distortion
no_dist

##Estimating two-point recombination fractions
LOD_sug <- suggest_lod(bins_no_redudancy)
LOD_sug

twopts <- rf_2pts(bins_no_redudancy, LOD = LOD_sug, max.rf = 0.5)

##Assigning markers to linkage groups
mark_no_dist <- make_seq(twopts, c(no_dist))

marker_type(mark_no_dist)

LGs <- group(mark_no_dist)

print(LGs, detailed = FALSE)

CHR2L <- make_seq(twopts, "chr2L")
CHR2S <- make_seq(twopts, "chr2S")
CHR1S <- make_seq(twopts, "chr1S")
CHR1L <- make_seq(twopts, "chr1L")

CHR_mks <- group_seq(input.2pts = twopts, seqs = list(CHR1L = CHR1L,CHR2L = CHR2L, CHR1S = CHR1S, CHR2S=CHR2S),unlink.mks = mark_no_dist, repeated = FALSE)
CHR2L_frame <- mds_onemap(CHR_mks$sequences$CHR2L)

```
