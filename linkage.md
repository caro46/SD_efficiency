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
Below is a test. Need to make a better and optimized script.

```R
#libraries
library(onemap)
library(vcfR)

#working directory
setwd("/home/cauretc/scratch/clivii_radseq_unknown_sex")

#converting vcf to onemap format
vcfR.object_1_2 <- read.vcfR("mpileup_MQ20_varsOnly_cliviiFamily_chr1-2.vcf.gz")
vcfR.object_3_5 <- read.vcfR("mpileup_MQ20_varsOnly_cliviiFamily_chr3-5.vcf.gz")
vcfR.object_6_910 <- read.vcfR("mpileup_MQ20_varsOnly_cliviiFamily_chr6-910.vcf.gz")

vcf_example_out_1_2 <- onemap_read_vcfR(vcfR.object = vcfR.object_1_2, parent1="all_Xc_BJE4530", parent2="all_Xc_BJE4531",cross = "outcross")
vcf_example_out_3_5 <- onemap_read_vcfR(vcfR.object = vcfR.object_3_5, parent1="all_Xc_BJE4530", parent2="all_Xc_BJE4531",cross = "outcross")
vcf_example_out_6_910 <- onemap_read_vcfR(vcfR.object = vcfR.object_6_910, parent1="all_Xc_BJE4530", parent2="all_Xc_BJE4531",cross = "outcross")

write_onemap_raw(vcf_example_out_1_2, file.name = "Xclivii_family_unknown_sex_chr1-2.raw", cross="outcross")
write_onemap_raw(vcf_example_out_3_5, file.name = "Xclivii_family_unknown_sex_chr3-5.raw", cross="outcross")
write_onemap_raw(vcf_example_out_6_910, file.name = "Xclivii_family_unknown_sex_chr6-910.raw", cross="outcross")

onemap_example_out_geno_1_2 <- read_onemap(inputfile = "Xclivii_family_unknown_sex_chr1-2.raw")
onemap_example_out_geno_3_5 <- read_onemap(inputfile = "Xclivii_family_unknown_sex_chr3-5.raw")
onemap_example_out_geno_6_910 <- read_onemap(inputfile = "Xclivii_family_unknown_sex_chr6-910.raw")

#combine onemap
comb_geno_all_chr <- combine_onemap(onemap_example_out_geno_1_2, onemap_example_out_geno_3_5,onemap_example_out_geno_6_910)

write_onemap_raw(comb_geno_all_chr, file.name = "Xclivii_family_unknown_sex_all_chr.raw", cross="outcross")

onemap_example_out_geno <- read_onemap(inputfile = "Xclivii_family_all_info_all_chr.raw")

```
```
cat Xclivii_family_unknown_sex_all_chr.raw Xclivii_family_dmw_info_only_pheno.raw >Xclivii_family_all_info_all_chr.raw #and editing the ref to pheno (0 to 1)
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

CHR1S <- make_seq(twopts, "chr1S")
CHR2S <- make_seq(twopts, "chr2S")
CHR1L <- make_seq(twopts, "chr1L")
CHR2L <- make_seq(twopts, "chr2L")

CHR3S <- make_seq(twopts, "chr3S")
CHR4S <- make_seq(twopts, "chr4S")
CHR3L <- make_seq(twopts, "chr3L")
CHR4L <- make_seq(twopts, "chr4L")

CHR5S <- make_seq(twopts, "chr5S")
CHR6S <- make_seq(twopts, "chr6S")
CHR5L <- make_seq(twopts, "chr5L")
CHR6L <- make_seq(twopts, "chr6L")

CHR7S <- make_seq(twopts, "chr7S")
CHR8S <- make_seq(twopts, "chr8S")
CHR7L <- make_seq(twopts, "chr7L")
CHR8L <- make_seq(twopts, "chr8L")

CHR9_10S <- make_seq(twopts, "chr9_10S")
CHR9_10L <- make_seq(twopts, "chr9_10L")

CHR_mks <- group_seq(input.2pts = twopts, seqs = list(CHR1L = CHR1L,CHR2L = CHR2L, CHR3L = CHR3L, CHR4L = CHR4L, CHR5L = CHR5L, CHR6L = CHR6L, CHR7L = CHR7L, CHR8L = CHR8L, CHR9_10L = CHR9_10L, CHR1S = CHR1S, CHR2S=CHR2S, CHR3S = CHR3S, CHR4S = CHR4S, CHR5S = CHR5S, CHR6S = CHR6S, CHR7S = CHR7S, CHR8S = CHR8S, CHR9_10S = CHR9_10S),unlink.mks = mark_no_dist, repeated = FALSE)

#One or more of the provided marker sequences from list(CHR1L = CHR1L, CHR2L = CHR2L, CHR3L = CHR3L, CHR4L = CHR4L,      CHR5L = CHR5L, CHR6L = CHR6L, CHR7L = CHR7L, CHR8L = CHR8L,      CHR9_10L = CHR9_10L, CHR1S = CHR1S, CHR2S = CHR2S, CHR3S = CHR3S,      CHR4S = CHR4S, CHR5S = CHR5S, CHR6S = CHR6S, CHR7S = CHR7S,      CHR8S = CHR8S, CHR9_10S = CHR9_10S) do not form single linkage groups. The group with the highest number of markers belonging to the sequence will be considered.

CHR1L_frame <- mds_onemap(CHR_mks$sequences$CHR1L) 
CHR1S_frame <- mds_onemap(CHR_mks$sequences$CHR1S) 
CHR2L_frame <- mds_onemap(CHR_mks$sequences$CHR2L) #not working
CHR2S_frame <- mds_onemap(CHR_mks$sequences$CHR2S) #WORKING!!!
CHR3S_frame <- mds_onemap(CHR_mks$sequences$CHR3S) 
CHR3L_frame <- mds_onemap(CHR_mks$sequences$CHR3L)
CHR4S_frame <- mds_onemap(CHR_mks$sequences$CHR4S)
CHR4L_frame <- mds_onemap(CHR_mks$sequences$CHR4L)
CHR5S_frame <- mds_onemap(CHR_mks$sequences$CHR5S)
CHR5L_frame <- mds_onemap(CHR_mks$sequences$CHR5L)
CHR6S_frame <- mds_onemap(CHR_mks$sequences$CHR6S)
CHR6L_frame <- mds_onemap(CHR_mks$sequences$CHR6L)
CHR7S_frame <- mds_onemap(CHR_mks$sequences$CHR7S)
CHR7L_frame <- mds_onemap(CHR_mks$sequences$CHR7L)
CHR8S_frame <- mds_onemap(CHR_mks$sequences$CHR8S)
CHR8L_frame <- mds_onemap(CHR_mks$sequences$CHR8L)
CHR9_10S_frame <- mds_onemap(CHR_mks$sequences$CHR9_10S)
CHR9_10L_frame <- mds_onemap(CHR_mks$sequences$CHR9_10L)

rf_graph_table(CHR2S_frame) # graphic not showed

CHR2S_test_map <- map(CHR2S_frame)
#52 markers            log-likelihood: -1097.49

CHR1L_test_map <- map(CHR1L_frame) 
CHR1S_test_map <- map(CHR1S_frame) 
CHR2L_test_map <- map(CHR2L_frame) 
CHR3S_test_map <- map(CHR3S_frame) 
CHR3L_test_map <- map(CHR3L_frame)
CHR4S_test_map <- map(CHR4S_frame)
CHR4L_test_map <- map(CHR4L_frame)
CHR5S_test_map <- map(CHR5S_frame)
CHR5L_test_map <- map(CHR5L_frame)
CHR6S_test_map <- map(CHR6S_frame)
CHR6L_test_map <- map(CHR6L_frame)
CHR7S_test_map <- map(CHR7S_frame)
CHR7L_test_map <- map(CHR7L_frame)
CHR8S_test_map <- map(CHR8S_frame)
CHR8L_test_map <- map(CHR8L_frame)
CHR9_10L_test_map <- map(CHR9_10L_frame)
CHR9_10S_test_map <- map(CCHR9_10S_frame)

map1 <- list(CHR1S_test_map, CHR2S_test_map, CHR3S_test_map,CHR4S_test_map,CHR5S_test_map,CHR6S_test_map,CHR7S_test_map,CHR8S_test_map,CHR9_10S_test_map,CHR1L_test_map, CHR2L_test_map, CHR3L_test_map,CHR4L_test_map,CHR5L_test_map,CHR6L_test_map,CHR7L_test_map,CHR8L_test_map,CHR9_10L_test_map)

write_map(map1, "all_test_map.onemap.map")
```
`R/qtl`
```R
raw.file<-paste(system.file("example",package="onemap"),"fake.f2.onemap.raw", sep="/")
fake.f2.qtl <- read.cross("mm", file=raw.file, mapfile="fake.f2.onemap.map")
newmap <- est.map(fake.f2.qtl, tol=1e-6, map.function="kosambi")
plot.map(fake.f2.qtl, newmap)

fake.f2.qtl <- calc.genoprob(fake.f2.qtl, step=2)
out.em <- scanone(fake.f2.qtl, method="em")
out.hk <- scanone(fake.f2.qtl, method="hk")
plot(out.em, out.hk, col=c("blue","red"))

plot(out.hk, bandcol="gray70")

```
