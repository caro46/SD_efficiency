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
To work interactively on cedar:
```
salloc --time 2:0:0 --ntasks=1 --account=[account_name] --mem 2G
```

# Onemap Linkage
Our crosses are usually outcrosses so need to use `Onemap` followed by `R/qtl`.
Below is a test. Need to make a better and optimized script.

```R
#libraries
library(onemap)
library(vcfR)

#working directory
setwd("[path_to_scratch]/scratch/clivii_radseq_unknown_sex")

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

#CHR1L_frame <- mds_onemap(CHR_mks$sequences$CHR1L) 
#CHR1S_frame <- mds_onemap(CHR_mks$sequences$CHR1S) 
#CHR2L_frame <- mds_onemap(CHR_mks$sequences$CHR2L) #not working
#CHR2S_frame <- mds_onemap(CHR_mks$sequences$CHR2S) #WORKING!!!
#CHR3S_frame <- mds_onemap(CHR_mks$sequences$CHR3S) 
#CHR3L_frame <- mds_onemap(CHR_mks$sequences$CHR3L)
#CHR4S_frame <- mds_onemap(CHR_mks$sequences$CHR4S)
#CHR4L_frame <- mds_onemap(CHR_mks$sequences$CHR4L)
#CHR5S_frame <- mds_onemap(CHR_mks$sequences$CHR5S)
#CHR5L_frame <- mds_onemap(CHR_mks$sequences$CHR5L)
#CHR6S_frame <- mds_onemap(CHR_mks$sequences$CHR6S)
#CHR6L_frame <- mds_onemap(CHR_mks$sequences$CHR6L)
#CHR7S_frame <- mds_onemap(CHR_mks$sequences$CHR7S)
#CHR7L_frame <- mds_onemap(CHR_mks$sequences$CHR7L)
#CHR8S_frame <- mds_onemap(CHR_mks$sequences$CHR8S)
#CHR8L_frame <- mds_onemap(CHR_mks$sequences$CHR8L)
#CHR9_10S_frame <- mds_onemap(CHR_mks$sequences$CHR9_10S)
#CHR9_10L_frame <- mds_onemap(CHR_mks$sequences$CHR9_10L)

#Since issue with mds_onemap for 2L, tried the other way and worked so let's try for everybody

CHR1L_ord <- order_seq(CHR_mks$sequences$CHR1L) 
CHR1L_frame <- make_seq(CHR1L_ord, "force")
CHR1S_ord <- order_seq(CHR_mks$sequences$CHR1S) 
CHR1S_frame <- make_seq(CHR1S_ord, "force")

CHR2L_ord <- order_seq(CHR_mks$sequences$CHR2L) 
CHR2L_frame <- make_seq(CHR2L_ord, "force")
CHR2S_ord <- order_seq(CHR_mks$sequences$CHR2S) 
CHR2S_frame <- make_seq(CHR2S_ord, "force")

CHR3L_ord <- order_seq(CHR_mks$sequences$CHR3L) 
CHR3L_frame <- make_seq(CHR3L_ord, "force")
CHR3S_ord <- order_seq(CHR_mks$sequences$CHR3S) 
CHR3S_frame <- make_seq(CHR3S_ord, "force")

CHR4L_ord <- order_seq(CHR_mks$sequences$CHR4L) 
CHR4L_frame <- make_seq(CHR4L_ord, "force")
CHR4S_ord <- order_seq(CHR_mks$sequences$CHR4S) 
CHR4S_frame <- make_seq(CHR4S_ord, "force")

CHR5L_ord <- order_seq(CHR_mks$sequences$CHR5L) 
CHR5L_frame <- make_seq(CHR5L_ord, "force")
CHR5S_ord <- order_seq(CHR_mks$sequences$CHR5S) 
CHR5S_frame <- make_seq(CHR5S_ord, "force")

CHR6L_ord <- order_seq(CHR_mks$sequences$CHR6L) 
CHR6L_frame <- make_seq(CHR6L_ord, "force")
CHR6S_ord <- order_seq(CHR_mks$sequences$CHR6S) 
CHR6S_frame <- make_seq(CHR6S_ord, "force")

CHR7L_ord <- order_seq(CHR_mks$sequences$CHR7L) 
CHR7L_frame <- make_seq(CHR7L_ord, "force")
CHR7S_ord <- order_seq(CHR_mks$sequences$CHR7S) 
CHR7S_frame <- make_seq(CHR7S_ord, "force")

CHR8L_ord <- order_seq(CHR_mks$sequences$CHR8L) 
CHR8L_frame <- make_seq(CHR8L_ord, "force")
CHR8S_ord <- order_seq(CHR_mks$sequences$CHR8S) 
CHR8S_frame <- make_seq(CHR8S_ord, "force")

CHR9_10L_ord <- order_seq(CHR_mks$sequences$CHR9_10L) 
CHR9_10L_frame <- make_seq(CHR9_10L_ord, "force")
CHR9_10S_ord <- order_seq(CHR_mks$sequences$CHR9_10S) 
CHR9_10S_frame <- make_seq(CHR9_10S_ord, "force")

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
# Plinks

## plink2

```bash
bcftools query -l mpileup_MQ20_varsOnly_cliviiFamily_chr3-5.vcf.gz
bcftools concat mpileup_MQ20_varsOnly_cliviiFamily_chr1-2.vcf.gz mpileup_MQ20_varsOnly_cliviiFamily_chr3-5.vcf.gz mpileup_MQ20_varsOnly_cliviiFamily_chr6-910.vcf.gz -o mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat.vcf
gzip mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat.vcf
sbatch [path_to_project]project/[account]/scripts/plink.sh mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat.vcf.gz pheno_binary.plink all_chr_plink
```
Example of beginning of `.plink` file
```
#IID    dmw_e2
all_Xc_BJE4530  1
all_Xc_BJE4531  2
Xc_BE10_boy.fastq.gz    2
```
```R
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

setwd("[path_to_scratch]/scratch/clivii_radseq_unknown_sex")

Xcl_dmw_gwas_plinks <- read.table("./all_chr_plink.dmw_e2.glm.logistic",h=F,sep="\t")
head(Xcl_dmw_gwas_plinks)

colnames(Xcl_dmw_gwas_plinks) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "OR", "LOG(OR)_SE", "Z_STAT", "LOG10_P")

main_ggplot <- ggplot(data = Xcl_dmw_gwas_plinks) + 
  geom_point(aes(x=POS/1000000, #to put the axis into Mb instead of exponents
                 y=LOG10_P))+
  # black & white color theme 
  theme(axis.text.x = element_text(colour = "black",size=15),
        axis.text.y = element_text(colour = "black",size=15),
        axis.title=element_text(size=16,face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())+
  labs(x = "position (Mb)") +
  facet_grid(cols=vars(CHROM))

ggsave("Rplot_Xcl_dmw_gwas_plinks_all_chr_logP.png", plot = main_ggplot, dpi = 200)
```
## plink1
```
bcftools reheader --samples header_change_ID_XCl.txt -o mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat_no_underscoreID.vcf ../mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat.vcf

plink --vcf mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat_no_underscoreID.vcf --make-bed --out mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat --allow-extra-chr
```
The previous command creates `.bed`, `.bim` and `.fam` files.
```
plink --bed mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat.bed --fam Pedigree_pheno_no_underscore.txt --bim mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat.bim --geno 0.05 --maf 0.05 --hwe 0.000001 --tdt --ci 0.95 --allow-extra-chr --allow-no-sex

plink --bed mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat.bed --fam Pedigree_pheno_no_underscore.txt --bim mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat.bim --geno 0.05 --maf 0.1 --tdt --ci 0.95 --allow-extra-chr --allow-no-sex
```
Kept only the plot from the last command. Both similar. Signal on multiple unexpected chromosomes.

Beginning of `header_change_ID_XCl.txt` to rename without the underscores.
```
all_Xc_BJE4530  allXcBJE4530
all_Xc_BJE4531  allXcBJE4531
Xc_BE10_boy.fastq.gz    XcBE10boy.fastq.gz
Xc_BE11_boy.fastq.gz    XcBE11boy.fastq.gz
```
`Pedigree_pheno_no_underscore.txt`
```
Family1 allXcBJE4530    0       0       2       2
Family1 allXcBJE4531    0       0       1       1
Family1 XcBE10boy.fastq.gz      allXcBJE4531    allXcBJE4530    0       1
```
# Fst
```bash
module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.16vcftools --remove-indv Xc_BE3_girl.fastq.gz --gzvcf mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat.vcf.gz --recode --out mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat_no_BE3.vcf.gz
vcftools --remove-indv Xc_BE3_girl.fastq.gz --remove-indv Xc_U4_ukn.fastq.gz --remove-indv Xc_BE9_girl.fastq.gz --remove-indv Xc_BE12_girl.fastq.gz --gzvcf mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat.vcf.gz --recode --out mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat_no_BE3_U4_BE9_BE12.vcf.gz
vcf-to-tab < mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat_no_BE3_U4_BE9_BE12.vcf.gz.recode.vcf > mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat_no_BE3_U4_BE9_BE12.vcf.gz.recode.tab
sbatch ~/additional_scripts/vcftools.sh ../mpileup_MQ20_varsOnly_cliviiFamily_all_chr_concat_no_BE3_U4_BE9_BE12.vcf.gz.recode.vcf dmw_pop.txt no_dmw_pop.txt vcftools_fst_Xcliv_uknSex_dmw_no_dmw_50kb_wind.weir.fst
```
```R
setwd("[path_to_scratch]/scratch/clivii_radseq_unknown_sex/fst")

Xcl_dmw_fst <- read.table("./vcftools_fst_Xcliv_uknSex_dmw_no_dmw.weir.fst",h=T,sep="\t")

cleanPlot2 <- function() {theme_bw()+
    theme(
      panel.spacing = unit(0.2, "lines"),
      strip.background = element_blank(),
      strip.text.y = element_text(angle = 0, size = 16),
      strip.text.x = element_text(size = 16),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      #legend.position = 'top',
      axis.line = element_blank(),
      #axis.ticks = element_blank(),
      panel.border = element_rect(fill = NA,colour = "grey20"),
      panel.grid = element_blank(),
      #panel.grid.major = element_line(color = 'grey85'),
      #panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
}

main_ggplot <- ggplot(data = Xcl_dmw_fst) + 
  geom_point(aes(x=POS/1000000, #to put the axis into Mb instead of exponents
                 y=WEIR_AND_COCKERHAM_FST),size=1)+
  labs(x = "position (Mb)", y = "Fst (dmw/no dmw)") +
  facet_grid(rows=vars(CHROM)) +
  cleanPlot2()

ggsave("Rplot_Xcl_dmw_no_dmw_fst.png", plot = main_ggplot, dpi = 200)
```
# Conclusion
Not enough data to make ccl.
