.libPaths("/home/cauretc/project/cauretc/R_libs")
#.libPaths()
#install.packages('dplyr',lib='/home/cauretc/project/cauretc/R_libs')
#install.packages('tidyr',lib='/home/cauretc/project/cauretc/R_libs')
#install.packages('statmod',lib='/home/cauretc/project/cauretc/R_libs')
library(edgeR)
library(dplyr)
library(tidyr)
library(statmod)

#in addition to the manual
#good info: https://research.stowers.org/cws/CompGenomics/Projects/edgeR.html

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}

#working directory
setwd("/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/star_mapping")

#Read the annotation files
laevis10gtf <- read.table("/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla_10_1_ref_genome/XENLA_10.1_GCF_XBmodels.gtf", skip = 9,
                          sep = "\t", h = F, col.names = c("chromosome","source","feature","start","end","score","strand","phase",
                                                           "gene_info"))
laevis10gtf$geneID <- unlist(lapply(laevis10gtf$gene_info , extract_attributes, "gene_id"))
laevis10gtf$geneName <- unlist(lapply(laevis10gtf$gene_info , extract_attributes, "gene_name"))

laevis10gtf_min_info <- laevis10gtf %>% select(geneID, geneName) %>%
  distinct()
#column_names = c("SCANW_T14_L002", "SCANW_T15_L002", "SCANW_T19_L002", "SCANW_T27_L002", "SCANW_T28_L002", "SCANW_T30_L002",
#                 "SCANW_T31_L002", "SCANW_T6_L002", "SCANW_T9_L002")

#counts_scanw <- read.table("SCANW_merged_all_ind_ReadsPerGene.tab", sep="\t", skip=1, row.names = 1) %>%
#  select(-V11) 
#colnames(counts_scanw) = column_names

#Read the count files (count obtained with STAR)
counts_scanw <- read.table("SCANW_merged_all_ind_ReadsPerGene.tab", sep="\t", h = TRUE, row.names = 1)
counts_dmw <- read.table("dmw_merged_all_ind_ReadsPerGene.tab", sep="\t", h = TRUE, row.names = 1)
#for ccdc we also have wt males, for now, exluce them
ccdc_female_KO_wt <- c("ccdc_14_Run1_L001L002", "ccdc_30_Run1_L001L002", "ccdc_32_Run1_L001L002", "ccdc_35_Run1_L001L002", 
                       "ccdc_36_Run1_L001L002", "ccdc_42_Run1_L001L002",
                       "ccdc_34_Run1_L001L002", "ccdc_3_Run1_L001L002") #only 2 wt_f
counts_ccdc_female_KO_wt <- read.table("ccdc_merged_all_ind_ReadsPerGene.tab", sep="\t", h = TRUE, row.names = 1) %>%
  select(all_of(ccdc_female_KO_wt))
  
head(counts_scanw)
dim(counts_scanw)

# Examine correlation among data sets
cor(counts_scanw)

# we can visualize the correlation
pdf("scan_first_look_heatmap.pdf")
heatmap(cor(counts_scanw))
dev.off()

#Based on https://research.stowers.org/cws/CompGenomics/Projects/edgeR.html - good to filter very low read counts
counts_scanw <- counts_scanw %>% filter(rowSums(.) > 5)
counts_dmw <- counts_dmw %>% filter(rowSums(.) > 5)
counts_ccdc_female_KO_wt <- counts_ccdc_female_KO_wt %>% filter(rowSums(.) > 5)

#from 
#dim(counts_scanw) #same number of genes for each set (scanw, ccdc and dmw)
#[1] 44478     9
#to
#dim(counts_scanw)
#[1] 32283     9
#dim(counts_ccdc_female_KO_wt)
#[1] 31475     8
#dim(counts_dmw)
#[1] 32528    12

#Making sure the wt_f is the ref
group_scanw = factor(c("wt_f","wt_f","scanwKO","scanwKO","scanwKO","scanwKO","scanwKO","wt_f","wt_f")) %>%
  relevel(., ref = "wt_f")
group_dmw = factor(c("wt_f","dmwKO","wt_f","dmwKO","wt_f","wt_f","wt_f","dmwKO","dmwKO", "dmwKO", "dmwKO", "wt_f"))%>%
  relevel(., ref = "wt_f")
batch_dmw <- factor(c("Run2","Run1","Run2","Run1","Run1","Run1","Run2","Run1","Run1", "Run1", "Run1", "Run2"))
group_ccdc = factor(c(rep("ccdc_KO",6),rep("wt_f", 2)))%>%
  relevel(., ref = "wt_f")

#Creating object with group indicator
dge_scanw = DGEList(counts_scanw, group = group_scanw)
dge_dmw = DGEList(counts_dmw, group = group_dmw)
dge_ccdc = DGEList(counts_ccdc_female_KO_wt, group = group_ccdc)
dge_scanw
dge_dmw
dge_ccdc

# normalize the data
dge_scanw <- calcNormFactors(dge_scanw)
dge_dmw <- calcNormFactors(dge_dmw)
dge_ccdc <- calcNormFactors(dge_ccdc)

#design matrix for dmw only since individuals on 2 separate runs
#design_dmw <- model.matrix(~batch_dmw+group_dmw)
design_dmw <- model.matrix(~0+batch_dmw+group_dmw)
rownames(design_dmw) <- colnames(dge_dmw)

# estimate dispersion
dge_scanw <- estimateDisp(dge_scanw)
#dge_dmw <- estimateDisp(dge_dmw)
dge_dmw <- estimateDisp(dge_dmw, design_dmw, robust=TRUE)
dge_ccdc <- estimateDisp(dge_ccdc)

# Examine samples
#cols=c(rep("black",2), rep("red",5),rep("black",2))
#pdf("scan_first_look_MDS.pdf")
#plotMDS(dge_scanw, col=cols, cex=0.5)
#dev.off()

# find genes differentially expressed between the two groups
et_scanw <- exactTest(dge_scanw)
#et_dmw <- exactTest(dge_dmw)
et_ccdc <- exactTest(dge_ccdc)

#Fit model for dmw
fit_dmw <- glmQLFit(dge_dmw, design_dmw, robust=TRUE)
qlf_dmw <- glmQLFTest(fit_dmw) #by default the test is for the last coefficient in the design matrix
#FDR <- p.adjust(qlf_dmw$table$PValue, method="BH")
#sum(FDR < 0.05)
#[1] 0

# get the top DE genes
topTags(et_scanw)
topTags(et_ccdc)
topTags(qlf_dmw)
summary(decideTests(qlf_dmw))
#group_dmwdmwKO
#Down                0
#NotSig          32528
#Up                  0
summary(decideTests(et_scanw))
#scanwKO-wt_f
#Down             23
#NotSig        32254
#Up                6
summary(decideTests(et_ccdc))
#ccdc_KO-wt_f
#Down              3
#NotSig        31466
#Up                6

### How many genes look significant?
sum(et_scanw$table$PValue < 0.05)
#[1] 1280
### How many genes show 2-fold enrichment?
sum(et_scanw$table$PValue < 0.05 & et_scanw$table$logFC > 1)
#[1] 413
sum(et_scanw$table$PValue < 0.05 & et_scanw$table$logFC < -1)
#[1] 527

#Check on dmrt1
#dmrt1.L XBXL10_1g2070 
#dmrt1.S XBXL10_1g4848

goi <- c(dmrt1.L="XBXL10_1g2070",dmrt1.S="XBXL10_1g4848",foxl2.L="XBXL10_1g24241", foxl2.S="XBXL10_1g26060")
#goi <- c(dmrt1.L="XBXL10_1g2070",dmrt1.S="XBXL10_1g4848")
iv_scanw <- match(goi, rownames(et_scanw$table))
et_scanw$table[iv_scanw,]
#logFC    logCPM     PValue
#XBXL10_1g2070   1.0078066 -0.2314961 0.2431423
#NA                     NA         NA        NA
#XBXL10_1g24241  3.5422747 -2.1836900 0.0212517
#XBXL10_1g26060 -0.5087939 -3.3017051 0.8465553

# adjust p-values and assign the result to our table
et_scanw$table$padj <- p.adjust(et_scanw$table$PValue, method="BH")
#et_dmw$table$padj <- p.adjust(et_dmw$table$PValue, method="BH")
et_ccdc$table$padj <- p.adjust(et_ccdc$table$PValue, method="BH")
sum(et_scanw$table$padj < 0.05 & et_scanw$table$logFC > 1)
#[1] 6
sum(et_scanw$table$padj < 0.05 & et_scanw$table$logFC < -1)
#[1] 23

#enriched mutant
#before relevel
#scanw_DE_genes_sig_enriched <- et_scanw$table[et_scanw$table$padj < 0.05 & et_scanw$table$logFC < -1,]
#dmw_DE_genes_sig_enriched <- et_dmw$table[et_dmw$table$padj < 0.05 & et_dmw$table$logFC < -1,]
#ccdc_DE_genes_sig_enriched <- et_ccdc$table[et_ccdc$table$padj < 0.05 & et_ccdc$table$logFC < -1,]
scanw_DE_genes_sig_enriched <- et_scanw$table[et_scanw$table$padj < 0.05 & et_scanw$table$logFC > 0,]
ccdc_DE_genes_sig_enriched <- et_ccdc$table[et_ccdc$table$padj < 0.05 & et_ccdc$table$logFC > 0,]
#scanw_DE_genes_sig_enriched_two_fold <- et_scanw$table[et_scanw$table$padj < 0.05 & et_scanw$table$logFC > 1,] #same without 2fold
#ccdc_DE_genes_sig_enriched_two_fold <- et_ccdc$table[et_ccdc$table$padj < 0.05 & et_ccdc$table$logFC > 1,] #same without 2fold

#depleted mutants
#before relevel
#scanw_DE_genes_sig_depleted <- et_scanw$table[et_scanw$table$padj < 0.05 & et_scanw$table$logFC > 1,]
#dmw_DE_genes_sig_depleted <- et_dmw$table[et_dmw$table$padj < 0.05 & et_dmw$table$logFC > 1,]
#ccdc_DE_genes_sig_depleted <- et_ccdc$table[et_ccdc$table$padj < 0.05 & et_ccdc$table$logFC > 1,]
scanw_DE_genes_sig_depleted <- et_scanw$table[et_scanw$table$padj < 0.05 & et_scanw$table$logFC < 0,]
ccdc_DE_genes_sig_depleted <- et_ccdc$table[et_ccdc$table$padj < 0.05 & et_ccdc$table$logFC < 0,]
#scanw_DE_genes_sig_depleted_two_fold <- et_scanw$table[et_scanw$table$padj < 0.05 & et_scanw$table$logFC < -1,] #same without 2fold
#ccdc_DE_genes_sig_depleted_two_fold <- et_ccdc$table[et_ccdc$table$padj < 0.05 & et_ccdc$table$logFC < -1,] #same without 2fold

#combining DE analysis and gene names
scanw_DE_genes_sig <- rbind(scanw_DE_genes_sig_enriched, scanw_DE_genes_sig_depleted) 
scanw_DE_genes_sig <- scanw_DE_genes_sig %>%
  mutate(geneID = rownames(scanw_DE_genes_sig))
#laevis10gtf_min_info %>% filter(geneID %in% rownames(DE_genes_sig_depleted))
scanw_DE_genes_sig <- left_join(scanw_DE_genes_sig, laevis10gtf_min_info, by = "geneID") 

#dmw_DE_genes_sig <- rbind(dmw_DE_genes_sig_enriched, dmw_DE_genes_sig_depleted) 
#dmw_DE_genes_sig <- dmw_DE_genes_sig %>%
#  mutate(geneID = rownames(dmw_DE_genes_sig))
#dmw_DE_genes_sig <- left_join(dmw_DE_genes_sig, laevis10gtf_min_info, by = "geneID") 

ccdc_DE_genes_sig <- rbind(ccdc_DE_genes_sig_enriched, ccdc_DE_genes_sig_depleted) 
ccdc_DE_genes_sig <- ccdc_DE_genes_sig %>%
  mutate(geneID = rownames(ccdc_DE_genes_sig))
ccdc_DE_genes_sig <- left_join(ccdc_DE_genes_sig, laevis10gtf_min_info, by = "geneID") 

head(scanw_DE_genes_sig)
#head(dmw_DE_genes_sig)
head(ccdc_DE_genes_sig)
dim(scanw_DE_genes_sig)
#[1] 29  6
dim(ccdc_DE_genes_sig)
#[1] 9 6
#Saving 
write.csv(ccdc_DE_genes_sig,"/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/EdgeR_analysis/DE_genes_ccdc_KO.csv", 
          row.names = FALSE, quote=FALSE)
write.csv(scanw_DE_genes_sig,"/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/EdgeR_analysis/DE_genes_scanw_KO.csv", 
          row.names = FALSE, quote=FALSE)

#checking if any gene from Piprek et al. 2018
genes_interest <- c("gata4.L","sox9.L","dmrt1.L","amh.L","fgf9.L","ptgds.L","fshr.L","cyp17a1.L","xdm-w","fst.L","foxl2.L","cyp19a1.L","gata4.S","sox9.S","dmrt1.S","amh.S","fgf9.S","ptgds.S","fshr.S","cyp17a1.S","xdm-w","fst.S","foxl2.S","cyp19a1.S")
ccdc_DE_genes_sig %>% filter(geneName %in% genes_interest)
scanw_DE_genes_sig %>% filter(geneName %in% genes_interest)
#>nothing

#anything common between ccdc and scanw
scanw_ccdc_shared <- inner_join(scanw_DE_genes_sig, ccdc_DE_genes_sig, by = "geneID")
#ccdc_DE_genes_sig %>% filter(geneName %in% scanw_DE_genes_sig$geneName)
#scanw_DE_genes_sig %>% filter(geneName %in% ccdc_DE_genes_sig$geneName)
#<0 rows> (or 0-length row.names)
#>nothing

####OLD
#Some results before filtering for a min of 5 reads for a gene:

#dge_scanw <- calcNormFactors(dge_scanw)

#dge_scanw$samples
#                 group lib.size norm.factors
#SCANW_T14_L002    wt_f 36342950    0.9753863
#SCANW_T15_L002    wt_f 31114374    0.8065953
#SCANW_T19_L002 scanwKO 36868149    0.8095117
#SCANW_T27_L002 scanwKO 33647428    1.1214668
#SCANW_T28_L002 scanwKO 32650831    1.0690811
#SCANW_T30_L002 scanwKO 33260050    1.1912418
#SCANW_T31_L002 scanwKO 40925387    0.8812362
#SCANW_T6_L002     wt_f 33130843    1.1225273
#SCANW_T9_L002     wt_f 38793411    1.1113689

#topTags(et_scanw)

#Comparison of groups:  scanwKO- wt_f 
#logFC   logCPM       PValue          FDR
#XBXL10_1g34381 -6.386599 3.145556 2.609778e-21 1.160777e-16 #cel.2.L, NM_001086315.2
#XBXL10_1g21885 -6.499912 4.728659 6.134821e-21 1.364323e-16 #amy2b.S, NM_001086441.1
#XBXL10_1g10780 -3.433434 4.975135 6.975193e-20 1.034142e-15 #endoul2.S, NM_001110753.1
#XBXL10_1g30134 -3.418415 5.890410 2.483507e-19 2.761535e-15 #acta2.L, NM_001091337.1
#XBXL10_1g12880 -6.196385 5.613784 4.004002e-19 3.561800e-15 #cpa1.L, NM_001095030.2
#XBXL10_1g30046 -6.660248 3.765290 7.068575e-18 5.239935e-14 #LOC100036900, XM_018224202.2
#XBXL10_1g29967 -5.841578 2.036482 3.647646e-16 2.317714e-12 #LOC108695981, XM_018224940.2
#XBXL10_1g19493 -6.201585 3.262667 2.502314e-15 1.240078e-11 #amy2a.L, NM_001096169.1
#XBXL10_1g25971 -6.599612 5.370523 2.509263e-15 1.240078e-11 #cpb1.S, NM_001095031.2
#XBXL10_1g30041 -6.184325 3.286144 1.983284e-14 8.821250e-11 #MGC64575, NM_001086279.1

#sum(et_scanw$table$PValue < 0.05)
#[1] 1926
### How many genes show 2-fold enrichment?
#sum(et_scanw$table$PValue < 0.05 & et_scanw$table$logFC > 1)
#[1] 610
#sum(et_scanw$table$PValue < 0.05 & et_scanw$table$logFC < -1)
#[1] 1165

#goi <- c(dmrt1.L="XBXL10_1g2070",dmrt1.S="XBXL10_1g4848",foxl2.L="XBXL10_1g24241", foxl2.S="XBXL10_1g26060")
#iv_scanw <- match(goi, rownames(et_scanw$table))
#et_scanw$table[iv_scanw,]
#logFC    logCPM     PValue
#XBXL10_1g2070 -0.4479935 -0.231144 0.56892950
#XBXL10_1g4848 -4.5947434 -3.795928 0.06952234

#sum(et_scanw$table$padj < 0.05 & et_scanw$table$logFC > 1)
#[1] 17
#sum(et_scanw$table$padj < 0.05 & et_scanw$table$logFC < -1)
#[1] 192

#scanw_DE_genes_sig %>% filter(geneName %in% genes_interest)
#1 -2.263684 2.733762 6.306143e-05 0.01908059 XBXL10_1g37486  ptgds.S #not after filtering based on 5 reads for each gene
