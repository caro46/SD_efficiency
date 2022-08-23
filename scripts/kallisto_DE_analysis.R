.libPaths("/home/cauretc/project/cauretc/R_libs")

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("rhdf5")

library(dplyr)
library(tidyr)
library(tximport)
library(DESeq2)
library(stringr)
library("rhdf5")

# table containing individual "names" and KO status
ind_KO_status <- read.table("/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/ind_KO_status.txt", sep="\t", header = FALSE)
colnames(ind_KO_status) <- c("name","KO")

# filter by KO gene
ind_KO_status_ccdc <- ind_KO_status %>% 
	filter(grepl("ccdc",ind_KO_status$name),!grepl("m",ind_KO_status$KO)) %>% 
	mutate(KO_status = factor(KO, levels = c("wt_f","ccdc_KO"))) %>%
	arrange(KO)

ind_KO_status_scan <- ind_KO_status %>% 
	filter(grepl("SCANW",ind_KO_status$name)) %>%
	mutate(KO_status = factor(KO, levels = c("wt_f","scanwKO"))) %>%
	arrange(KO) 

ind_KO_status_dmw <- ind_KO_status %>% 
	filter(grepl("dmw",ind_KO_status$name)) %>%
	separate(name, c("gene", "ind_numb", "run", "lane"), remove = FALSE) %>%
	select(name, run, KO) %>%
	mutate(run = factor(run, levels = c("Run1","Run2")), KO = factor(KO, levels = c("wt_f","dmwKO"))) %>%
	arrange(KO, run) 

# annotation from transcript header
laevis10transcript_info <- read.table("/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla10_1_seq/XENLA_10.1_GCF.transcripts_headers_only.txt", sep = ";", h = F, fill = T ) %>% 
	select(V1:V4) %>% dplyr::rename(transcript_ID = V1, gene_ID = V2, ref = V3, product= V4) %>%
	#mutate(gene_ID = str_remove(gene_ID, "gene=")) %>%
	select(transcript_ID, gene_ID)

# load the count files - transcript level
files_ccdc <- file.path("/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/kallisto_dir", paste(ind_KO_status_ccdc$name, "kallisto_bout_out", sep ="_"), "abundance.h5")
names(files_ccdc) <- ind_KO_status_ccdc$name
txi.kallisto_ccdc <- tximport(files_ccdc, type = "kallisto", txOut = TRUE)
head(txi.kallisto_ccdc$counts)

files_dmw <- file.path("/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/kallisto_dir", paste(ind_KO_status_dmw$name, "kallisto_bout_out", sep ="_"), "abundance.h5")
names(files_dmw) <- ind_KO_status_dmw$name
txi.kallisto_dmw <- tximport(files_dmw, type = "kallisto", txOut = TRUE)
head(txi.kallisto_dmw$counts)

files_scan <- file.path("/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/kallisto_dir", paste(ind_KO_status_scan$name, "kallisto_bout_out", sep ="_"), "abundance.h5")
names(files_scan) <- ind_KO_status_scan$name
txi.kallisto_scan <- tximport(files_scan, type = "kallisto", txOut = TRUE)
head(txi.kallisto_scan$counts)

# load the count files - gene level
txi.kallisto_dmw_gene <- tximport(files_dmw, type = "kallisto", tx2gene = laevis10transcript_info)
txi.kallisto_ccdc_gene <- tximport(files_ccdc, type = "kallisto", tx2gene = laevis10transcript_info)
txi.kallisto_scan_gene <- tximport(files_scan, type = "kallisto", tx2gene = laevis10transcript_info)

#######################
### DESeq2 analysis ###
#######################

ind_KO_status_ccdc <- ind_KO_status_ccdc %>%
	`rownames<-`(.$name) %>%
	select(KO)

ind_KO_status_scan <- ind_KO_status_scan %>%
	`rownames<-`(.$name) %>%
	select(KO)

ind_KO_status_dmw <- ind_KO_status_dmw %>%
	`rownames<-`(.$name) %>%
	select(KO, run)

## Transcript level

# checking if samples are in the same order in the metadata and the count tables
all(rownames(ind_KO_status_dmw) == colnames(txi.kallisto_dmw$counts)) 
all(rownames(ind_KO_status_ccdc) == colnames(txi.kallisto_ccdc$counts)) 
all(rownames(ind_KO_status_scan) == colnames(txi.kallisto_scan$counts)) 
#TRUE for all - good!

dds_scan <- DESeqDataSetFromTximport(txi.kallisto_scan, ind_KO_status_scan, ~KO) 
dds_scan <- dds_scan[rowSums(counts(dds_scan)) > ncol(counts(dds_scan)), ]
#from dim: 83536 9 to dim: 61536 9 
dds_ccdc <- DESeqDataSetFromTximport(txi.kallisto_ccdc, ind_KO_status_ccdc, ~KO) 
dds_ccdc <- dds_ccdc[rowSums(counts(dds_ccdc)) > ncol(counts(dds_ccdc)), ]
#from  dim: 83536 8 to dim: 60290 8
dds_dmw <- DESeqDataSetFromTximport(txi.kallisto_dmw, ind_KO_status_dmw, ~ run + KO) 
dds_dmw <- dds_dmw[rowSums(counts(dds_dmw)) > ncol(counts(dds_dmw)), ]
#from dim: 83536 12 to dim: 61035 12

dds_scan <- DESeq(dds_scan)
dds_dmw <- DESeq(dds_dmw)
dds_ccdc <- DESeq(dds_ccdc)

res05_scanw <- results(dds_scan, alpha=0.05) 
summary(res05_scanw)
sum(res05_scanw$padj < 0.05, na.rm=TRUE)
res05_scanw_Sig <- subset(res05_scanw, padj < 0.05) 
res05_scanw_Sig$KOgene = "scanw"

res05_ccdc <- results(dds_ccdc, alpha=0.05)
summary(res05_ccdc)
sum(res05_ccdc$padj < 0.05, na.rm=TRUE)
res05_ccdc_Sig <- subset(res05_ccdc, padj < 0.05)
res05_ccdc_Sig$KOgene = "ccdc69w"

res05_dmw <- results(dds_dmw, alpha=0.05)
summary(res05_dmw)
sum(res05_dmw$padj < 0.05, na.rm=TRUE)
res05_dmw_Sig <- subset(res05_dmw, padj < 0.05)
res05_dmw_Sig$KOgene = "dmw"

res05_scanw_dmw_ccdc <- rbind(res05_scanw_Sig,res05_ccdc_Sig,res05_dmw_Sig)
res05_scanw_dmw_ccdc_df <- as.data.frame(res05_scanw_dmw_ccdc) %>%
add_rownames(., var = "transcript_ID")
res05_scanw_dmw_ccdc_gene_ann <- left_join(res05_scanw_dmw_ccdc_df, laevis10transcript_info, by = "transcript_ID") 

write.csv(as.data.frame(res05_scanw_dmw_ccdc_gene_ann), 
          file="/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/kallisto_dir/DESeq2_wtVSKO_all_genes_05sign_transcript_level.csv") 

## Gene level

# checking if samples are in the same order in the metadata and the count tables
all(rownames(ind_KO_status_dmw) == colnames(txi.kallisto_dmw_gene$counts)) 
all(rownames(ind_KO_status_ccdc) == colnames(txi.kallisto_ccdc_gene$counts)) 
all(rownames(ind_KO_status_scan) == colnames(txi.kallisto_scan_gene$counts)) 
#TRUE!

dds_scan_gene <- DESeqDataSetFromTximport(txi.kallisto_scan_gene, ind_KO_status_scan, ~KO) 
dds_scan_gene <- dds_scan_gene[rowSums(counts(dds_scan_gene)) > ncol(counts(dds_scan_gene)), ]
#from dim: 19306 9 to dim: 16036 9
dds_ccdc_gene <- DESeqDataSetFromTximport(txi.kallisto_ccdc_gene, ind_KO_status_ccdc, ~KO) 
dds_ccdc_gene <- dds_ccdc_gene[rowSums(counts(dds_ccdc_gene)) > ncol(counts(dds_ccdc_gene)), ]
#from dim: 19306 8 to dim: 15999 8
dds_dmw_gene <- DESeqDataSetFromTximport(txi.kallisto_dmw_gene, ind_KO_status_dmw, ~ run + KO) 
dds_dmw_gene <- dds_dmw_gene[rowSums(counts(dds_dmw_gene)) > ncol(counts(dds_dmw_gene)), ]
#from dim: 19306 12 to dim: 15983 12 

dds_scan_gene <- DESeq(dds_scan_gene)
dds_dmw_gene <- DESeq(dds_dmw_gene)
dds_ccdc_gene <- DESeq(dds_ccdc_gene)

res05_scanw_gene <- results(dds_scan_gene, alpha=0.05) 
summary(res05_scanw_gene)
res05_scanw_Sig_gene <- subset(res05_scanw_gene, padj < 0.05) 
res05_scanw_Sig_gene$KOgene = "scanw"

res05_ccdc_gene <- results(dds_ccdc_gene, alpha=0.05)
summary(res05_ccdc_gene)
res05_ccdc_Sig_gene <- subset(res05_ccdc_gene, padj < 0.05)
res05_ccdc_Sig_gene$KOgene = "ccdc69w"

res05_dmw_gene <- results(dds_dmw_gene, alpha=0.05)
summary(res05_dmw_gene)
res05_dmw_Sig_gene <- subset(res05_dmw_gene, padj < 0.05)
res05_dmw_Sig_gene$KOgene = "dmw"
res05_scanw_dmw_ccdc_gene <- rbind(res05_scanw_Sig_gene,res05_ccdc_Sig_gene,res05_dmw_Sig_gene)

write.csv(as.data.frame(res05_scanw_dmw_ccdc_gene), 
          file="/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/kallisto_dir/DESeq2_wtVSKO_all_genes_05sign_gene_level.csv")

#######################
### EdgeR analysis - not run !###
#######################

#Creating object with group indicator
dge_scanw = DGEList(txi.kallisto_scan_gene$counts, group = ind_KO_status_scan$KO)
dge_ccdc = DGEList(txi.kallisto_ccdc_gene$counts, group = ind_KO_status_ccdc$KO)
dge_dmw <- DGEList(txi.kallisto_dmw_gene$counts, group = ind_KO_status_dmw$KO)

# normalize the data
dge_scanw <- calcNormFactors(dge_scanw)
dge_ccdc <- calcNormFactors(dge_ccdc)
dge_dmw <- calcNormFactors(dge_dmw)

# design matrix for dmw only since individuals on 2 separate runs
design_dmw <- model.matrix(~0+batch_dmw+group_dmw)
rownames(design_dmw) <- colnames(dge_dmw)

# estimate dispersion
dge_scanw <- estimateDisp(dge_scanw)
dge_ccdc <- estimateDisp(dge_ccdc)
dge_dmw <- estimateDisp(dge_dmw, design_dmw, robust=TRUE)

# find genes differentially expressed between the two groups
et_scanw <- exactTest(dge_scanw)
et_ccdc <- exactTest(dge_ccdc)

#Fit model for dmw
fit_dmw <- glmQLFit(dge_dmw, design_dmw, robust=TRUE)
qlf_dmw <- glmQLFTest(fit_dmw) #by default the test is for the last coefficient in the design matrix

# Res summary
summary(decideTests(et_scanw))
summary(decideTests(et_ccdc))
summary(decideTests(qlf_dmw))

et_scanw$table$padj <- p.adjust(et_scanw$table$PValue, method="BH")
et_ccdc$table$padj <- p.adjust(et_ccdc$table$PValue, method="BH")
qlf_dmw$table$padj <- p.adjust(qlf_dmw$table$PValue, method="BH")

#enriched mutants
scanw_DE_genes_sig_enriched <- et_scanw$table[et_scanw$table$padj < 0.05 & et_scanw$table$logFC > 0,]
scanw_DE_genes_sig_enriched$change = "+"
ccdc_DE_genes_sig_enriched <- et_ccdc$table[et_ccdc$table$padj < 0.05 & et_ccdc$table$logFC > 0,]
ccdc_DE_genes_sig_enriched$change = "+"
dmw_DE_genes_sig_enriched <- qlf_dmw$table[qlf_dmw$table$padj < 0.05 & qlf_dmw$table$logFC > 0,]
dmw_DE_genes_sig_enriched$change = "+"

#depleted mutants
scanw_DE_genes_sig_depleted <- et_scanw$table[et_scanw$table$padj < 0.05 & et_scanw$table$logFC < 0,]
scanw_DE_genes_sig_depleted$change = "-"
ccdc_DE_genes_sig_depleted <- et_ccdc$table[et_ccdc$table$padj < 0.05 & et_ccdc$table$logFC < 0,]
ccdc_DE_genes_sig_depleted$change = "-"
dmw_DE_genes_sig_depleted <- qlf_dmw$table[qlf_dmw$table$padj < 0.05 & qlf_dmw$table$logFC < 0,]
dmw_DE_genes_sig_depleted$change = "-"

scanw_DE_genes_sig <- rbind(scanw_DE_genes_sig_enriched, scanw_DE_genes_sig_depleted) 
scanw_DE_genes_sig <- scanw_DE_genes_sig %>%
  mutate(geneID = rownames(scanw_DE_genes_sig))

dmw_DE_genes_sig <- rbind(dmw_DE_genes_sig_enriched, dmw_DE_genes_sig_depleted) 
dmw_DE_genes_sig <- dmw_DE_genes_sig %>%
  mutate(geneID = rownames(dmw_DE_genes_sig))

ccdc_DE_genes_sig <- rbind(ccdc_DE_genes_sig_enriched, ccdc_DE_genes_sig_depleted) 
ccdc_DE_genes_sig <- ccdc_DE_genes_sig %>%
  mutate(geneID = rownames(ccdc_DE_genes_sig))

ccdc_DE_genes_sig <- left_join(ccdc_DE_genes_sig, laevis10transcript_info, by = "transcript_ID") 
dmw_DE_genes_sig <- left_join(dmw_DE_genes_sig, laevis10transcript_info, by = "transcript_ID") 
scanw_DE_genes_sig <- left_join(scanw_DE_genes_sig, laevis10transcript_info, by = "transcript_ID") 
