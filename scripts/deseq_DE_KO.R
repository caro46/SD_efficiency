.libPaths("/home/cauretc/project/cauretc/R_libs")
#BiocManager::install("apeglm")
#BiocManager::install("tximport")
library(DESeq2)
library(dplyr)
library(tidyr)
library("apeglm")

# Using the results from STAR (need to use tximport for res from kallisto)
#setwd("/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/kallisto_dir/trinity_matrix_format/")
setwd("/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/star_mapping")

#counts_scanw <- read.table("scanw.isoform.counts.matrix", sep="\t", h = TRUE, row.names = 1)
#colnames(counts_scanw)<-gsub("_kallisto_bout_out","",colnames(counts_scanw))
#counts_dmw <- read.table("dmw.isoform.counts.matrix", sep="\t", h = TRUE, row.names = 1)
#colnames(counts_dmw)<-gsub("_kallisto_bout_out","",colnames(counts_dmw))

#for ccdc we also have wt males, for now, exluce them
#ccdc_female_KO_wt <- c("ccdc_14_Run1_L001L002", "ccdc_30_Run1_L001L002", "ccdc_32_Run1_L001L002", "ccdc_35_Run1_L001L002", 
#                       "ccdc_36_Run1_L001L002", "ccdc_42_Run1_L001L002",
#                       "ccdc_34_Run1_L001L002", "ccdc_3_Run1_L001L002") #only 2 wt_f
#counts_ccdc_female_KO_wt <- read.table("ccdc69w.isoform.counts.matrix", sep="\t", h = TRUE, row.names = 1)
#colnames(counts_ccdc_female_KO_wt)<-gsub("_kallisto_bout_out","",colnames(counts_ccdc_female_KO_wt))

#counts_ccdc_female_KO_wt <- counts_ccdc_female_KO_wt %>%
#  select(all_of(ccdc_female_KO_wt))
#counts_ccdc_female_KO_wt <- as.matrix(counts_ccdc_female_KO_wt)

counts_scanw <- read.table("SCANW_merged_all_ind_ReadsPerGene.tab", sep="\t", h = TRUE, row.names = 1)
counts_dmw <- read.table("dmw_merged_all_ind_ReadsPerGene.tab", sep="\t", h = TRUE, row.names = 1)
#for ccdc we also have wt males, for now, exluce them
ccdc_female_KO_wt <- c("ccdc_14_Run1_L001L002", "ccdc_30_Run1_L001L002", "ccdc_32_Run1_L001L002", "ccdc_35_Run1_L001L002", 
                       "ccdc_36_Run1_L001L002", "ccdc_42_Run1_L001L002",
                       "ccdc_34_Run1_L001L002", "ccdc_3_Run1_L001L002") #only 2 wt_f
counts_ccdc_female_KO_wt <- read.table("ccdc_merged_all_ind_ReadsPerGene.tab", sep="\t", h = TRUE, row.names = 1) %>%
  select(all_of(ccdc_female_KO_wt))
  

# Metadata containing info about individuals, KO status, run and lanes
metadata_KO <- read.table("/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/ind_KO_status.txt", sep="\t", h = FALSE) 
colnames(metadata_KO) <- c("ind_info", "KO_status") 

## Separate by gene and rearrange (wt = first)
scanw_meta <- metadata_KO %>% filter(., grepl('SCANW', ind_info)) %>%
	mutate(KO_status = factor(KO_status, levels = c("wt_f","scanwKO"))) %>%
	arrange(KO_status) %>%
	`rownames<-`(.$ind_info) %>%
	select(KO_status)

dmw_meta <- metadata_KO %>% filter(., grepl('dmw', ind_info)) %>%
	separate(ind_info, c("gene", "ind_numb", "run", "lane"), remove = FALSE) %>%
	select(ind_info, run, KO_status) %>%
	mutate(run = factor(run, levels = c("Run1","Run2")), KO_status = factor(KO_status, levels = c("wt_f","dmwKO"))) %>%
	arrange(KO_status, run) %>%
	`rownames<-`(.$ind_info) %>%
	select(KO_status, run)

ccdc_meta <- metadata_KO %>% filter(., grepl('ccdc', ind_info)) %>%
  filter(ind_info %in% ccdc_female_KO_wt) %>% 
	mutate(KO_status = factor(KO_status, levels = c("wt_f","ccdc_KO"))) %>%
	arrange(KO_status) %>%
	`rownames<-`(.$ind_info) %>%
	select(KO_status)

## Same order for metadata and expression data
counts_scanw <- counts_scanw[, rownames(scanw_meta)]
all(rownames(scanw_meta) == colnames(counts_scanw)) 
#[1] TRUE
counts_dmw <- counts_dmw[, rownames(dmw_meta)]
all(rownames(dmw_meta) == colnames(counts_dmw))
#[1] TRUE
counts_ccdc_female_KO_wt <- counts_ccdc_female_KO_wt[, rownames(ccdc_meta)]
all(rownames(ccdc_meta) == colnames(counts_ccdc_female_KO_wt))
#[1] TRUE

# Filtering low expressed genes
counts_scanw <- counts_scanw %>% filter(rowSums(.) > ncol(.))
counts_dmw <- counts_dmw %>% filter(rowSums(.) > ncol(.))
counts_ccdc_female_KO_wt <- counts_ccdc_female_KO_wt %>% filter(rowSums(.) > ncol(.))

counts_ccdc_female_KO_wt <- as.matrix(counts_ccdc_female_KO_wt)
counts_scanw <- as.matrix(counts_scanw) 
counts_dmw <- as.matrix(counts_dmw)

# DESeqDataSetFromMatrix expect non negative integer
deseq_scanw <- DESeqDataSetFromMatrix(countData = counts_scanw,
                              colData = scanw_meta,
                              design = ~ KO_status)
deseq_scanw
deseq_ccdc <- DESeqDataSetFromMatrix(countData = counts_ccdc_female_KO_wt,
                              colData = ccdc_meta,
                              design = ~ KO_status)
deseq_ccdc
deseq_dmw <- DESeqDataSetFromMatrix(countData = counts_dmw,
                              colData = dmw_meta,
                              design = ~ run + KO_status)
deseq_dmw

deseq_scanw <- DESeq(deseq_scanw)
deseq_dmw <- DESeq(deseq_dmw)
deseq_ccdc <- DESeq(deseq_ccdc)

#res_scanw <- results(deseq_scanw)
#res_scanw
#schrincage
#resultsNames(deseq_scanw)
#resLFC_scanw <- lfcShrink(deseq_scanw, coef="KO_status_scanwKO_vs_wt_f", type="apeglm")
#resLFC_scanw

# FDR cutoff of .05 (default .1)
## Significant res

#resOrdered_scanw <- res_scanw[order(res_scanw$pvalue),]
#summary(res_scanw)
#sum(res_scanw$padj < 0.1, na.rm=TRUE)
res05_scanw <- results(deseq_scanw, alpha=0.05) 
summary(res05_scanw)
sum(res05_scanw$padj < 0.05, na.rm=TRUE)
res05_scanw_Sig <- subset(res05_scanw, padj < 0.05) 
res05_scanw_Sig$KOgene = "scanw"

res05_ccdc <- results(deseq_ccdc, alpha=0.05)
summary(res05_ccdc)
sum(res05_ccdc$padj < 0.05, na.rm=TRUE)
res05_ccdc_Sig <- subset(res05_ccdc, padj < 0.05)
res05_ccdc_Sig$KOgene = "ccdc69w"

res05_dmw <- results(deseq_dmw, alpha=0.05)
summary(res05_dmw)
sum(res05_dmw$padj < 0.05, na.rm=TRUE)
res05_dmw_Sig <- subset(res05_dmw, padj < 0.05)
res05_dmw_Sig$KOgene = "dmw"

# Saving significant result for all KO together
res05_scanw_dmw_ccdc <- rbind(res05_scanw_Sig,res05_ccdc_Sig,res05_dmw_Sig)
write.csv(as.data.frame(res05_scanw_dmw_ccdc), 
          file="/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/DESeq2_analysis/DESeq2_wtVSKO_all_genes_05sign.csv")

res05_scanw_dmw_ccdc_df <- as.data.frame(res05_scanw_dmw_ccdc) %>%
add_rownames(., var = "geneID")

# Add annotation from gtf file

## function to parse gtf files
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}

## laevis v10 gtf file
laevis10gtf <- read.table("/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/Xenla_10_1_ref_genome/XENLA_10.1_GCF_XBmodels.gtf", skip = 9,
                          sep = "\t", h = F, col.names = c("chromosome","source","feature","start","end","score","strand","phase",
                                                           "gene_info"))
laevis10gtf$geneID <- unlist(lapply(laevis10gtf$gene_info , extract_attributes, "gene_id"))
laevis10gtf$geneName <- unlist(lapply(laevis10gtf$gene_info , extract_attributes, "gene_name"))

laevis10gtf_min_info <- laevis10gtf %>% select(geneID, geneName) %>%
  distinct()

## left join so we add info but don't lose locus with no annotation 
res05_scanw_dmw_ccdc_df <- left_join(res05_scanw_dmw_ccdc_df, laevis10gtf_min_info, by = "geneID") %>%
mutate(upordown = case_when (log2FoldChange < 0 ~ "-",
	TRUE ~ "+") )

# Save result keeping all previous info
write.csv(as.data.frame(res05_scanw_dmw_ccdc_df), 
          file="/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/DESeq2_analysis/DESeq2_wtVSKO_all_genes_05sign_annotation.csv")

# keeping only column requested by BJE
res05_scanw_dmw_ccdc_df_for_BJE <- res05_scanw_dmw_ccdc_df %>% select(KOgene, geneName, upordown)
write.csv(as.data.frame(res05_scanw_dmw_ccdc_df_for_BJE), 
          file="/home/cauretc/projects/rrg-ben/cauretc/2021_KO_rnaseq/STAR_dir/DESeq2_analysis/DESeq2_wtVSKO_all_genes_05sign_annotation_for_BJE.csv")