## ---------------------------
##
## Script name: 01_gene_counts.R
##
## Purpose of script: 
## This script performs the analysis and generates the figures 
## for comparing the number of genes identified between tape-strips
## and biopsies
## 
##
## Author: Blaine Fritz and Isabel Diaz-Pines Cort 
##
## Email: bgfritz@sund.ku.dk



#-------------------------------------------------------------------------------
# Setup ------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Clear workspace --------------------------------------------------------------
rm(list = ls())


# Check that output directories exist and, if not, create them -----------------

out_dir <- "./results"

if (!dir.exists(out_dir)) {dir.create(out_dir)}

dirs <- c("generated_data", 
          "generated_figures",
          "generated_rds") 

for(dir in file.path(out_dir, dirs)){
  if (!dir.exists(dir)) {dir.create(dir)}
}


# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(reshape2)
library(ggpubr)
# Define Functions -------------------------------------------------------------

#Function to count genes above cutoff X in each sample

gene_count_vs_cutoff <- function(counts, metadata, x){
  count_cutoff <- tibble(gene_count = as.numeric(colSums(counts >= x)),
                         Sample = names(colSums(counts >= x))) %>%
    inner_join(metadata, by = "Sample")
  count_cutoff$cutoff <- rep(x, times=length(count_cutoff$Sample))
  return(count_cutoff)
}

# Function to perform a one sample t-test for the difference in number of identified 
# genes between methods

pval_vs_cutoff <- function(gene_count_diff, x){
  diff <- gene_count_diff %>%
    filter(cutoff == x) %>%
    dplyr::select(diff)
  t_test <- t.test(diff, mu=0)
  n_gene.diff <- tibble (cutoff = x,
                         pval = t_test$p.value,
                         mean_diff = t_test$estimate)
  return(n_gene.diff)
}

## Function to identify how many genes are shared between each tape-biopsy pair

find_common_genes <- function(metadata, counts, cutoff, x){ # x is Participant nr
  # find  tape and biopsy samples for current Participant
  tape <- metadata$Sample[metadata$Participant == x & metadata$Sample_Type == "tape_strip"]
  biopsy <- metadata$Sample[meta$Participant == x & metadata$Sample_Type == "biopsy"]
  # find genes detected by each method above certain cutoff and intersect 
  detected_tape <- rownames (counts[counts[tape] >= cutoff,])
  detected_biopsy <- rownames (counts[counts[biopsy] >= cutoff,])
  common_genes <- intersect(detected_tape, detected_biopsy)
  n_common_genes <- tibble(n_common = length(common_genes),
                           n_tape = length(detected_tape),
                           n_biopsy = length(detected_biopsy),
                           Participant = as.character(x),
                           cutoff = cutoff)
  return(n_common_genes)
}

# Function to find genes detected by both methods in each patient (using raw counts)
# and then calculate correlation of expression levels with the norm counts

get_corr_common_genes <- function(meta, raw_counts, norm_counts, cutoff, x){ # x is patient nr
  
  samples_by_patient <- meta %>% 
    dplyr::select(Participant, Sample) %>%
    group_by(Participant) %>%
    dplyr::mutate(row = row_number()) %>% # create a unique id row for each patient (avoids getting list-cols)
    pivot_wider(names_from = Participant, values_from = Sample) %>%
    dplyr::select(-row)
  
  # for each patient find tape and biopsy samples
  sample_1 <- as.character(samples_by_patient[1, x])
  sample_2 <- as.character(samples_by_patient[2, x])
  # which genes are detected above 10 in each of the samples?
  detected_1 <- rownames (raw_counts[raw_counts[sample_1] >= cutoff,])
  detected_2 <- rownames (raw_counts[raw_counts[sample_2] >= cutoff,])
  # find common genes
  common_genes <- intersect(detected_1, detected_2)
  # necessary bc after filtering lowly expressed genes (less than 70% samples with count >10) 
  # some of the common genes are removed
  common_genes <- intersect(common_genes, rownames(norm_counts))
  # get correlation coeff from normalized counts
  corr <- cor(norm_counts[common_genes, sample_1],
              norm_counts[common_genes, sample_2],
              method = "pearson")
  corr_df <- tibble(corr = corr,
                    n_common_genes = length(common_genes),
                    patient_nr = as.factor(x),
                    cutoff = cutoff)
  return(corr_df)
}


getCodingGenes <- function() {
  
  require(biomaRt)
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                  dataset = "hsapiens_gene_ensembl", 
                  host = 'https://www.ensembl.org')
  genes <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                         "external_gene_name",
                                         "chromosome_name",
                                         "transcript_biotype"),
                          filters = c("transcript_biotype","chromosome_name"),
                          values = list("protein_coding",c(1:22)),
                          mart = mart)
  
  return(genes)
  
}



#-------------------------------------------------------------------------------
# Import the data --------------------------------------------------------------
#-------------------------------------------------------------------------------

meta <- read_tsv("data/Metadata_Final.tsv", 
                 col_types = list(Sample_Type = col_factor(levels = c("tape_strip", "biopsy")))
                 )

count <- read_tsv("data/count_matrix.tsv")

#Import and reorder counts
count <- column_to_rownames(count, "gene")
count <- count[, meta$Sample]

# remove genes with all zero counts
count <- count[rowSums(count)>0,]


#-------------------------------------------------------------------------------
# Data Sets --------------------------------------------------------------------
#-------------------------------------------------------------------------------

## Quantify # of genes with counts over X for each patient ---------------------
x <- seq(1,500)
df_gene_count <- lapply(x, gene_count_vs_cutoff, counts=count, metadata=meta)
df_gene_count <- bind_rows(df_gene_count)

## Calculate difference between Biopsy and Tape-Strip for each participant with 
## each cutoff value (0,10,100)

df_gene_count_diff <- df_gene_count %>%
  dplyr::select(gene_count, Participant, Sample_Type, cutoff, Eczema_dorsalhand, AD) %>%
  pivot_wider(names_from = Sample_Type, values_from = gene_count) %>%
  mutate(diff = biopsy - tape_strip) 

## Test whether the difference in number of genes between Biopsy and Tape-Strip 
## at each cutoff is significantly different from zero (one-sample t-test)

x <- seq(1,500)
df_pval_vs_cutoff <- lapply(x, pval_vs_cutoff, gene_count_diff=df_gene_count_diff)
df_pval_vs_cutoff <- bind_rows(df_pval_vs_cutoff) 
df_pval_vs_cutoff$pval_adj <- p.adjust(df_pval_vs_cutoff$pval, method = "fdr") #Adjust for multiple testing


# Identify how many genes are shared between the methods for cutoff values 
# of 10, 75, and 200 

x <- meta$Participant
n_common_genes_10 <- lapply(x, find_common_genes,
                            metadata = meta,
                            counts = count,
                            cutoff = 10)
n_common_genes_75 <- lapply(x, find_common_genes,
                             metadata = meta,
                             counts = count,
                             cutoff = 75)
n_common_genes_200 <- lapply(x, find_common_genes,
                             metadata = meta,
                             counts = count,
                             cutoff = 200)

df_n_common_genes <- bind_rows(n_common_genes_10, n_common_genes_75, n_common_genes_200)
df_n_common_genes <- melt(df_n_common_genes, id.vars = c("Participant", "cutoff"))


# VST transformed raw count data 

# Pre-filtering 
smallestGroupSize <- 5 
keep <- rowSums(count >= 10) >= smallestGroupSize
count_filtered <- count[keep,]


df_vst_count <- DESeq2::vst(as.matrix(count_filtered))



# Calculate Correlations between genes shared among each tape-biopsy pair

x <- seq(1,12)
df_corr <- lapply(x, get_corr_common_genes,
                  meta = meta,
                  raw_counts=count_filtered,
                  norm_counts=df_vst_count,
                  cutoff=10)

df_corr <- bind_rows(df_corr)



#Export all of the data sets 
write_csv(df_gene_count, paste0(out_dir, "/generated_data/01_df_n_genes.csv"))
write_csv(rownames_to_column(count_filtered, "gene"), paste0(out_dir, "/generated_data/01_df_count_filtered.csv"))
write_csv(rownames_to_column(data.frame(df_vst_count), "gene"), 
          paste0(out_dir, "/generated_data/01_df_count_vst.csv"))
write_csv(df_gene_count_diff, paste0(out_dir, "/generated_data/01_df_gene_count_diff.csv"))
write_csv(df_pval_vs_cutoff, paste0(out_dir, "/generated_data/01_df_pval_vs_cutoff.csv"))
write_csv(df_n_common_genes, paste0(out_dir, "/generated_data/01_df_n_common_genes.csv"))
write_csv(df_corr, paste0(out_dir, "/generated_data/01_df_corr.csv"))
write_csv(filter(df_n_common_genes, cutoff == 10), paste0(out_dir, "/generated_data/01_n_common_genes_cutoff10.csv"))


#-------------------------------------------------------------------------------
# Figures ----------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Define colors for plots ------------------------------------------------------

sample_type.cols <- c("steelblue1", "indianred3")
names(sample_type.cols) <- levels(meta$Sample_Type)

AD.cols <- c("tomato1", "forestgreen")
names(AD.cols) <- levels(meta$AD)

eczema.cols <- c("darkorange", "lightgreen", "lightpink1")
names(eczema.cols) <- levels(meta$Eczema_dorsalhand)



# Boxplot of number of detected genes above 10 and 100 count cutoff ------------

fig_n_genes_method <- df_gene_count %>%
  filter(cutoff %in% c(10, 75, 200)) %>%
ggplot(mapping = aes(Sample_Type, 
                     gene_count,
                              fill = Sample_Type)) +
  geom_boxplot(width = 0.5, outlier.colour = NULL) +
  geom_jitter(size=2, alpha=0.5, width = 0.1) +
  geom_line(aes(group=Participant), colour="grey", linetype="11") +
  facet_wrap(~ cutoff, labeller = labeller(cutoff = c("10" = "Read Count Cut-off:10",
                                                      "75"="Read Count Cut-off: 75",
                                                      "200"="Read Count Cut-off:200"))) +
  labs(title = "Number of Detected Genes with each Sampling Method",
       colour = "Sample Type") +
  ylab("# detected genes") +
  xlab("Sample Type") +
  scale_color_discrete(labels = c("Tape Strip", "Biopsy")) +
  theme_bw()
  

ggsave(file = paste0(out_dir, "/generated_figures/01_n_genes_vs_method.png"),
       height = 6,
       width = 9)

saveRDS(fig_n_genes_method, paste0(out_dir, "/generated_rds/01_n_genes_vs_method.rds"))



# Plot of P-value for one-sample t-testing whether the difference in # of genes 
# between method != 0 

fig_PvalDiff <-ggplot(df_pval_vs_cutoff,
       mapping = aes(cutoff, pval_adj)) +
  geom_point(size=0.3) +
  ylim(0, 1) +
  geom_hline(yintercept = 0.05, 
             colour = "black", linetype = "dashed") +
  ylab("fdr adjusted p-value") +
  xlab("Count Cut-off") +
  labs(title = "p-value vs Count Cut-off") +
  theme(plot.title = element_text(size = 10))

fig_PvalDiff

ggsave(filename = paste0(out_dir, "/generated_figures/01_pvalue_diff.png"),
       height = 6,
       width = 10) 

saveRDS(fig_PvalDiff, paste0(out_dir, "/generated_rds/01_PvalVsCutoff.rds"))


# Plot the mean difference in number of genes between methods for each pair by cutoff 

fig_MeanDiff <-ggplot(df_pval_vs_cutoff,
       mapping = aes(cutoff, mean_diff)) +
  geom_point(size=0.3) +
  ylab("Mean Difference") +
  xlab("Count Cut-off") +
  labs(title = "Mean Difference in # of Detected Genes (Biopsies - Tapes) vs \nCount Cut-off") +
  theme(plot.title = element_text(size = 10))

fig_MeanDiff

ggsave(filename = paste0(out_dir, "/generated_figures/01_DiffVsCutoff.png"),
       height = 6,
       width = 10)  

saveRDS(fig_MeanDiff, paste0(out_dir, "/generated_rds/01_DiffVsCutoff.rds"))


# Barplot of number of genes detected by each method and the number of genes in 
# common

common_genes_boxplot<-
  ggplot(df_n_common_genes,
         mapping = aes(x=variable,
                       y=value)) +
  geom_boxplot(aes(fill = variable)) +
  geom_jitter(colour="black",fill="grey", size=1, width = 0.1) +
  facet_wrap(~ cutoff, labeller = labeller(cutoff = c("10"="Read Count Cut-off: 10",
                                                      "75"="Read Count Cut-off: 75",
                                                      "200" = "Read Count Cut-off: 200"))) +
  theme_bw() +
  xlab("Patient Number") +
  ylab("# genes") +
  labs(title = "# Genes Detected by Both Methods",
       fill = "Detected Genes") +
  scale_fill_brewer(palette = "Set3",
                    labels = c("Common", "Tape Strip", "Biopsy"))

common_genes_boxplot  

ggsave(file=paste0(out_dir, "/generated_figures/01_n_common_genesBOX.png"),
       height = 6,
       width = 11)

saveRDS(common_genes_boxplot, paste0(out_dir, "/generated_rds/01_CommonGenesBoxplot.rds"))

