## ---------------------------
##
## Script name: 02_global_read_patterns.R
##
## Purpose of script: 
## This script performs the analysis and generates the figures 
## for exploring the gene expression patterns (PCA, etc). 
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
library(biomaRt)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(factoextra)


# Define Functions 

# filter out non-coding genes
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



# Load data --------------------------------------------------------------------

meta <- read_tsv("data/Metadata_Final.tsv", 
                 col_types = list(Sample_Type = col_factor(levels = c("tape_strip", "biopsy")))
)




count_vst <- 
  read_csv("./results/generated_data/01_df_count_vst.csv") %>%
    column_to_rownames("gene")          
 


#Fix sample names so all biops start with B and tapes with T

if(all(colnames(count_vst)==meta$Sample)){
  
  #For visualization, going to create a better Sample column to specify biops and ts 
  
  meta$Sample2 <- 
    paste0(
      toupper(sub('(^.).*', '\\1', meta$Sample_Type)),
      gsub("P", "", meta$Participant))
  
  meta$Sample <- meta$Sample2
  colnames(count_vst) <- meta$Sample
  
} else { 
  stop("ERROR: Something is wrong with sample names in metadata and counts!")
}


# Get coding genes and subset 

coding_genes <- getCodingGenes()

count_vst_coding <- count_vst[rownames(count_vst) %in% coding_genes$ensembl_gene_id,]



# Define colors for plots ------------------------------------------------------
sample_type.cols <- c("steelblue1", "indianred3")
names(sample_type.cols) <- levels(meta$Sample_Type)

AD.cols <- c("forestgreen", "tomato1")
names(AD.cols) <- levels(meta$AD)

eczema.cols <- c("darkorange", "lightgreen", "lightpink1")
names(eczema.cols) <- levels(meta$Eczema_dorsalhand)


## PCA of all samples together Coding genes) -----------------------------------

mat <- t(count_vst_coding)

# Calculate PCA
pca <- prcomp(mat)

res.pca <- 
  data.frame(pca$x) %>%
    rownames_to_column("Sample") %>% 
    left_join(meta, by = "Sample")



# scree plot
scree_all <- fviz_eig(pca,
                      title = "Scree Plot: All Samples",
                      barfill = "darkseagreen",
                      barcolor = "black",
                      gg_theme = theme_minimal())
scree_all
ggsave(scree_all, 
       file = "./results/generated_figures/02_scree_plot.png",
       height = 5,
       width = 8)

saveRDS(scree_all, file.path(out_dir, "generated_rds/02_scree_plot.rds"))
        
# % variance explained by each component
var_explained <- get_eig(pca)$variance.percent 

# PCA Biplot - PC1 and PC2 

biplot_12_all <- ggplot(res.pca,
       aes(x = PC1, y = PC2)) +
  geom_point(aes(color = AD, shape = Sample_Type),
             size = 3.5,
             alpha = .8) +
  xlab(paste("PC1 (", round(var_explained[1], 2), " % var explained)",
             sep = "")) + 
  ylab(paste("PC2 (", round(var_explained[2], 2), " % var explained)",
             sep = "")) +
  labs(title = "PC1 vs PC2 biplot: All samples",
       shape = "Sample Type", color = "AD") +
  scale_shape_discrete(labels = c("Tape Strip", "Biopsy")) +
  scale_color_manual(values = AD.cols) +
  theme_bw()

saveRDS(biplot_12_all, file.path(out_dir, "generated_rds/02_PCA_plot_AllSamples.rds"))

biplot_12_all

ggsave(biplot_12_all,
       file =  file.path(out_dir, "generated_figures/02_PCA_plot_AllSamples.png"),
       height = 5,
       width = 6)


## Distance matrix of all samples (Coding Genes) -------------------------------

dist_all <- dist(t(count_vst_coding),
                 method = "euclidean", 
                 diag = T, 
                 upper = T)



heatmap <- pheatmap(as.matrix(dist_all),
         main = "Euclidean Distance Matrix",
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = column_to_rownames(meta[,c("Sample", "AD", "Sample_Type", "Eczema_dorsalhand")], 
                                             "Sample"),
         annotation_colors = list(AD = c(Yes = "#ff826c", No = "#4ea24e"),
                                  Sample_Type = c(tape_strip = "#ffffb3", biopsy = "#bebada")),
         cluster_rows = T) 
heatmap

saveRDS(heatmap, file.path(out_dir, "generated_rds/02_distance_heatmap.rds"))

ggsave(heatmap,
       file = file.path(out_dir, "generated_figures/02_heatmap_sub.png"),
       height = 6,
       width = 7)


