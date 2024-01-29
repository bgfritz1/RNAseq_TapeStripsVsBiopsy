## -----------------------------------------------------------------------------
##
## Script name: 13_EnrichmentAnalysis.R
##
## Purpose of script: 
## This script performs the testing for whether the DEGs identified 
## in 11_ImprovedDESeqAnalysis.R represent enrichment for groups 
## of genes associated with skin battier, itch, immune response, and different
## cell types 
## 
##
## Author: Blaine Fritz 
##
## Email: bgfritz@sund.ku.dk
##------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Setup ------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(gt)
library(cowplot)
library(gtable)
library(grid)
library(gridtext)
library(gridExtra)
library(ggtext)
library(ggpubr)


# Functions --------------------------------------------------------------------

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

# Define enrichment function 

calcEnrichment <- function(de_genes, n_genes, gene_set, dirs = NULL){
  overlap_genes <- de_genes[de_genes %in% gene_set]
  n_overlap <- sum(de_genes %in% gene_set)
  n_total <- n_genes
  n_pathway <- length(gene_set)
  n_DE <- length(de_genes)
  
  
  
  res <- phyper(q = n_overlap - 1,
                m = n_pathway,
                n = n_total-n_pathway,
                k = length(de_genes),
                lower.tail = F)
  
  if(length(dirs)>0){
    overlap_genes <- paste0(de_genes[de_genes %in% gene_set],
                            dirs[de_genes %in% gene_set])
  }
  
  
  return(list("q" =  n_overlap - 1, 
              "m" =  n_pathway, 
              "n" =  n_total-n_pathway,
              "k" = n_DE,
              "p_val" = res, 
              "n_DE" = n_DE, 
              "n_overlap" = n_overlap, 
              "n_gene_set" = n_pathway, 
              "overlap_genes" = toString(unlist(overlap_genes))))
}


# Check that output directories exist and, if not, create them -----------------

out_dir <- "./results"

if (!dir.exists(out_dir)) {dir.create(out_dir)}

dirs <- c("generated_data", 
          "generated_figures",
          "generated_rds") 

for(dir in file.path(out_dir, dirs)){
  if (!dir.exists(dir)) {dir.create(dir)}
}

#-------------------------------------------------------------------------------
# Import Pathways and Data -----------------------------------------------------
#-------------------------------------------------------------------------------

# Import Path Genes ------------------------------------------------------------

skin_barrier <- read.table("./data/BarrierGenes_DelDucaEtAl.txt", header = T)
immune_response <- read.table("data/Immune_Response.txt", header = T)
itch <- read.table("./data/ItchGenes_DelDucaEtAl.txt", header = T)

#Combined Data Set 

pathways <- list(skin_barrier, immune_response, itch)
names(pathways) <- c("skin_barrier", "immune_response", "itch")

pathways <- bind_rows(pathways, .id = "gene_set")

# Read in the DESeq data -------------------------------------------------------

categories <- c(
  "TapeVsBiopsy_NoADorEcz",
  "TapeVsBiopsy_ADandEcz",
  "Tape_nonADVsADEcz",
  "Tape_nonADVsADNoEcz",
  "Biops_nonADVsADEcz",
  "Biops_nonADVsADNoEcz")

data <- lapply(paste0(out_dir, "/generated_data/11_res_", categories, "_sig.csv"), function(x){

  df <- 
    read_delim(x) %>%
      mutate("dir" = ifelse(log2FoldChange>0, "(+)", "(-)"))
  
  return(df)
})
               


names(data) <- categories


#Import Coding Genes

coding_genes <- getCodingGenes()


#Check DelDuca Genes in my data  - A few missing 
length(intersect(skin_barrier$gene_symbol, coding_genes$external_gene_name))
length(intersect(itch$gene_symbol, coding_genes$external_gene_name))
length(intersect(immune_response$gene_symbol, coding_genes$external_gene_name))


# A few missing, but not too many. Will continue

# Enrichment Analysis ----------------------------------------------------------

# Tapes - AD Eczema vs Non-AD

counts <- 17975 #Number of genes used in DEseq data 

enrich_tapes_ADEczvsnonAD <- 
     lapply(unique(pathways$gene_set), function(x){
       enrich <- calcEnrichment(data[["Tape_nonADVsADEcz"]]$external_gene_name, 
                                counts, 
                                pathways$gene_symbol[pathways$gene_set == x],
                                dir = data[["Tape_nonADVsADEcz"]]$dir)
       enrich <- unlist(enrich)
       enrich$name <- x 
       return(unlist(enrich))
     }) %>% bind_rows()


enrich_tapes_ADNoEczvsnonAD  <- 
  lapply(unique(pathways$gene_set), function(x){
    enrich <- calcEnrichment(data[["Tape_nonADVsADNoEcz"]]$external_gene_name, 
                             counts, 
                             pathways$gene_symbol[pathways$gene_set == x],
                             dir = data[["Tape_nonADVsADNoEcz"]]$dir)
    enrich <- unlist(enrich)
    enrich$name <- x 
    return(unlist(enrich))
  }) %>% bind_rows()

enrich_biops_ADEczvsnonAD  <- 
  lapply(unique(pathways$gene_set), function(x){
    enrich <- calcEnrichment(data[["Biops_nonADVsADEcz"]]$external_gene_name, 
                             counts, 
                             pathways$gene_symbol[pathways$gene_set == x],
                             dir = data[["Biops_nonADVsADEcz"]]$dir)
    enrich <- unlist(enrich)
    enrich$name <- x 
    return(unlist(enrich))
  }) %>% bind_rows()

enrich_biops_ADNoEczvsnonAD  <- 
  lapply(unique(pathways$gene_set), function(x){
    enrich <- calcEnrichment(data[["Biops_nonADVsADNoEcz"]]$external_gene_name, 
                             counts, 
                             pathways$gene_symbol[pathways$gene_set == x],
                             dir = data[["Biops_CnonADVsADNoEcz"]]$dir)
    enrich <- unlist(enrich)
    enrich$name <- x 
    return(unlist(enrich))
  }) %>% bind_rows()


enrich_TapesVSBiops_ADEcz  <- 
  lapply(unique(pathways$gene_set), function(x){
    enrich <- calcEnrichment(data[["TapeVsBiopsy_ADandEcz"]]$external_gene_name, 
                             counts, 
                             pathways$gene_symbol[pathways$gene_set == x],
                             dir = data[["TapeVsBiopsy_ADandEcz"]]$dir)
    enrich <- unlist(enrich)
    enrich$name <- x 
    return(unlist(enrich))
  }) %>% bind_rows()


enrich_TapesVSBiops_NoADorEcz  <- 
  lapply(unique(pathways$gene_set), function(x){
    enrich <- calcEnrichment(data[["TapeVsBiopsy_NoADorEcz"]]$external_gene_name, 
                             counts, 
                             pathways$gene_symbol[pathways$gene_set == x],
                             dir = data[["TapeVsBiopsy_NoADorEcz"]]$dir)
    enrich <- unlist(enrich)
    enrich$name <- x 
    return(unlist(enrich))
  }) %>% bind_rows()


# For Tapes Vs Biopsies - Only considering overlapping genes--------------------

up_tapes_ecz <- data[["TapeVsBiopsy_ADandEcz"]]$external_gene_name[data[["TapeVsBiopsy_ADandEcz"]]$log2FoldChange>1]
up_tapes_NoEcz <- data[["TapeVsBiopsy_NoADorEcz"]]$external_gene_name[data[["TapeVsBiopsy_NoADorEcz"]]$log2FoldChange>1]

up_biops_ecz <- data[["TapeVsBiopsy_ADandEcz"]]$external_gene_name[data[["TapeVsBiopsy_ADandEcz"]]$log2FoldChange<1]
up_biops_NoEcz <- data[["TapeVsBiopsy_NoADorEcz"]]$external_gene_name[data[["TapeVsBiopsy_NoADorEcz"]]$log2FoldChange<1]

intersect(up_tapes_ecz, up_tapes_NoEcz)


# Visualize --------------------------------------------------------------------

#Format input data for highlight table

makeHighlightTable <- function(de_data, pathway_genes, color_up, color_down, n_cutoff=20){
  
  overlapping_genes <- intersect(de_data$external_gene_name, pathway_genes)
  
  if(length(overlapping_genes)==0){
    warning("NO OVERLAPPING GENES")
    return(gt(data.frame("Error" = "No Overlapping Genes")))
  }
  
  top_positive <- 
    de_data %>% 
    filter(external_gene_name %in% overlapping_genes) %>%
    filter(log2FoldChange>0) %>%
    slice_min(padj, n = 10)
  
  top_negative <- 
    de_data %>% 
    filter(external_gene_name %in% overlapping_genes) %>%
    filter(log2FoldChange<0) %>%
    slice_min(padj, n = 10)  
  
  df <- 
    bind_rows(top_positive, top_negative) %>% 
    dplyr::select(external_gene_name, baseMean, log2FoldChange, padj) %>%
    mutate("padj" = signif(padj, digits = 3)) %>%
    mutate("log2FoldChange" = signif(log2FoldChange, digits = 3)) %>%
    mutate("baseMean" = signif(baseMean, 0)) %>%
    arrange(desc(log2FoldChange))
  
  colnames(df) <- c("Gene", "Mean\nExpression", "Log2FC", "Adjusted\nP-value")
  return(df)
}


# Make the Tables 

# Tapes Vs. Biopsies -----------------------------------------------------------

# Inflammation - AD + Eczema 

df_immune_AD<-
  makeHighlightTable(data[["TapeVsBiopsy_ADandEcz"]], immune_response$gene_symbol)

# Inflammation - Non-AD 

df_immune_nonAD<-
  makeHighlightTable(data[["TapeVsBiopsy_NoADorEcz"]], immune_response$gene_symbol)
# Skin Barrier - AD + Eczema 
df_barrier_AD<-
  makeHighlightTable(data[["TapeVsBiopsy_ADandEcz"]], skin_barrier$gene_symbol)

# Skin Barrier - Non-AD   

df_barrier_nonAD<-
  makeHighlightTable(data[["TapeVsBiopsy_NoADorEcz"]], skin_barrier$gene_symbol)

# Itch - AD + Eczema 
df_itch_AD<-
  makeHighlightTable(data[["TapeVsBiopsy_ADandEcz"]], itch$gene_symbol)

# Itch  - Non-AD 

df_itch_nonAD<-
  makeHighlightTable(data[["TapeVsBiopsy_NoADorEcz"]], itch$gene_symbol)


# Make Figure 3 ----------------------------------------------------------


#Import Venns 

venn_ADEczVsnonAD <- read_rds("./results/generated_rds/12_venn_ADEczVsnonAD.rds")
venn_ADNoEczVsnonAD <- read_rds("./results/generated_rds/12_venn_ADNoEczVsnonAD.rds")
Venn_TapesVsBiopsy_NoAD <- read_rds("./results/generated_rds/12_Venn_TapesVsBiopsy_NoAD.rds")
Venn_SharedDEgenes_tapes <- read_rds("./results/generated_rds/12_Venn_SharedDEgenes_tapes.rds")


#Justification function bc not built into tableGrob

justify <- function(x, just = 1) {
  height <- as.numeric(grid::convertHeight(sum(x$heights), 'npc'))
  x$vp <- grid::viewport(y = just + height * (0.5 - just))
  x
}

#coloring function
pal <- function(x, color_up, color_down) {
  f_neg <- scales::col_numeric(
    palette = c(color_down,'white' ), 
    domain= c(0,-15),
    na.color = "black")
  
  f_pos <- scales::col_numeric(
    palette = c('white', color_up),
    domain = c(0,15))
  
  ifelse(x < 0, f_neg(x), f_pos(x))
}

# Format Grobs 

formatTables <- function(x, title){
  
  #Generate the background colors for the LFC
  bg_cols <- pal(x$Log2FC, color_up = "red", color_down="#7065fb")
  
  bg_mat <- matrix("white", nrow(x), ncol(x))
  bg_mat[,3] <- bg_cols
  
  theme <- ttheme_default(base_size = 12,
                          core = list(
                          bg_params = list(fill=bg_mat)
                          ))
  
  t1 <- tableGrob(x, rows = NULL, theme = theme) 
  title <- richtext_grob(title, gp = gpar(fontsize=20))
  padding <- unit(5,"mm")
  
  table <- gtable_add_rows(
    t1, 
    heights = grobHeight(title) + padding,
    pos = 0)
  
  table <- gtable_add_grob(
    table, 
    title, 
    1, 1, 1, ncol(table))
  
  return(table)
  
}


# Titles 

fig3B_title1 <- c("Immune Response: <br>
<span style='color:#f48f73'>Tape-Strips</span> 
vs. 
<span style='color:#9067fc'>Biopsy</span> <br>(Active AD)")

fig3B_title2 <- c("Immune Response: <br>
<span style='color:#f48f73'>Tape-Strips</span> 
vs. 
<span style='color:#9067fc'>Biopsy</span> <br>(No AD)")

fig3C_title1 <- c("Itch:<br>
<span style='color:#f48f73'>Tape-Strips</span> 
vs. 
<span style='color:#9067fc'>Biopsy</span> <br>(Active AD)")

fig3C_title2 <- c("Itch:<br>
<span style='color:#f48f73'>Tape-Strips</span> 
vs. 
<span style='color:#9067fc'>Biopsy</span> <br>(No AD)")

fig3D_title1 <- c("Skin Barrier: <br>
<span style='color:#f48f73'>Tape-Strips</span> 
vs. 
<span style='color:#9067fc'>Biopsy</span> <br>(Active AD)")

fig3D_title2 <- c("Skin Barrier: <br>
<span style='color:#f48f73'>Tape-Strips</span> 
vs. 
<span style='color:#9067fc'>Biopsy</span> <br>(No AD)")

Fig3B_1 <-formatTables(df_immune_AD, title = fig3B_title1)
Fig3B_2 <- formatTables(df_immune_nonAD, title = fig3B_title2)
Fig3B <- arrangeGrob(justify(Fig3B_1, just = 0.5), justify(Fig3B_2, just = 0.58), ncol = 2)

Fig3C_1 <- formatTables(df_itch_AD, title = fig3C_title1)
Fig3C_2 <- formatTables(df_itch_nonAD, title = fig3C_title2)
Fig3C <- arrangeGrob(justify(Fig3C_1, just = 1.03), justify(Fig3C_2, just = 1), ncol = 2)

Fig3D_1 <- formatTables(df_barrier_AD, title = fig3D_title1)
Fig3D_2 <- formatTables(df_barrier_nonAD, title = fig3D_title2)
Fig3D <- arrangeGrob(justify(Fig3D_1, just = 0), justify(Fig3D_2, just = 0.03), ncol = 2)

g <- list(Venn_TapesVsBiopsy_NoAD + theme(plot.margin = margin(0.5,2,0.5,2, "cm")), Fig3B, Fig3C, Fig3D)

figure_mat <- rbind(c(1,1,2,2), c(3,3,4,4),c(3,3,4,4))

Fig3<-
  arrangeGrob(grobs = g, layout_matrix = figure_mat) |>
  as_ggplot() +
  draw_plot_label(label = c("A.", "B.", "C.", "D."), 
                  x = c(0, 0.48, 0 , 0.48),
                  y = c(1, 1, 0.62, 0.62),
                  size = 25)


ggsave("./results/generated_figures/Figure3.png", Fig3, width = 15, height = 15)


# Cell Type Enrichent Analysis 

#Get immune cell signatures
#From: https://github.com/ba306/immune-cell-signature-discovery-classification-paper/blob/main/gene_lists_papers/all_immune_cell_signatures_curated.txt

signatures <- 
  read_delim("./data/all_immune_cell_signatures.tsv") %>% 
  dplyr::select(Cell_type, Gene) %>% distinct()


enrich_CellTypes_TapesVsBiopsAD  <- 
  lapply(unique(signatures$Cell_type), function(x){
    enrich <- calcEnrichment(data[["TapeVsBiopsy_ADandEcz"]]$external_gene_name, 
                             counts, 
                             signatures$Gene[signatures$Cell_type == x])
    enrich <- unlist(enrich)
    enrich$name <- x 
    return(unlist(enrich))
  }) %>% bind_rows() %>% mutate("padj" = p.adjust(p_val, method = "fdr")) %>% filter(padj<0.05)


enrich_CellTypes_TapesVsBiopsnoAD  <- 
  lapply(unique(signatures$Cell_type), function(x){
    enrich <- calcEnrichment(data[["TapeVsBiopsy_NoADorEcz"]]$external_gene_name, 
                             counts, 
                             signatures$Gene[signatures$Cell_type == x])
    enrich <- unlist(enrich)
    enrich$name <- x 
    return(unlist(enrich))
  }) %>% bind_rows() %>% mutate("padj" = p.adjust(p_val, method = "fdr")) %>% filter(padj<0.05)


sig_out <- 
  bind_rows(enrich_CellTypes_TapesVsBiopsAD, enrich_CellTypes_TapesVsBiopsnoAD, .id="type") %>% 
  mutate("type" = ifelse(type == 1, "TapeVsBiopsy_ADandEcz", "TapeVsBiopsy_NoADorEcz"))
  
write_delim(sig_out, paste0(out_dir, "/generated_data/13_CellTypes.tsv"), delim = "\t")
