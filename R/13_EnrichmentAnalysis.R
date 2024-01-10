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

# Individual Plots -------------------------------------------------------------

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
    mutate("baseMean" = signif(baseMean, 3)) %>%
    arrange(desc(log2FoldChange))
  
  pal <- function(x) {
    f_neg <- scales::col_numeric(
      palette = c(color_down,'white' ), 
      domain= c(0,-15),
      na.color = "black")
    
    f_pos <- scales::col_numeric(
      palette = c('white', color_up),
    domain = c(0,15))

ifelse(x < 0, f_neg(x), f_pos(x))
  }

  out <- 
  gt(df) %>%
    data_color(columns = log2FoldChange, fn = pal) %>%
    cols_label(
      external_gene_name = "Gene",
      padj = "Adjusted p-value"
    ) %>% 
    cols_align(align = "center")
}




# Tapes Vs. Biopsies -----------------------------------------------------------

# Inflammation - AD + Eczema 


makeHighlightTable(data[["TapeVsBiopsy_ADandEcz"]], immune_response$gene_symbol,
                   color_up = "red", color_down = "blue") %>%
  tab_header(title = md("Immune Reponse Genes By Method (Active AD)"),
             subtitle = md("<span style='color:#f48f73'> Tape-strip </span>
                           vs 
                           <span style='color:#9067fc'> Biopsy </span>")) %>%
  tab_options(table.font.size = 75) %>%
  gtsave(paste0(out_dir, "/generated_figures/13_tbl_Immune_ADandEcz.png"), vwidth = 2000)


# Inflammation - Non-AD 


makeHighlightTable(data[["TapeVsBiopsy_NoADorEcz"]], immune_response$gene_symbol,
                   color_up = "red", color_down = "blue") %>%
  tab_header(title = md("Immune Reponse Genes By Method (Non-AD)"),
             subtitle = md("<span style='color:#f48f73'> Tape-strip </span>
                           vs 
                           <span style='color:#9067fc'> Biopsy </span>")) %>%
  tab_options(table.font.size = 75) %>%
  gtsave(paste0(out_dir, "/generated_figures/13_tbl_Immune_NoADorEcz.png"), vwidth = 2000)


# Skin Barrier - AD + Eczema 

makeHighlightTable(data[["TapeVsBiopsy_ADandEcz"]], skin_barrier$gene_symbol,
                   color_up = "red", color_down = "blue") %>%
  tab_header(title = md("Skin Barrier Genes By Method (Active AD)"),
             subtitle = md("<span style='color:#f48f73'> Tape-strip </span>
                           vs 
                           <span style='color:#9067fc'> Biopsy </span>")) %>%
  tab_options(table.font.size = 75) %>%
  gtsave(paste0(out_dir, "/generated_figures/13_tbl_Barrier_ADandEcz.png"), vwidth = 2000)


# Skin Barrier - Non-AD   


makeHighlightTable(data[["TapeVsBiopsy_NoADorEcz"]], skin_barrier$gene_symbol,
                   color_up = "red", color_down = "blue") %>%
  tab_header(title = md("Skin Barrier Genes By Method (Non-AD)"),
             subtitle = md("<span style='color:#f48f73'> Tape-strip </span>
                           vs 
                           <span style='color:#9067fc'> Biopsy </span>")) %>%
  tab_options(table.font.size = 75) %>%
  gtsave(paste0(out_dir, "/generated_figures/13_tbl_Barrier_NoADorEcz.png"), vwidth = 2000)



# Itch - AD + Eczema 

makeHighlightTable(data[["TapeVsBiopsy_ADandEcz"]], itch$gene_symbol,
                   color_up = "red", color_down = "blue") %>%
  tab_header(title = md("Itch Genes By Method (Active AD)"),
             subtitle = md("<span style='color:#f48f73'> Tape-strip </span>
                           vs 
                           <span style='color:#9067fc'> Biopsy </span>")) %>%
  tab_options(table.font.size = 75)  %>%
  gtsave(paste0(out_dir, "/generated_figures/13_tbl_Itch_ADandEcz.png"), vwidth = 2000)



# Itch  - Non-AD 


makeHighlightTable(data[["TapeVsBiopsy_NoADorEcz"]], itch$gene_symbol,
                   color_up = "red", color_down = "blue") %>%
  tab_header(title = md("Itch Genes By Method (Healthy Non-AD)"),
             subtitle = md("<span style='color:#f48f73'> Tape-strip </span>
                           vs 
                           <span style='color:#9067fc'> Biopsy </span>")) %>%
  tab_options(table.font.size = 75) %>%
  gtsave(paste0(out_dir, "/generated_figures/13_tbl_Itch_NoADorEcz.png"), vwidth = 2000)

# Make Combined Plots ----------------------------------------------------------

plotPathways <- function(pathway){
  plt_adecz<-
  ggdraw() + draw_image(paste0(out_dir, "/generated_figures/13_tbl_",pathway,"_ADandEcz.png"), valign = 1)+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
  
  plt_NOadecz<-
    ggdraw() + draw_image(paste0(out_dir, "/generated_figures/13_tbl_",pathway,"_NoADorEcz.png"), valign = 1)+
      theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

  plot_grid(plt_adecz, plt_NOadecz)
}

pw_plots_Immune <- plotPathways("Immune")
pw_plots_Barrier <- plotPathways("Barrier")
pw_plots_Itch <- plotPathways("Itch")

ggsave(paste0(out_dir, "/generated_figures/13_ImmuneAll.png"), pw_plots_Immune)
ggsave(paste0(out_dir, "/generated_figures/13_BarrierAll.png"), pw_plots_Barrier)
ggsave(paste0(out_dir, "/generated_figures/13_Itch.png"), pw_plots_Itch)


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
