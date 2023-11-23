## -----------------------------------------------------------------------------
##
## Script name: 11_ImprovedDESeqAnalysis.R
##
## Purpose of script: 
## This script performe the DESeq analysis to identify 
## differentially expressed genes between method and between 
## AD and non-AD. 
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
library(ggplot2)


# Functions 

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

EnsemblToSymbol <- function(ensembl_id, gene_list){
  gene_list$external_gene_name[gene_list$ensembl_gene_id == ensembl_id]
}

#-------------------------------------------------------------------------------
# Import Data, Get Coding Genes, and Create "type_AD_ecz" Variable -------------
#-------------------------------------------------------------------------------


# Load Counts and metadata -----------------------------------------------------

meta <- 
  read_tsv("data/Metadata_Final.tsv", 
           col_types = list(Sample_Type = col_factor(levels = c("tape_strip", "biopsy")))
  )

count <- 
  read_csv("results/generated_data/01_df_count_filtered.csv") %>%
  column_to_rownames("gene")


## Create "type_AD_ecz" Variable -----------------------------------------------

meta_fixed <- meta %>% 
  #column_to_rownames("Sample") %>%
  mutate("Eczema" =  factor(ifelse(meta$Eczema_dorsalhand=="Yes", "Yes", "No"))) %>%
  mutate("AD" = factor(AD, levels = c("No", "Yes"))) %>%
  mutate("type_AD_ecz"  = paste0("TYPE", Sample_Type, "AD", AD, "ECZ", Eczema)) %>%
  dplyr::select(Sample, Eczema, Participant, Sample_Type, AD, type_AD_ecz) %>%
  arrange(Sample) %>% 
  column_to_rownames("Sample")
  



## Count matrix - only coding genes --------------------------------------------

coding_genes <- getCodingGenes()

count_coding <- count[row.names(count) %in% coding_genes$ensembl_gene_id,]




## Create deseq2 object --------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count_coding,
                              colData = meta_fixed,
                              design = ~ type_AD_ecz) # Will make new model matrix


dds <- DESeq(dds)


# Generate the Results 

#Tape vs Biopsy - Healthy Controls
res_TapeVsBiopsy_NoADorEcz <- results(dds, contrast = c("type_AD_ecz", "TYPEtape_stripADNoECZNo", "TYPEbiopsyADNoECZNo"), alpha = 0.05) 

#Tape Vs Biopsy - Eczema of AD patients 
res_TapeVsBiopsy_ADandEcz <- results(dds, contrast = c("type_AD_ecz", "TYPEtape_stripADYesECZYes", "TYPEbiopsyADYesECZYes"), alpha = 0.05)

#Tapes Only - AD+Eczema vs Non-AD Control
res_Tape_ControlVsADEcz <- results(dds, contrast = c("type_AD_ecz", "TYPEtape_stripADYesECZYes", "TYPEtape_stripADNoECZNo"), alpha = 0.05)

#Tapes Only - AD+No Eczema vs Non-AD Control
res_Tape_ControlVsADNoEcz <- results(dds, contrast = c("type_AD_ecz", "TYPEtape_stripADYesECZNo", "TYPEtape_stripADNoECZNo"), alpha = 0.05)

#Biops Only - AD+Eczema vs Non-AD Control
res_Biops_ControlVsADEcz <- results(dds, contrast = c("type_AD_ecz", "TYPEbiopsyADYesECZYes", "TYPEbiopsyADNoECZNo"), alpha = 0.05)

#Biops Only - AD+No Eczema vs Non-AD Control
res_Biops_ControlVsADNoEcz <- results(dds, contrast = c("type_AD_ecz", "TYPEbiopsyADYesECZNo", "TYPEbiopsyADNoECZNo"), alpha = 0.05)



# Filter DESeq results 


formatDESeqOut<-function(res){
  data.frame(res) %>%
    filter(padj<0.05, abs(log2FoldChange) > 1) %>%
    rownames_to_column("gene") %>%
    left_join(coding_genes[c("ensembl_gene_id", "external_gene_name")],
              by = c("gene" = "ensembl_gene_id")) %>%
    relocate(external_gene_name, .after = gene)
}

#Filtered Results

res_TapeVsBiopsy_NoADorEcz_sig  <- formatDESeqOut(res_TapeVsBiopsy_NoADorEcz)
res_TapeVsBiopsy_ADandEcz_sig <- formatDESeqOut(res_TapeVsBiopsy_ADandEcz)
res_Tape_ControlVsADEcz_sig <- formatDESeqOut(res_Tape_ControlVsADEcz)
res_Tape_ControlVsADNoEcz_sig <- formatDESeqOut(res_Tape_ControlVsADNoEcz)
res_Biops_ControlVsADEcz_sig <- formatDESeqOut(res_Biops_ControlVsADEcz)
res_Biops_ControlVsADNoEcz_sig <- formatDESeqOut(res_Biops_ControlVsADNoEcz)


# Export CSV

write_csv(res_TapeVsBiopsy_NoADorEcz_sig, file.path(out_dir, "generated_data/11_res_TapeVsBiopsy_NoADorEcz_sig.csv"))
write_csv(res_TapeVsBiopsy_ADandEcz_sig,  file.path(out_dir, "generated_data/11_res_TapeVsBiopsy_ADandEcz_sig.csv"))
write_csv(res_Tape_ControlVsADEcz_sig,  file.path(out_dir, "generated_data/11_res_Tape_ControlVsADEcz_sig.csv"))
write_csv(res_Tape_ControlVsADNoEcz_sig,  file.path(out_dir, "generated_data/11_res_Tape_ControlVsADNoEcz_sig.csv"))
write_csv(res_Biops_ControlVsADEcz_sig,  file.path(out_dir, "generated_data/11_res_Biops_ControlVsADEcz_sig.csv"))
write_csv(res_Biops_ControlVsADNoEcz_sig,  file.path(out_dir, "generated_data/11_res_Biops_ControlVsADNoEcz_sig.csv"))






