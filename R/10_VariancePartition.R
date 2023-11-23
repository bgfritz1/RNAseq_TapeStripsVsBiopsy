## ---------------------------
##
## Script name: 10_VariancePartition.R
##
## Purpose of script: 
## This script utilizes the "variancepartition" R package 
## to partition the variance and assess the contribution of 
## different metadata factors to the overall variance in 
## the data. 
## 
##
## Author: Blaine Fritz
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
library(variancePartition)

# Functions --------------------------------------------------------------------

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
# Import -----------------------------------------------------------------------
#-------------------------------------------------------------------------------


# Load data --------------------------------------------------------------------
meta <-  read_tsv("data/Metadata_Final.tsv", 
                         col_types = list(Sample_Type = col_factor(levels = c("tape_strip", "biopsy")))
)

count_vst <- read_csv("./results/generated_data/01_df_count_vst.csv")

coding_genes <- getCodingGenes()


#-------------------------------------------------------------------------------
# Analysis---------------------------------------------------------------------- 
#-------------------------------------------------------------------------------

# Only use coding genes 

count_vst_coding <- 
  count_vst %>% 
    filter(gene %in% coding_genes$ensembl_gene_id) %>%
    column_to_rownames("gene") 

# Get samples of each 
tapes <- meta$Sample[meta$Sample_Type == "tape_strip"]
biops <- meta$Sample[meta$Sample_Type == "biopsy"]


# Subset counts 
count_vst_coding_tape <- count_vst_coding[,tapes]
count_vst_coding_biops <- count_vst_coding[,biops] 


# there are some genes with zero variance for the biops data. Need to remove these (its 42 genes)

idx <- which(apply(count_vst_coding_biops[,-1], 1, function(x) var(x) == 0))

count_vst_coding_biops <- count_vst_coding_biops[-c(idx),]


# Run Variance Partition -------------------------------------------------------

# All data 

form <- ~ (1|Sample_Type) + (1|Eczema_dorsalhand) + (1|AD)
form2 <- ~ (1|AD) + (1|Eczema_dorsalhand)

varPart <- fitExtractVarPartModel(count_vst_coding, form, meta)

vp <- sortCols(varPart)
vp_plot_all <- plotVarPart(vp)

saveRDS(varPart, file.path(out_dir, "generated_rds/10_fitExtractVarPartModel_AllSamples.RDS"))
saveRDS(vp_plot_all, file.path(out_dir, "generated_rds/10_vp_plot_all.RDS"))

# Only Tapes

varPart_tapes <- fitExtractVarPartModel(count_vst_coding_tape, form2, meta[meta$Sample %in% tapes,])
vp_tapes <- sortCols(varPart_tapes)
vp_plot_tape <- plotVarPart(vp_tapes)

saveRDS(varPart_tapes, file.path(out_dir, "generated_rds/10_fitExtractVarPartModel_Tapes.RDS"))
saveRDS(vp_plot_tape, file.path(out_dir, "generated_rds/10_vp_plot_tape.RDS"))

# Only Biopsies

varPart_biops <- fitExtractVarPartModel(count_vst_coding_biops, form2, meta[meta$Sample %in% biops,])
vp_biops <- sortCols(varPart_biops)
vp_plot_biops <- plotVarPart(vp_biops)

saveRDS(varPart_biops, file.path(out_dir, "generated_rds/10_fitExtractVarPartModel_Biops.RDS"))
saveRDS(vp_plot_biops, file.path(out_dir, "generated_rds/10_vp_plot_biops.RDS"))






