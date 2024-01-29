## ---------------------------
##
## Script name: 08_metadataTable.R
##
## Purpose of script: 
## This script summarized the metadata and then makes a
## a gtable summary and saves is as an RDS for use in the 
## manuscript
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
library(gtsummary)

# Load data --------------------------------------------------------------------

meta <- 
  read_tsv("data/Metadata_Final.tsv", 
                         col_types = list(Sample_Type = col_factor(levels = c("tape_strip", "biopsy")),
                                          Sex = col_factor(levels = c("F", "M")),
                                          Fillagrin_mutation = col_factor(levels = c("No", "Homozygous", "Heterozygous")),
                                          Eczema_dorsalhand = col_factor(levels = c("No", "Yes")))) %>% 
  select(-c(Sample_Type, Sample)) %>%
  distinct()




# Clean up metadata for output stats table -------------------------------------

out_table <- 
  meta %>%
    tbl_summary(
      by = AD,
      include = c(Age_years, Sex, BMI, Fillagrin_mutation, Eczema_dorsalhand, EASI_total, TLSS_total),
      label = list(Age_years ~ "Age (years)", 
                   BMI ~ "BMI", 
                   Sex ~ "Sex",
                   Fillagrin_mutation ~ "Fillagrin Gene Mutation",
                   EASI_total ~ "Eczema Area and Severity Index (EASI)",
                   TLSS_total ~ "Total Lesion Severity Score (TLSS)\n(Dorsal Hand)",
                   Eczema_dorsalhand ~ "Eczema (Dorsal Hand)"),
      statistic = list(all_continuous() ~ "{mean} (± {sd})"),
      digits = list(EASI_total ~ c(0),
                    TLSS_total ~ c(0)),
      type = list(EASI_total ~ "continuous",
                  TLSS_total ~ "continuous",
                  Eczema_dorsalhand ~ "categorical"),
      missing = "no") %>%
     modify_spanning_header(c("stat_1", "stat_2") ~ "Atopic Dermatitis") %>%
     modify_table_styling(
       columns = label,
       rows = label == "TLSS (Dorsal Hand)",
       footnote = "Only for AD participants with eczema on dorsal hand"
     ) %>%
  modify_caption("Summary of Participant Metadata") %>%
  modify_footnote(all_stat_cols() ~ "Mean (± SD) for Age, BMI, EASI, TLSS; n(%) for Fillagrin Gene Mutation, Eczema on dorsal hand")

out_table$table_body$stat_1[out_table$table_body$stat_1 == "NA (± NA)"] <- "NA"
out_table$table_body$stat_2[out_table$table_body$stat_2 == "NA (± NA)"] <- "NA"

out_table


saveRDS(out_table, file.path(out_dir, "generated_rds/08_summary_table.rds"))





