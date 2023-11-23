## ---------------------------
##
## Script name: 06_multiQC.R
##
## Purpose of script: 
## This script parses and analyzes sopem of the results from MultiQC,  
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
library(readr)
library(tidyverse)
library(ggplot2)


#-------------------------------------------------------------------------------
# Import the data --------------------------------------------------------------
#-------------------------------------------------------------------------------


# Load data --------------------------------------------------------------------
stats <- read_tsv("data/multiqc_general_stats.txt")

featureCounts <- read_tsv("data/multiqc_featureCounts.txt")

# Import the metadata 

meta <- read_tsv("data/Metadata_Final.tsv", 
                         col_types = list(Sample_Type = col_factor(levels = c("tape_strip", "biopsy")))
)


#-------------------------------------------------------------------------------
# Analysis ---------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Data wrangling ---------------------------------------------------------------
stats <- meta %>% 
  left_join(stats, by=c("Sample"))

featureCounts <- 
  meta %>% 
    left_join(featureCounts, by=c("Sample")) %>%
    dplyr::select(c(Participant, Sample, Sample_Type, Subject_type, Total, Assigned, Unassigned_Unmapped,
             Unassigned_NoFeatures, Unassigned_Overlapping_Length, percent_assigned))

write_csv(stats, file.path(out_dir, "generated_data/06_seq_stat_summary.csv"))
write_csv(featureCounts, file.path(out_dir, "generated_data/06_featureCounts_summary.csv"))




#For visualization, going to create a better Sample column to specify biops and ts 
  
  featureCounts$Sample <- 
    paste0(
      toupper(sub('(^.).*', '\\1', featureCounts$Sample_Type)),
      gsub("P", "", featureCounts$Participant))


#-------------------------------------------------------------------------------
# Plots ------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# Define colors for plots ------------------------------------------------------
sample_type.cols <- c("steelblue1", "indianred3")
names(sample_type.cols) <- levels(meta$Sample_Type)

AD.cols <- c("tomato1", "forestgreen")
names(AD.cols) <- levels(meta$AD)

eczema.cols <- c("darkorange", "lightgreen", "lightpink1")
names(eczema.cols) <- levels(meta$Eczema_dorsalhand)



## Initial number of sequences--------------------------------------------------

initial_reads <- stats %>% 
  ggplot(mapping = aes(x=Sample_Type,
                       y=`FastQC_mqc-generalstats-fastqc-total_sequences`,
                       fill=Sample_Type)) +
  geom_boxplot(alpha=.7) +
  #geom_point(size=2, alpha=0.5) +
  #geom_jitter(colour="black", size=1) +
  labs(title="Initial # of Read Pairs",
       fill="Sample Type") +
  xlab(NULL)+
  ylab("# Read Pairs") +
  scale_fill_manual(labels = c("Tape Strip", "Biopsy"),
                    values = sample_type.cols) +
  theme_bw()

ggsave(initial_reads,
       file= file.path(out_dir, "generated_figures/06_init_reads.png"),
       height = 5,
       width=7)

# Feature Count Statistics plot 

featureCounts_melt <- 
  featureCounts %>%
  dplyr::select(-c(Total, percent_assigned)) %>%
  melt(id.vars = c("Participant", "Sample", "Sample_Type", "Subject_type"))

assigned <- featureCounts_melt %>%
  ggplot(mapping = aes(x=Sample,
                       y=value/1000000,
                       fill=variable)) +
  geom_bar(stat = "identity",
           position = "stack") +
  facet_wrap(~ Sample_Type,
             scales = "free_x",
             labeller = labeller(Sample_Type = c("tape_strip" = "Sampling Method: Tape Strip",
                                            "biopsy" = "Sampling Method: Biopsy"))) +
  ylab("# Read Pairs (M)") +
  xlab("Patient nr") +
  labs(title = "featureCounts: Read Assignments",
       fill = NULL) +
  theme_bw()

assigned

ggsave(assigned,
       file = file.path(out_dir, "generated_figures/06_assigned_reads.png"),
       height = 5,
       width=9)


saveRDS(assigned, file.path(out_dir, "generated_rds/06_featureCounts_figure.rds"))
