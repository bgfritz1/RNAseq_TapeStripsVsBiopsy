## ---------------------------
##
## Script name: 03_bacterial_reads.R
##
## Purpose of script: 
## This script parses the results from the Kraken output files 
## to assess the relative proportion of bacterial reads present
## in the samples. 
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
# Import and Analyze------------------------------------------------------------
#-------------------------------------------------------------------------------

# Import and Format the Kraken Data

# Load kraken data -------------------------------------------------------------
kraken_dir = "data/kraken_output"

kraken_data <-
  list.files(kraken_dir, full.names = T) %>%
  set_names(., nm = map(.x = ., ~gsub(".kraken.out", "", basename(.x)))) %>%
  map(function(x) {
    x <- read.delim(x, header = F)
    x$V6 <- trimws(x$V6, which = "left")
    x <- subset(x, V6 %in% c("root", "Bacteria", "Homo sapiens"))
    return(x)}) %>%
  map_df(~as.data.frame(.x), .id = "Sample_ID") %>%
  dplyr::select(Sample_ID, "Percent_assigned" = V1, "Total_assigned" = V2, "Organism" = V6) %>%
  group_by(Sample_ID) %>%
  summarise(
    "Total_reads" = Total_assigned[Organism=="root"],
    "Bac_reads" = Total_assigned[Organism=="Bacteria"],
    "Human_reads" = Total_assigned[Organism=="Homo sapiens"],
    "Bac_prcnt" = Percent_assigned[Organism=="Bacteria"],
    "Human_prcnt" = Percent_assigned[Organism=="Homo sapiens"])


# Import the metadata 

meta <- meta <- read_tsv("data/Metadata_Final.tsv", 
                         col_types = list(Sample_Type = col_factor(levels = c("tape_strip", "biopsy")))
)

#Select only Bacterial species and exclude the other domains 

kraken_bac <- 
  list.files(kraken_dir, full.names = T) %>%
  set_names(., nm = map(.x = ., ~gsub(".kraken.out", "", basename(.x)))) %>%
  map(function(x) {
    #### Read in the files and format a bit 
    x <- read.delim(x, header = F)
    x$V6 <- trimws(x$V6, which = "left")
    x <- x[, c(1,2,4,6)] #"Subset columns of interest"
    colnames(x)<-c("Prcnt_root_taxon", "nr_reads_root_taxon", "Rank_code", "ID")
    
    #### We need to extract the bacteria info, but the order of the domains in the file are not consistent
    #### So here is some ugly code to get only bacterial species
    
    domain_indices <- subset(x, Rank_code == "D") #get indexes of domain ranks 
    domain_indices$index <-c(1:nrow(domain_indices)) #make a vector for indexing
    
    
    bac_index <- domain_indices$index[domain_indices$ID =="Bacteria"] #Find which row is bacteria 
    next_domain <- bac_index+1
    
    bac_index <- as.numeric(row.names(domain_indices)[bac_index])
    
    if(next_domain == nrow(domain_indices)){
      stop_index <- nrow(x)
    }
    else{ 
      stop_index <-as.numeric(row.names(domain_indices)[next_domain])-1
    }
    x <- x[bac_index:stop_index,]
    x <- subset(x, Rank_code %in% c("S"))
    return(x)}) %>%
  bind_rows(., .id = "Sample_ID")


# Data wrangling ---------------------------------------------------------------
kraken_data <- 
  meta %>%
  left_join(kraken_data, by = c("Sample" = "Sample_ID"))

#Fix sample names so all biops start with B and tapes with T
  kraken_data$Sample <- 
    paste0(
      toupper(sub('(^.).*', '\\1', kraken_data$Sample_Type)),
      gsub("P", "", kraken_data$Participant))
  
  

# Calculate Abundance ----------------------------------------------------------

kraken_bac_abund <-
  kraken_bac %>%
  group_by(Sample_ID) %>%
  mutate("Rel_abund" = nr_reads_root_taxon/sum(nr_reads_root_taxon)*100) %>%
  left_join(meta, by = c("Sample_ID" = "Sample"))



#-------------------------------------------------------------------------------
# Make the Plots ---------------------------------------------------------------
#-------------------------------------------------------------------------------


# Define colors for plots ------------------------------------------------------
sample_type.cols <- c("steelblue1", "indianred3")
names(sample_type.cols) <- levels(meta$Sample_Type)

AD.cols <- c("tomato1", "forestgreen")
names(AD.cols) <- levels(meta$AD)

eczema.cols <- c("darkorange", "lightgreen", "lightpink1")
names(eczema.cols) <- levels(meta$Eczema_dorsalhand)


# Plots-------------------------------------------------------------------------

bac_percent <- kraken_data %>%
  ggplot(mapping=aes(x=Sample,
                     y=Bac_prcnt, 
                     fill=Sample_Type)) +
  geom_bar(stat="identity") +
  #geom_point(size=2, alpha=0.5) +
  labs(title="Bacterial Content",
       fill="Sample Type") +
  xlab(NULL) +
  facet_grid(~Sample_Type, scale = "free_x")+
  ylab("Percentage of Bacterial Reads (%)") +
  scale_fill_manual(labels = c("Tape Strip", "Biopsy"),
                    values = sample_type.cols) +
  theme_bw()

bac_percent


bac_percent_boxplot <- kraken_data %>%
  ggplot(mapping=aes(x=Sample_Type,
                     y=Bac_prcnt, 
                     fill=Sample_Type)) +
  geom_boxplot(alpha=.7, width = 0.5) +
  #geom_point(size=2, alpha=0.5) +
  geom_jitter(colour="black", size=1, width = 0.08) +
  labs(title="Bacterial Content",
       fill="Sample Type") +
  xlab(NULL) +
  ylab("Percentage of Bacterial Reads (%)") +
  scale_fill_manual(labels = c("Tape Strip", "Biopsy"),
                    values = sample_type.cols) +
  theme_bw()

bac_percent_boxplot

ggsave(bac_percent_boxplot,
       filename = file.path(out_dir, "generated_figures/05_bac_percent.png"),
       height = 4,
       width = 6)


# T-test between percentage of bacterial reads in tapes and biopsies

test <- 
  kraken_data %>%
    dplyr::select(Sample_Type, Participant, Bac_prcnt) %>%
    group_by(Participant) %>%
    arrange(Sample_Type) %>%
    summarize("diff" = diff(Bac_prcnt))


t.test(test$diff, mu = 0)


# Species Analysis -------------------------------------------------------------

top_species<-
kraken_bac_abund %>%
  group_by(Sample_ID) %>%
  filter(Rel_abund>5 & Sample_ID %in% meta$Sample) %>%
  slice_max(nr_reads_root_taxon, n=10)



ggplot(top_species, aes(x=Sample_ID, y = Rel_abund, fill = ID))+
  geom_bar(position="stack", stat="identity")+
  facet_grid(~Sample_Type, scales = "free_x")

ggsave(filename = file.path(out_dir, "generated_figures/05_bac_species_boxplot.png"),
       height = 4,
       width = 6)

  
# Quick Test of DE species 

species <- unique(top_species$ID)

diff <- 
  lapply(species, function(x) {
    AD <- kraken_bac_abund %>%
            filter(AD == "Yes" & ID == x & Sample_Type == "tape_strip" ) %>%
            pull(Rel_abund)
    
    No_AD <- kraken_bac_abund %>%
      filter(AD == "No" & ID == x & Sample_Type == "tape_strip" ) %>%
      pull(Rel_abund)

    out <- wilcox.test(AD, No_AD)
    
    return(data.frame("x" = x, "pval" = out$p.value))
  })


test <- bind_rows(diff) 
test$padj <- p.adjust(test$pval)



# Top Mean Abundance -----------------------------------------------------------

top_relabund_type<-
  kraken_bac_abund %>%
  filter(!is.na(Sample_Type)) %>%
  group_by(Sample_Type, ID) %>%
  summarize(mean_rel_abund = mean(Rel_abund), 
            sd_rel_abund = sd(Rel_abund)) %>%
  slice_max(mean_rel_abund, n=10)

#Write out the plots 
saveRDS(bac_percent, file.path(out_dir, "generated_rds/05_bac_prcnt_barplot.rds"))
saveRDS(bac_percent_boxplot, file.path(out_dir, "generated_rds/05_bac_prcnt_boxplot.rds")) 

# write out kraken data 
write_csv(top_relabund_type, file.path(out_dir, "generated_data/05_top_species.csv"))
write_csv(kraken_data, file.path(out_dir, "generated_data/05_kraken_results.csv"))
