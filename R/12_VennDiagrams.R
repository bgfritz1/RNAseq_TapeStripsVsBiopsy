## -----------------------------------------------------------------------------
##
## Script name: 12_VennDiagrams.R
##
## Purpose of script: 
## This script generates the RDS file that are used 
## to make the venn diagram figures, which illustrate the DEGs
## that are shared/distinct between different comparisons. 
## 
## 
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
library(ggVennDiagram)
library(ggtext)

#-------------------------------------------------------------------------------
# Import the Data --------------------------------------------------------------
#-------------------------------------------------------------------------------

# Read in the data -------------------------------------------------------------

categories <- c(
  "TapeVsBiopsy_NoADorEcz",
  "TapeVsBiopsy_ADandEcz",
  "Tape_nonADVsADEcz",
  "Tape_nonADVsADNoEcz",
  "Biops_nonADVsADEcz",
  "Biops_nonADVsADNoEcz")

data <- lapply(paste0(out_dir, "/generated_data/11_res_", categories, "_sig.csv"), read_delim)


names(data) <- categories


# Venn Diagrams ----------------------------------------------------------------

#Define a function to make the plots

plotVennDiagrams <- function(set1, set2, labels=c("set1", "set2"), colours = c("#fdfab4", "#dc6653", "#6fa0a3")){
  require(ggVennDiagram)

  set1_up <- filter(set1, log2FoldChange>0) %>% pull(gene)
  set2_up <- filter(set2, log2FoldChange>0) %>% pull(gene)
  set1_down <- filter(set1, log2FoldChange<0) %>% pull(gene)
  set2_down <- filter(set2, log2FoldChange<0) %>% pull(gene)
  
  up <- list(set1_up, set2_up)
  down <- list(set1_down, set2_down)
  
  names(up) <- labels
  names(down) <- labels
  
  venn_up <- Venn(up)
  venn_down <- Venn(down)

  
  data_venn_up <- process_data(venn_up)
  data_venn_down <- process_data(venn_down)
  
  region_data_up <- venn_region(data_venn_up)
  region_data_down <- venn_region(data_venn_down)
  
  region_data_up$item <- unlist(lapply(region_data_up$item, length))
  region_data_down$item <- unlist(lapply(region_data_down$item, length))
  
  region_data_combined <- bind_rows(region_data_up, region_data_down, .id = "set")
  
  plt <- 
    ggplot() +
    geom_sf(aes(fill=name), data = venn_region(data_venn_up)) +
    scale_fill_manual(values = colours)+
    geom_sf(size = 2,color = "black", data = venn_setedge(data_venn_up), show.legend = F) +
    geom_sf_text(aes(label = name), data = venn_setlabel(data_venn_up)) +
    geom_sf_label(aes(label = item), color = "red", fontface = "bold", nudge_y = 50, data =  region_data_up) +
    geom_sf_label(aes(label = item), color = "blue", fontface = "bold", data = region_data_down)+
    theme_void()+
    theme(legend.position = "none")
  
  
  print(paste0("Shape input set1: ", dim(set1)))
  print(paste0("Shape input set2: ", dim(set2)))
  
  print(paste0("Set1_up= ", length(set1_up)))
  print(paste0("Set1_down= ", length(set1_down)))
  print(paste0("Set2_up= ", length(set2_up)))
  print(paste0("Set2_down= ", length(set2_down)))

  
  
  return(plt)
}


# Lesional AD vs Non-AD

venn_ADEczVsnonAD <- 
plotVennDiagrams(data[["Tape_nonADVsADEcz"]], 
                 data[["Biops_nonADVsADEcz"]],
                 labels = c("Tape Strip", "Biopsy"))+
  labs(title = "DEGs: 
       **<span style='color:#e0433f;'>Active AD</span>** vs.
       **<span style='color:#3246fb'>Non-AD</span>**")+
  theme(plot.title = element_markdown(hjust = 0.5))

# Non-Lesional AD vs Non-AD

venn_ADNoEczVsnonAD <- 
  plotVennDiagrams(data[["Tape_nonADVsADNoEcz"]], 
                   data[["Biops_nonADVsADNoEcz"]],
                   labels = c("Tape Strip", "Biopsy"))+
  labs(title = "DEGs: 
       **<span style='color:#e0433f;'>Active AD</span>** vs.
       **<span style='color:#3246fb;'>Non-AD</span>**")+
  theme(plot.title = element_markdown(hjust = 0.5))

# Tapes Vs. Biopsy 

Venn_TapesVsBiopsy_NoAD <- 
  plotVennDiagrams(data[["TapeVsBiopsy_NoADorEcz"]], 
                   data[["TapeVsBiopsy_ADandEcz"]],
                   labels = c("No AD", "AD"))+
  labs(title = "DEGs: 
       **<span style='color:#e0433f;'>Tape-strips</span>** vs.
       **<span style='color:#3246fb;'>Biopsies</span>**")+
  theme(plot.title = element_markdown(hjust = 0.5))


# Shared DE genes - Lesional and non-lesional AD - Tapes 

Venn_SharedDEgenes_tapes <- 
  plotVennDiagrams(data[["Tape_nonADVsADNoEcz"]], 
                   data[["Tape_nonADVsADEcz"]],
                   labels = c("EczemaVSnonAD", "Non-EczemaVSnonAD"))+
   labs(title = "DEGs (Tapes): 
        **<span style='color:#e0433f;'>Active AD</span>** vs.
        **<span style='color:#3246fb;'>Non-AD</span>**")+
  theme(plot.title = element_markdown(hjust = 0.5))


out <- 
cowplot::plot_grid(venn_ADEczVsnonAD, 
                   venn_ADNoEczVsnonAD, 
                   Venn_TapesVsBiopsy_NoAD,
                   Venn_SharedDEgenes_tapes, nrow = 2)


ggsave(file.path(out_dir, "./generated_figures/12_VennDiagrams.png"), out, width = 10)

#Export Venndiagrams

write_rds(venn_ADEczVsnonAD, file.path(out_dir, "./generated_rds/12_venn_ADEczVsnonAD.rds"))
write_rds(venn_ADNoEczVsnonAD,  file.path(out_dir, "./generated_rds/12_venn_ADNoEczVsnonAD.rds"))
write_rds(Venn_TapesVsBiopsy_NoAD, file.path(out_dir, "./generated_rds/12_Venn_TapesVsBiopsy_NoAD.rds"))
write_rds(Venn_SharedDEgenes_tapes, file.path(out_dir, "./generated_rds/12_Venn_SharedDEgenes_tapes.rds"))


