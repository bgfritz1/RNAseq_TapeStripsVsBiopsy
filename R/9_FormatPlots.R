#This will import all the plots so I can make some nice summarized, faceted
#figures

# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggtext)
library(cowplot)
library(gridExtra)
library(gt)
library(ggpubr)
library(gridtext)
library(grid)
library(gtable)


# Import the Figures 

fig1_featureCounts <- readRDS("./results/generated_rds/06_featureCounts_figure.rds")
fig1_GenesPerMethod <- readRDS("./results/generated_rds/01_n_genes_vs_method.rds")
fig1_PvalvsCutoff <- readRDS("./results/generated_rds/01_PvalVsCutoff.rds")
fig1_DiffvsCutoff <- readRDS("./results/generated_rds/01_DiffVsCutoff.rds")
fig1_CommonGenesBoxplot <- readRDS("./results/generated_rds/01_CommonGenesBoxplot.rds")

fig2_pca <- readRDS("./results/generated_rds/02_PCA_plot_AllSamples.rds")
fig2_heatmap <- readRDS("./results/generated_rds/02_distance_heatmap.rds")
fig2_scree <- readRDS("./results/generated_rds/02_scree_plot.rds")
fig2_varPartALL<- readRDS("./results/generated_rds/10_vp_plot_all.rds")
fig2_varPartTape<- readRDS("./results/generated_rds/10_vp_plot_tape.rds")
fig2_varPartBiops<- readRDS("./results/generated_rds/10_vp_plot_biops.rds")

figSupp_Eczpca <- readRDS("./results/generated_rds/02_PCA_plot_AllSamples_ecz.rds")

# Define Color Schemes 

color_theme <- c("biopsy" = "#bebada", 
            "tape_strip" = "#ffffb3",
            "n_common" = "#8dd3c7",
            "n_biopsy" = "#bebada", 
            "n_tape" = "#ffffb3",
            "Residuals" = "grey",
            "Assigned" = "#145291",
            "Eczema_dorsalhand" = "#ff8c00",
            "AD" = "#dc6653",
            "nonAD" = "#fdfab4",
            "Sample_Type" = "#00bfc4",
            "Unassigned_NoFeatures" = "#ffeb75",
            "Unassigned_Unmapped" = "#ffa04d", 
            "Unassigned_Overlapping_Length" = "#a70000") 

#Figure 1 ----------------------------------------------------------------------

#Figure 1A

fig1_featureCounts_formatted<-
fig1_featureCounts +
  facet_grid(~Sample_Type, labeller = as_labeller(c("tape_strip" = "Tape-Strips",
                                                    "biopsy" = "Biopsy")),
             scales = "free_x")+
  theme_classic()+
  scale_y_continuous(limits = c(0,105), expand = c(0, 0))+
  theme(legend.position = c(0.5, 0.85),
        plot.margin = margin(0.2,0.8,0.2,0.2, "cm"),
        legend.title = element_blank(),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(1, "lines"),
        legend.box.background = element_rect(colour = "black", fill = "white"),
        legend.margin=margin(c(1,5,5,5)),
        axis.text.x = element_text(angle = 45,vjust = 0.8, hjust = 1))+
  guides(fill = guide_legend(nrow = 2)) +
  scale_fill_manual(values = color_theme, 
                    labels = c("Unassigned_NoFeatures" = "No Feature",
                               "Unassigned_Unmapped" = "Unmapped",
                               "Unassigned_Overlapping_Length" = "Multiple Features"))+
  labs(title="Read Assignments",
       y = "Number of Read Pairs (M)",
       x = "Sample")


#Figure 1B

#Ugly way to add asterisk 
dat_txt<-data.frame(count=200,Sample_Type="biopsy",gene_count=43000, lab="*")

fig1_GenesPerMethod_formatted <- 
  fig1_GenesPerMethod +
     theme(legend.position = c(0.86, 0.7),
          plot.margin = margin(0.2,0.8,0.2,0.2, "cm"),
          legend.title = element_text(size = 8), 
          legend.text  = element_text(size = 8),
          legend.key.size = unit(1, "lines"),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black", fill = "white"),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
     facet_grid(~cutoff, labeller = as_labeller(c("10" = "Read Count Cutoff:\n10",
                                                  "75" = "Read Count Cutoff:\n75",
                                                  "200" = "Read Count Cutoff:\n200")))+
  labs(title = "Identified Genes By Method",
       y = "Number of Detected Genes")+
  ylim(0,45000)+
  guides(fill = guide_legend(nrow= 2, title = "Sample Type", override.aes = list(shape=NA)))+
  scale_fill_manual(values = color_theme, labels = c("tape_strip" = "Tape-Strips", "biopsy" = "Biopsy"))+
  stat_compare_means(method = "wilcox.test", paired = T, label = "p.format", label.y.npc = 0.95)+
  geom_text(data=dat_txt, label = c("","*",""), nudge_x = -0.2, nudge_y = 800)


fig1_GenesPerMethod_formatted 


#Figure 1C


fig1_DiffvsCutoff_formatted <- 
  fig1_DiffvsCutoff+
    labs(title="Mean Difference in Number of Detected Genes \n by cutoff", 
         y="Number of Gene Difference \n(Biopsy - Tape-Strips)")+
    theme_classic()+
    theme(plot.margin = margin(0.2,0.8,0.2,0.2, "cm"))+
    scale_x_continuous(breaks = c(10,75,100,200,300,400,500))+
    ylim(-10000,10000)+
    geom_vline(xintercept = 75, color = "black", linetype="dashed")+
    geom_vline(xintercept = 4, color = "black", linetype="dashed")+
    annotate("text", x = 140, y = 10000, label = "FDR<=0.05", color = "black", size = 3)+
    annotate("text", x = 35, y = 10000, label = "NS", color = "black", size = 3)+
    annotate("rect", xmin = 4, xmax = 75, ymin = -Inf, ymax = Inf, alpha = 0.3)
  
fig1_DiffvsCutoff_formatted

#Figure 1D

fig1_PvalvsCutoff_formatted <- 
  fig1_PvalvsCutoff +
    theme_classic()+
    ylim(0,1)+
    xlim(10,525)+
    geom_vline(xintercept = 75, color = "black", linetype="dashed", alpha = 0.5)+
    geom_vline(xintercept = 4, color = "black", linetype="dashed", alpha = 0.5)+
    scale_x_continuous(breaks = c(10,75,100,200,300,400,500))+
    theme(plot.margin = margin(0.2,0.8,0.2,0.2, "cm"))+
    labs(y = "False Discovery Rate (FDR)",
         title="P-value by cutoff*",
         caption="*P. values from paired t-test and FDR adjusted")+
    geom_text(label = "FDR = 0.05", y = 0.1, x = 450, size = 3, color = "red")+
    annotate("rect", xmin = 4, xmax = 75, ymin = -Inf, ymax = Inf, alpha = 0.3)+
    annotate("text", x = 35, y = 1, label = "NS", color = "black", size = 3)


fig1_PvalvsCutoff_formatted

#Figure 1E

fig1_CommonGenesBoxplot_formatted <- 
  fig1_CommonGenesBoxplot +
    #theme_classic() + 
    theme(legend.position = c(0.8355, 0.80),
          plot.margin = margin(0.2,0.8,0.2,0.2, "cm"),
          legend.title = element_text(size = 8), 
          legend.text  = element_text(size = 8),
          legend.key.size = unit(1, "lines"),
          legend.box.background = element_rect(colour = "black", fill = "white"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank())+
    labs(y = "Number of Genes")+
    guides(fill = guide_legend(nrow= 3, title = "Sample Type"))+
    scale_fill_manual(values = color_theme, labels = c("n_common" = "Common \n(Tape-Strip + Biopsy)",
                                                       "n_tape" = "Tape-Strip",
                                                       "n_biopsy" = "Biopsy"))

fig1_CommonGenesBoxplot_grid <- plot_grid(fig1_CommonGenesBoxplot_formatted, NULL, rel_widths = c(0.6,0.3))
  
Figure1_ABCD<-
plot_grid(fig1_featureCounts_formatted, 
          fig1_GenesPerMethod_formatted,
          fig1_PvalvsCutoff_formatted,
          fig1_DiffvsCutoff_formatted,
          ncol = 2, labels = c("A.", "B.", "C.", "D."))

Figure1_final <- plot_grid(Figure1_ABCD, fig1_CommonGenesBoxplot_grid , 
                           rel_heights = c(0.66, 0.33),
                           labels = c("","E."),
                           nrow = 2) 




save_plot("./results/generated_figures/Figure1.png", Figure1_final, base_height = 10, base_width = 10)


 
# Figure 2 ---------------------------------------------------------------------


varplot_theme <- 
  theme(#legend.position = c(0.86, 0.7),
    legend.position = "bottom",
    plot.margin = margin(0.2,0.2,0.2,0.2, "cm"),
    legend.title = element_text(size = 8), 
    legend.text  = element_text(size = 8),
    legend.key.size = unit(1, "lines"),
    legend.background = element_blank(),
    #legend.box.background = element_rect(colour = "black", fill = "white"),
    legend.box.background = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank())
  

fig2_pca_formatted <- 
  fig2_pca+
    stat_ellipse(aes(group=Sample_Type))+
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(),
        plot.margin = margin(2,1,2,1, "cm"))+
  guides(color=guide_legend(override.aes = list(shape = 15)))


fig2_varPartALL_formatted <- 
  fig2_varPartALL+
  labs(title = "Variance Explained:\nAll Samples")+
  
  scale_fill_manual(values = color_theme, labels = c("Sample_Type" = "Sample Type\n(Tape-Strip/Biopsy)",
                                                     "AD" = "AD (AD/non-AD)",
                                                     "Eczema_dorsalhand" = "Eczema on \nDorsal Hand (Yes/No)"))+
  guides(fill = guide_legend(nrow= 1, title = "Variance Source",title.position = "top"))+
  varplot_theme 

fig2_varPartTape_formatted<-
  fig2_varPartTape +
    labs(title = "Variance Explained:\nTape-Strips")+
    scale_fill_manual(values = color_theme, labels = c("Sample_Type" = "Sample Type\n((Tape-Strip/Biopsy)",
                                                       "AD" = "AD (AD/non-AD)",
                                                       "Eczema_dorsalhand" = "Eczema on \nDorsal Hand (Yes/No)"))+
    guides(fill = guide_legend(nrow= 2, title = "Variance Source",title.position = "top"))+
    varplot_theme+
  theme(legend.position = "none")

fig2_varPartBiops_formatted <- 
  fig2_varPartBiops+
    labs(title = "Variance Explained:\nBiopsies")+
    scale_fill_manual(values = color_theme, labels = c("Sample_Type" = "Sample Type\n((Tape-Strip/Biopsy)",
                                                       "AD" = "AD (AD/non-AD)",
                                                       "Eczema_dorsalhand" = "Eczema on \nDorsal Hand (Yes/No)"))+
    guides(fill = guide_legend(nrow= 2, title = "Variance Source",title.position = "top"))+
    varplot_theme +
    theme(legend.position = "none")

fig2_varLegend <- get_legend(fig2_varPartALL_formatted)

fig2_scree_formatted <- 
  fig2_scree + 
    theme(plot.margin = margin(0.2,1.5,0.2,1.5, "cm"))


fig2_varpart <- plot_grid(fig2_varPartALL_formatted + theme(legend.position = "none"), 
                          fig2_varPartTape_formatted + theme(legend.position = "none"), 
                          fig2_varPartBiops_formatted + theme(legend.position = "none"),
                          labels = c("C.","D.","E."), ncol = 3)

fig2_varfinal <- plot_grid(fig2_varpart, fig2_varLegend, nrow = 2, rel_heights = c(0.8, 0.2))

fig2a<- plot_grid(fig2_heatmap, fig2_pca_formatted, rel_widths = c(0.6,0.4),
                  labels = c("A.", "B."))

fig2_final<- plot_grid(fig2a, fig2_varfinal, nrow = 2, rel_heights = c(0.65, 0.4))             
                 
fig2_final


save_plot("./results/generated_figures/Figure2.png", fig2_final, base_height = 9, base_width = 10, bg = "white")



# Figure 3 ---------------------------------------------------------------------

############################################
# SEE 13_EnrichmentAnalysis for Figure 3 ###
############################################


# Figure 4 ---------------------------------------------------------------------

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
    domain= c(0,-20),
    na.color = "black")
  
  f_pos <- scales::col_numeric(
    palette = c('white', color_up),
    domain = c(0,20))
  
  ifelse(x < 0, f_neg(x), f_pos(x))
}

# Format Grobs 

formatTables <- function(x, title){
  
  #Generate the background colors for the LFC
  bg_cols <- pal(x$Log2FC, color_up = "red", color_down= "#7065fb")
  
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



#Import Venns 

venn_ADEczVsnonAD <- read_rds("./results/generated_rds/12_venn_ADEczVsnonAD.rds")
venn_ADNoEczVsnonAD <- read_rds("./results/generated_rds/12_venn_ADNoEczVsnonAD.rds")
Venn_TapesVsBiopsy_NoAD <- read_rds("./results/generated_rds/12_Venn_TapesVsBiopsy_NoAD.rds")
Venn_SharedDEgenes_tapes <- read_rds("./results/generated_rds/12_Venn_SharedDEgenes_tapes.rds")


cats<-c(
  "Tape_nonADVsADEcz",
  "Biops_nonADVsADEcz",
  "Tape_nonADVsADNoEcz",
  "Biops_nonADVsADNoEcz")



de_tables <- lapply(cats, function(x){
  
  title_type <- ifelse(grepl("Tape", x), "Tape-Strips", "Biopsy")
  title_condition1 <- ifelse(grepl("nonADVsADEcz", x), "Active AD", "Inactive AD")
  title_condition2 <- "non-AD"
  
  
  data <- read_delim(paste0("./results/generated_data/11_res_",x,"_sig.csv"))   
  
  top_sig_up <- 
    data %>% 
    filter(log2FoldChange>0) %>% 
    slice_min(padj, n = 10) %>%
    pull(gene)
  
  top_sig_down <- 
    data %>% 
    filter(log2FoldChange<0) %>% 
    slice_min(padj, n = 10) %>%
    pull(gene)
  
  
  data_out<-
    data %>%
    filter(gene %in% c(top_sig_up, top_sig_down)) %>%
    arrange(desc(log2FoldChange)) %>%
    dplyr::select(external_gene_name, baseMean, log2FoldChange, padj) %>%
    mutate("baseMean" = signif(baseMean, digits = 0)) %>%
    mutate("padj" = signif(padj, digits = 3)) %>%
    mutate("log2FoldChange" = signif(log2FoldChange, digits = 3))
  
  
  
  colnames(data_out) <- c("Gene", "Mean\nExpression", "Log2FC", "Adjusted\nP-value")
  
  
  title = paste0(title_type,"<br><span style='color:#f48f73'>", title_condition1, "</span>
                   vs 
                 <span style='color:#9067fc'>", title_condition2,"</span>")
  
  out <- formatTables(data_out, title = title)
  
  return(out)
  
})

Fig4C <- arrangeGrob(justify(de_tables[[1]], just = -4), justify(de_tables[[2]], just = 1), ncol = 2)
Fig4D <- arrangeGrob(justify(de_tables[[3]], just = 0.05), justify(de_tables[[4]], just = 0.92), ncol = 2)


g <- list(venn_ADEczVsnonAD + theme(plot.margin = margin(0.5,2,0.5,2, "cm")),
          venn_ADNoEczVsnonAD + theme(plot.margin = margin(0.5,2,0.5,2, "cm")),
          Fig4C, Fig4D)

figure_mat <- rbind(c(1,1,2,2), c(3,3,4,4),c(3,3,4,4))

Fig4<-
  arrangeGrob(grobs = g, layout_matrix = figure_mat) |>
  as_ggplot() +
  draw_plot_label(label = c("A.", "B.", "C.", "D."), 
                  x = c(0, 0.48, 0 , 0.48),
                  y = c(1, 1, 0.62, 0.62),
                  size = 25)

ggsave("results/generated_figures/Figure4.png", width = 15.5, height = 15)


# Figure 5 ---------------------------------------------------------------------

plt_bacPrcnt <- readRDS("./results/generated_rds/05_bac_prcnt_barplot.rds")
plt_bacPrcntBoxplot <- readRDS("./results/generated_rds/05_bac_prcnt_boxplot.rds")
df_topSpecies <- read_csv("./results/generated_data/05_top_species.csv")


  
df_topSpecies_tapeNoAD <- 
  df_topSpecies %>%
  filter(Sample_Type == "tape_strip" & AD == "No") %>%
  dplyr::select(ID, AD, mean_rel_abund, sd_rel_abund) %>%
  mutate(mean_rel_abund = round(mean_rel_abund, 2)) %>%
  mutate(sd_rel_abund = round(sd_rel_abund,1)) %>%
  dplyr::rename("Mean\nRelative\nActivity (%)" = mean_rel_abund,
                "SD\nRelative\nActivity (%)" = sd_rel_abund)

df_topSpecies_tapeAD <- 
  df_topSpecies %>%
  filter(Sample_Type == "tape_strip" & AD == "Yes") %>%
  dplyr::select(ID, AD, mean_rel_abund, sd_rel_abund) %>%
  mutate(mean_rel_abund = round(mean_rel_abund, 2)) %>%
  mutate(sd_rel_abund = round(sd_rel_abund,1)) %>%
  dplyr::rename("Mean\nRelative\nActivity (%)" = mean_rel_abund,
                "SD\nRelative\nActivity (%)" = sd_rel_abund)


plt_topSpecies_tapeAD <- gridExtra::tableGrob(df_topSpecies_tapeAD, rows = NULL,
                                             theme = ttheme_default(base_size = 11,
                                                                    core=list(fg_params=list(lineheight=1),
                                                                              bg_params=list(fill="white")),
                                             padding = unit(c(3,3), "mm")))



plt_topSpecies_tapenoAD <- gridExtra::tableGrob(df_topSpecies_tapeNoAD, rows = NULL,
                                            theme = ttheme_default(base_size = 11,
                                                                   core=list(fg_params=list(lineheight=1),
                                                                             bg_params=list(fill="white")),
                                                                   padding = unit(c(3,3), "mm")))



plt_bacPrcnt_formatted <- 
  plt_bacPrcnt +  
  geom_bar(stat="identity", aes(x=Sample, y = Bac_prcnt, fill = AD), color = "black")+
  scale_fill_manual(values = c("Yes" = "#dc6653",
                               "No" = "#fdfab4"),
                    labels = c("Yes" = "AD",
                               "No" = "non-AD"))+
  theme_classic()+
  theme(text = element_text(size=12),
        legend.position = c(0.7, 0.7),
        plot.margin = margin(1,1,01,1, "cm"),
        axis.text.x = element_text(angle = 55, vjust = 0.8, hjust = 1))+
  labs(title = "Percentage of Bacterial Reads by Sample",
       y = "Percentage of total reads \nfrom bacteria (%)")
  


plt_bacPrcntBoxplot <- 
  plt_bacPrcntBoxplot +  
  scale_fill_manual(values = color_theme,
                    labels = c("tape_strip" = "Tape-Strips",
                               "biopsy" = "Biopsy"))+
  theme_classic()+
  theme(text = element_text(size=12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.7, 0.75),
        plot.margin = margin(1,1,01,1, "cm"),
        legend.title = element_text(size = 9), 
        legend.text  = element_text(size = 9),
        legend.key.size = unit(1, "lines"),
        legend.box.background = element_rect(colour = "black", fill = "white"))+
  labs(title = "Percentage of Bacterial Reads by Method",
       y = "Percentage of total reads \nfrom bacteria (%)")+
  stat_compare_means(paired = T, label.x = 1.5, label.y=17, label = "p.format")
    



# Building the figure 

fig5top<-plot_grid(plt_bacPrcnt_formatted, plt_bacPrcntBoxplot, ncol = 2, labels = c("A", "B"))
fig5bottom<-plot_grid(plt_topSpecies_tapeAD, plt_topSpecies_tapenoAD,  ncol = 2, labels = c("C. Top 10 Species: Tape-Strips AD", 
                                                                                      "D. Top 10 Species: Tape-Strips non-AD"), 
                      label_x = c(-0.21, -0.25),
                      label_y = c(0.95,0.95)
                      )

fig5 <- plot_grid(fig5top, fig5bottom, nrow = 2)

fig5

save_plot("./results/generated_figures/Figure5.png", fig5, base_height = 8, base_width = 10)


#=============================
# Supplementary Figures 
#============================= 


# PCA plot by AD and Lesion Status


figSupp_Eczpca_formatted <- 
  figSupp_Eczpca +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(),
        plot.margin = margin(2,1,2,1, "cm"))+
  guides(color=guide_legend(override.aes = list(shape = 15)))

figSupp_Eczpca_formatted 

save_plot("./results/generated_figures/FigureS1.png", figSupp_Eczpca_formatted, base_height = 8, base_width = 8)

# S_aureus counts and abundance 

plt_kraken_saureus_reads <- 
  readRDS("./results/generated_rds/05_Saureus_reads_barplot.rds")+
  facet_grid(~Sample_Type, labeller = as_labeller(c("tape_strip" = "Tape-Strips",
                                                    "biopsy" = "Biopsy")),
             scales = "free_x")+
  geom_bar(stat="identity", aes(fill = Subject_type), color = "black")+
  theme_classic()+
  theme(legend.position = c(0.86, 0.7),
        plot.margin = margin(0.2,0.8,0.2,0.2, "cm"),
        legend.title = element_text(size = 8), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(1, "lines"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black", fill = "white"),
        plot.title = element_markdown(),
        axis.text.x = element_text(angle = 55, vjust = 0.8, hjust = 1))+
scale_fill_manual(values = color_theme,
                  labels = c("AD" = "AD", 
                             "nonAD" = "non-AD"))+
  labs(y = "Number of Reads", title = "Assigned reads to *Staphylococcus aureus*",
       x = "Sample")+
  guides(fill = guide_legend(nrow = 2, title = "Subject Type"))

plt_kraken_saureus_relabund <- 
  readRDS("./results/generated_rds/05_Saureus_relabund_barplot.rds")+
  facet_grid(~Sample_Type, labeller = as_labeller(c("tape_strip" = "Tape-Strips",
                                                    "biopsy" = "Biopsy")),
             scales = "free_x")+
  geom_bar(stat="identity", aes(fill = Subject_type), color = "black")+
  theme_classic()+
  theme(legend.position = c(0.86, 0.7),
        plot.margin = margin(0.2,0.8,0.2,0.2, "cm"),
        legend.title = element_text(size = 8), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(1, "lines"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black", fill = "white"),
        plot.title = element_markdown(),
        axis.text.x = element_text(angle = 55, vjust = 0.8, hjust = 1))+
  scale_fill_manual(values = color_theme,
                    labels = c("AD" = "AD", 
                               "nonAD" = "non-AD"))+
  labs(y = "Relative activity (%)", title = "Relative activity of *Staphylococcus aureus*",
       x = "Sample")+
  guides(fill = guide_legend(nrow = 2, title = "Subject Type"))

fig_SuppSA <- plot_grid(plt_kraken_saureus_reads, plt_kraken_saureus_relabund , nrow = 2, align = "hv", labels = c("A.", "B."))

save_plot("./results/generated_figures/FigureS2.png", fig_SuppSA, base_height = 8, base_width = 8)
