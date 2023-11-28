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

# Define Color Schemes 

color_theme <- c("biopsy" = "#bebada", 
            "tape_strip" = "#ffffb3",
            "n_common" = "#8dd3c7",
            "n_biopsy" = "#bebada", 
            "n_tape" = "#ffffb3",
            "Assigned" = "#145291",
            "Unassigned_NoFeatures" = "#ffeb75",
            "Unassigned_Unmapped" = "#ffa04d", 
            "Unassigned_Overlapping_Length" = "#a70000") 

#Figure 1 ----------------------------------------------------------------------

#Figure 1A

fig1_featureCounts_formatted<-
fig1_featureCounts +
  facet_grid(~Sample_Type, labeller = as_labeller(c("tape_strip" = "Tape Strips",
                                                    "biopsy" = "Biopsy")),
             scales = "free_x")+
  theme_classic()+
  ylim(0,100)+
  theme(legend.position = c(0.4, 0.85),
        plot.margin = margin(0.2,0.8,0.2,0.2, "cm"),
        legend.title = element_text(size = 6), 
        legend.text  = element_text(size = 6),
        legend.key.size = unit(1, "lines"),
        axis.text.x = element_text(angle = 45,vjust = 0.8, hjust = 1))+
  guides(fill = guide_legend(nrow = 2)) +
  scale_fill_manual(values = color_theme)+
  labs(title="Read Assignments")


#Figure 1B

fig1_GenesPerMethod_formatted <- 
  fig1_GenesPerMethod +
     theme(legend.position = "none",
           plot.margin = margin(0.2,0.8,0.2,0.2, "cm"))+
     facet_grid(~cutoff, labeller = as_labeller(c("1" = "Read Count Cutoff:\n1",
                                                  "10" = "Read Count Cutoff:\n10",
                                                  "100" = "Read Count Cutoff:\n100")))+
  labs(title = "Identified Genes By Method")+
  scale_fill_manual(values = color_theme)


fig1_GenesPerMethod_formatted 


#Figure 1C


fig1_DiffvsCutoff_formatted <- 
  fig1_DiffvsCutoff+
    labs(title="Mean Difference in Number of Detected Genes \n by cutoff", 
         y="# of Gene Difference \n(Biopsy - Tape Strips)")+
    theme_classic()+
    theme(plot.margin = margin(0.2,0.8,0.2,0.2, "cm"))+
  ylim(-10000,10000)

#Figure 1D

fig1_PvalvsCutoff_formatted <- 
  fig1_PvalvsCutoff +
    theme_classic()+
    ylim(0,1)+
    theme(plot.margin = margin(0.2,0.8,0.2,0.2, "cm"))+
    labs(y = "False Discovery Rate (FDR)",
         title="P-value by cutoff*",
         caption="*P. values from paired t-test and FDR adjusted")+
    geom_text(label = "FDR = 0.05", y = 0.09, x = 470, size = 3, color = "red")

fig1_PvalvsCutoff_formatted

#Figure 1E

fig1_CommonGenesBoxplot <- 
  fig1_CommonGenesBoxplot +
    theme_classic() + 
    theme(legend.position = c(0.75, 0.80),
          plot.margin = margin(0.2,0.8,0.2,0.2, "cm"),
          legend.title = element_text(size = 6), 
          legend.text  = element_text(size = 6),
          legend.key.size = unit(1, "lines"),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"))+
    guides(fill = guide_legend(ncol= 2))+
    labs(x = "")
      
  
Figure1<-
plot_grid(fig1_featureCounts_formatted, 
          fig1_GenesPerMethod_formatted,
          fig1_DiffvsCutoff_formatted,
          fig1_PvalvsCutoff_formatted,
          fig1_CommonGenesBoxplot,
          ncol = 2, labels = c("A.", "B.", "C.", "D.", "E."))




save_plot("./results/generated_figures/Figure1.png", Figure1, base_height = 10, base_width = 10)


 
# Figure 2 ---------------------------------------------------------------------


fig2_pca_formatted <- 
  fig2_pca+
    stat_ellipse(aes(group=Sample_Type))

fig2_heatmap_formatted <- 
  fig2_heatmap
    

fig2_varPartALL_formatted <- 
  fig2_varPartALL+
  labs(title = "Variance Explained:\nAll Samples")+
  theme(plot.margin = margin(0.5,1,0.5,1, "cm"),
        axis.text.x = element_text(angle = 45))

fig2_varPartTape_formatted<-
  fig2_varPartTape +
    labs(title = "Variance Explained:\nTape-Strips")+
  theme(plot.margin = margin(0.5,1,0.5,1, "cm"),
        axis.text.x = element_text(angle = 45))

fig2_varPartBiops_formatted <- 
  fig2_varPartBiops+
    labs(title = "Variance Explained:\nBiopsies")+
  theme(plot.margin = margin(0.5,1,0.5,1, "cm"),
        axis.text.x = element_text(angle = 45))



fig2_scree_formatted <- 
  fig2_scree + 
    theme(plot.margin = margin(0.2,1.5,0.2,1.5, "cm"))


fig2_varpart<-plot_grid(fig2_varPartALL_formatted, fig2_varPartTape_formatted, fig2_varPartBiops_formatted, labels = c("D", NULL), ncol = 3)


fig2a<- plot_grid(fig2_heatmap_formatted$gtable, plot_grid(fig2_pca_formatted, fig2_scree_formatted, nrow =2, labels = c("B", "C")), rel_widths = c(0.6,0.4),
                  labels = c("A", NULL))

fig2_final<- plot_grid(fig2a, fig2_varpart, nrow = 2, rel_heights = c(0.65, 0.4))             
                 
fig2_final


save_plot("./results/generated_figures/Figure2.png", fig2_final, base_height = 10, base_width = 12)



# Figure 3 ---------------------------------------------------------------------

venn_ADEczVsControl <- read_rds("./results/generated_rds/12_venn_ADEczVsControl.rds")
venn_ADNoEczVsControl <- read_rds("./results/generated_rds/12_venn_ADNoEczVsControl.rds")
Venn_TapesVsBiopsy_NoAD <- read_rds("./results/generated_rds/12_Venn_TapesVsBiopsy_NoAD.rds")
Venn_SharedDEgenes_tapes <- read_rds("./results/generated_rds/12_Venn_SharedDEgenes_tapes.rds")

tbl_itch <- ggdraw() + draw_image("./results/generated_figures/13_Itch.png")
tbl_barrier <- ggdraw() + draw_image("./results/generated_figures/13_BarrierAll.png")
tbl_immune <- ggdraw() + draw_image("./results/generated_figures/13_ImmuneAll.png")

Fig3 <- plot_grid(tbl_itch, tbl_immune, tbl_barrier, Venn_TapesVsBiopsy_NoAD , nrow = 2,
                  labels = c("A", "B", 
                             "C", "D"),
                  label_x = c(-0.01, -0.04, -0.01, 0))


ggsave("./results/generated_figures/Figure3.png", Fig3)


# Figure 4 ---------------------------------------------------------------------

cats<-c(
"Tape_ControlVsADEcz",
"Biops_ControlVsADEcz",
"Tape_ControlVsADNoEcz",
"Biops_ControlVsADNoEcz")


de_tables <- lapply(cats, function(x){
  
  title_type <- ifelse(grepl("Tape", x), "Tape Strip", "Biopsy")
  title_condition1 <- ifelse(grepl("ControlVsADEcz", x), "Active AD", "Inactive AD")
  title_condition2 <- "Control"
  
  
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
  
  #print(read_delim(paste0("./results/11_res_",x,"_sig.csv")))
  data_out<-
    data %>%
        filter(gene %in% c(top_sig_up, top_sig_down)) %>%
        arrange(desc(log2FoldChange)) %>%
        dplyr::select(external_gene_name, baseMean, log2FoldChange, padj) %>%
        mutate("baseMean" = signif(baseMean, digits = 3)) %>%
        mutate("padj" = signif(padj, digits = 3)) %>%
        mutate("log2FoldChange" = signif(log2FoldChange, digits = 3))


  
  pal <- function(x) {
    f_neg <- scales::col_numeric(
      palette = c('blue','white' ),
      domain= c(0,-20),
      na.color = "black")

    f_pos <- scales::col_numeric(
      palette = c('white', 'red'),
      domain = c(0,20))

    ifelse(x < 0, f_neg(x), f_pos(x))
  }
  
  if(nrow(data_out)==0){
    warning("NO OVERLAPPING GENES")
    return(gt(data.frame("Error" = "No Overlapping Genes")))
  }

  out <-
    gt(data_out) %>%
    data_color(columns = log2FoldChange, fn = pal) %>%
    cols_label(
      external_gene_name = "Gene",
      padj = "Adjusted p-value"
    ) %>%
    cols_align(align = "center") %>%
    
    tab_header(title = md(paste0("DESeq - ", title_type)),
               subtitle = md(paste0("<span style='color:#f48f73'>", title_condition1, "</span>
                           vs 
                           <span style='color:#9067fc'>", title_condition2,"</span>"))) %>%
    tab_options(table.font.size = 75) %>%
    gtsave(paste0("results/generated_figures/09_tbl_Top10DEseq_", x,".png"), vwidth = 2000)
  
  })


plt_tables <- lapply(cats, function(x){
  plt <- ggdraw() + draw_image(paste0("results/generated_figures/09_tbl_Top10DEseq_", x, ".png"), valign = 1)
})


plt_tables <- plot_grid(plotlist=plt_tables, ncol = 4, labels = c("C","", "D", ""))

plt_venns <- plot_grid(venn_ADEczVsControl,
                       venn_ADNoEczVsControl, labels = c("A","B"))

plot_grid(plt_venns,
          plt_tables, nrow = 2, rel_heights = c(0.4, 0.7))


ggsave("results/generated_figures/Figure4.png")




# Figure 5 ---------------------------------------------------------------------

plt_bacPrcnt <- readRDS("./results/generated_rds/05_bac_prcnt_barplot.rds")
plt_bacPrcntBoxplot <- readRDS("./results/generated_rds/05_bac_prcnt_boxplot.rds")
df_topSpecies <- read_csv("./results/generated_data/05_top_species.csv")


df_topSpecies_tape <- 
  df_topSpecies %>%
    filter(Sample_Type == "tape_strip") %>%
    dplyr::select(ID, mean_rel_abund, sd_rel_abund) %>%
    mutate(mean_rel_abund = round(mean_rel_abund, 2)) %>%
    mutate(sd_rel_abund = round(sd_rel_abund,1)) %>%
    dplyr::rename("Mean\nRelative\nActivity (%)" = mean_rel_abund,
           "SD\nRelative\nActivity (%)" = sd_rel_abund)
  

df_topSpecies_biopsy <- 
  df_topSpecies %>%
    filter(Sample_Type == "biopsy") %>%
    dplyr::select(ID, mean_rel_abund, sd_rel_abund) %>%
    mutate(mean_rel_abund = round(mean_rel_abund, 2)) %>%
    mutate(sd_rel_abund = round(sd_rel_abund,1)) %>%
    dplyr::rename("Mean\nRelative\nActivity (%)" = mean_rel_abund,
           "SD\nRelative\nActivity (%)" = sd_rel_abund)


plt_topSpecies_biops <- gridExtra::tableGrob(df_topSpecies_biopsy, rows = NULL,
                                             theme = ttheme_default(base_size = 11,
                                                                    core=list(fg_params=list(lineheight=1)),
                                                                    padding = unit(c(3,3), "mm")))



plt_topSpecies_tape <- gridExtra::tableGrob(df_topSpecies_tape, rows = NULL,
                                            theme = ttheme_default(base_size = 11,
                                                                   core=list(fg_params=list(lineheight=1)),
                                                                   padding = unit(c(3,3), "mm")))



plt_bacPrcnt_formatted <- 
  plt_bacPrcnt +  
  geom_bar(stat="identity", aes(x=Sample, y = Bac_prcnt, fill = Sample_Type), color = "black")+
  scale_fill_manual(values = color_theme)+
  theme_classic()+
  theme(legend.position = c(0.7, 0.85),
        plot.margin = margin(1,1,01,1, "cm"),
        legend.title = element_text(size = 6), 
        legend.text  = element_text(size = 6),
        legend.key.size = unit(1, "lines"),
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 1))+
  labs(title = "Percentage of Bacterial Reads by Sample")


plt_bacPrcntBoxplot <- 
  plt_bacPrcntBoxplot +  
  scale_fill_manual(values = color_theme)+
  theme_classic()+
  theme(legend.position = c(0.7, 0.85),
        plot.margin = margin(1,1,01,1, "cm"),
        legend.title = element_text(size = 6), 
        legend.text  = element_text(size = 6),
        legend.key.size = unit(1, "lines"))+
  labs(title = "Percentage of Bacterial Reads by Method")
    



# Building the figure 

fig5top<-plot_grid(plt_bacPrcnt_formatted, plt_bacPrcntBoxplot, ncol = 2, labels = c("A", "B"))
fig5bottom<-plot_grid(plt_topSpecies_biops, plt_topSpecies_tape,  ncol = 2, labels = c("C. Top 10 Species: Biopsy", 
                                                                                      "D. Top 10 Species: Tape-Strips"), 
                      label_x = c(-0.21, -0.25),
                      label_y = c(0.95,0.95)
                      )

fig5 <- plot_grid(fig5top, fig5bottom, nrow = 2)

fig5

save_plot("./results/generated_figures/Figure5.png", fig5, base_height = 8, base_width = 10)




