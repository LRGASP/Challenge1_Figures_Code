#######################################
######### FIGURES CHALLENGE 1 #########
#######################################
## author: Francisco J. Pardo-Palacios, f.pardo.palacios@gmail.com
## author:Ana Conesa, ana.conesa@csic.es
## Last modified: April 25th 2023
#######################################

# LIBRARIES
#################
library(ggplot2)
library(ggrepel)
library(patchwork)
library(grid)
library(MetBrewer)
library(RColorConesa)
library(ggpubr)
library(scales)
library(huxtable)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(gridExtra)
library(reshape2)
library(gt)
library(ggpmisc)
library(UpSetR)
library(ComplexUpset)
library(ggplot2movies)
library(tidyr)
library(grid)
library(ggbreak) 
library(gg.gap)
library(ggcorrplot)
library(viridis)
library(fmsb)

outdir = "output/main"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# FUNCTIONS
#################
source("Functions_Supplementary_Figures_Challenge1_v4.R")

# FIGURES
#################

# Figure 2a. Detection
#########################
FirstPanelChl1 (data_sample = "WTC11", outdir = outdir, ylims = c(160000, 350000), xlims = c(12000, 28000))
# Plots are generated separately and composed in PowerPoint.

# Figure 2b. Overlap
#########################
Fig2b <- agreement.pipelines(data_sample = "WTC11_results/WTC11")
ggsave(file=paste0(outdir, "/Fig2b.svg"), plot=Fig2b, width=9, height=4)

# Figure 2c. SIRVs
#########################
data_sample = "WTC11_results/WTC11"
code <- read.csv("Challenge1_Figures_Data/code.csv", header = TRUE)
code$Label <- paste(code$Library_Preps, code$Platform, sep="-")
code$Lib_Plat <- paste(code$Library_Preps, code$Platform, sep="-")
spliced_SIRV_metrics <- read.csv(paste0("Challenge1_Figures_Data/", data_sample,".splicedSIRVS_metrics.csv"), sep=",", header=T) %>% t() %>% as.data.frame()
spliced_SIRV_metrics <- merge(spliced_SIRV_metrics, code, by.x=0, by.y="pipelineCode")
spliced_SIRV_metrics[, "F1 score"] <- apply(spliced_SIRV_metrics, 1, function(x){
  s=as.numeric(x["Sensitivity"])
  p=as.numeric(x["Precision"])
  (2*s*p)/(s+p)
})
pivoted_spliced <- pivot_longer(spliced_SIRV_metrics,
                                cols=c("Sensitivity","Precision", "F1 score"),
                                names_to = "metrics",
                                values_to = "value")

pivoted_spliced$metrics  <- pivoted_spliced$metrics %>% factor(levels = c("Sensitivity","Precision", "F1 score"),
                                                               labels = c("Sensitivity","Precision", "F1 score"))

pSIRVspliced <- ggplot(pivoted_spliced, aes(x=Label, y=value)) +
  geom_segment( aes(x=Label, xend=Label, y=0, yend=value, color=Lib_Plat), size=0.8) +
  geom_point( size=2, aes( shape=Data_Category, color=Lib_Plat ))  +
  facet_grid( metrics ~ Alias, scales = "free", space = "free_x", switch = "y"  ) +
  pub_theme +
  scale_color_manual(values = libplat.palette) +
  labs(x="", y="") +
  theme(legend.position="none") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),limits = c(0, NA), position="right")

ggsave(file=paste0(outdir, "/Fig2c_up.svg"), plot=pSIRVspliced, width=9, height=4)

#### same for unspliced

unspliced_SIRV_metrics <- read.csv(paste0("Challenge1_Figures_Data/", data_sample,".unsplicedSIRVS_metrics.csv"), sep=",", header=T) %>% t()
unspliced_SIRV_metrics <- merge(unspliced_SIRV_metrics, code, by.x=0, by.y="pipelineCode")

unspliced_SIRV_metrics[, "F1 score"] <- apply(unspliced_SIRV_metrics, 1, function(x){
  s=as.numeric(x["Sensitivity"])
  p=as.numeric(x["Precision"])
  (2*s*p)/(s+p)
})

pivoted_unspliced <- pivot_longer(unspliced_SIRV_metrics,
                                  cols=c("Sensitivity","Precision", "F1 score"),
                                  names_to = "metrics",
                                  values_to = "value")
pivoted_unspliced$metrics  <- pivoted_unspliced$metrics %>% factor(levels = c("Sensitivity","Precision", "F1 score"),
                                                                   labels = c("Sensitivity","Precision", "F1 score"))

pSIRVunspliced <- ggplot(pivoted_unspliced, aes(x=Label, y=value)) +
  geom_segment( aes(x=Label, xend=Label, y=0, yend=value, color=Lib_Plat), size=0.8) +
  geom_point( size=2, aes( shape=Data_Category, color=Lib_Plat ))  +
  facet_grid( metrics ~ Alias, scales = "free", space = "free_x", switch = "y"  ) +
  pub_theme +
  scale_color_manual(values = libplat.palette) +
  labs(x="", y="") +
  theme(legend.position="none") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),limits = c(0, NA), position="right")

ggsave(file=paste0(outdir, "/Fig2c_down.svg"), plot=pSIRVunspliced, width=9, height=4)


# Figure 2d. Radarplot
#########################
radar.simulation (species = "human", directory = "Challenge1_Figures_Data/Simulations/", pdf = paste0(outdir, "/Fig2d"))

radar.simulation (species = "mouse", directory = "Challenge1_Figures_Data/Simulations/", pdf = paste0(outdir, "/Fig2d_mouse"))


## Figure 2e. Evaluation against GENCODE
########################################
pa_GENCODE <- read.csv("Challenge1_Figures_Data/GENCODE_manualAnnot/presence_absence.GENCODE_loci_2.csv", sep=",", header = T) [,1:51] # Presence absence analysis of all transcripits of the 50 loci in pipelines evaluated against manual annotation. 
pa.WTC11 <- read.csv("Challenge1_Figures_Data/WTC11_results/WTC11_comparison.pa.csv", as.is = TRUE)
gencode_eval_results <- read.csv("Challenge1_Figures_Data/GENCODE_manualAnnot/new_GENCODE_manualAnnot_evaluation.csv", header = T) # evaluation result
code <- read.csv("Challenge1_Figures_Data/code.csv", header = TRUE)

genocode_eval_WTC11 <- performance.genecode (gencode.pa = pa_GENCODE, ID_UIC = NULL,
                                            pa = pa.WTC11,  code = code, 
                                            selection = NULL, evaluation = gencode_eval_results,
                                            mypattern = "SQ3_human",
                                            directory = "Challenge1_Figures_Data/GENCODE_manualAnnot/classifications/human/")

pivoted_gencode_gene <- pivot_longer(genocode_eval_WTC11, cols = c("Sensitivity.Genes", "Precision.Genes", "F1_score.Genes"))
pivoted_gencode_gene$name <- pivoted_gencode_gene$name %>% factor(levels = c("Sensitivity.Genes", "Precision.Genes", "F1_score.Genes"),
                                                                  labels = c("Sensitivity", "Precision", "F1-score"))

pivoted_gencode_known <- pivot_longer(genocode_eval_WTC11 , cols = c("Sensitivity_known", "Precision_known", "F1_known"))
pivoted_gencode_known$name <- pivoted_gencode_known$name %>% factor(levels = c("Sensitivity_known", "Precision_known", "F1_known"),
                                                                  labels = c("Sensitivity", "Precision", "F1-score"))

pivoted_gencode_novel <- pivot_longer(genocode_eval_WTC11 , cols = c("Sensitivity_novel", "Precision_novel", "F1_novel"))
pivoted_gencode_novel$name <- pivoted_gencode_novel$name %>% factor(levels = c("Sensitivity_novel", "Precision_novel", "F1_novel"),
                                                                  labels = c("Sensitivity", "Precision", "F1-score"))

genocode_eval_WTC11$FALSE_known <- genocode_eval_WTC11$Transcript_models_known - genocode_eval_WTC11$TRUE_known
genocode_eval_WTC11$FALSE_novel <- genocode_eval_WTC11$Transcript_models_novel - genocode_eval_WTC11$TRUE_novel
pivoted_gencode_TP <- pivot_longer(genocode_eval_WTC11 , cols = c("TRUE_known", "TRUE_novel", "FALSE_known", "FALSE_novel"))
pivoted_gencode_TP <- pivoted_gencode_TP[pivoted_gencode_TP$value > 0,]

#pivoted_gencode_TP <- pivoted_gencode_TP %>% separate(name, c("Ground_Truth", "novelty"))
pivoted_gencode_TP$name <- pivoted_gencode_TP$name %>% factor(levels = c("TRUE_known","TRUE_novel",
                                                                         "FALSE_known", "FALSE_novel"),
                                                           labels = c("TRUE\nknown", "TRUE\nnovel",
                                                                      "FALSE\nknown", "FALSE\nnovel"))


#pC.gene.left <- Performance_plot_left(pivoted_gencode_gene, main = "Gene level" )
pC.gene.left <- ggplot(pivoted_gencode_gene, aes(x=Label, y=value)) +
  geom_segment( aes(x=Label, xend=Label, y=0, yend=value, color=Sample_code), size=0.8) +
  geom_point( size=2, aes( shape=Data_Category, color=Sample_code ))  +
  facet_grid( name ~ Alias, scales = "free", space = "free_x", switch = "y"  ) +
  pub_theme +
  scale_color_manual(values = libplat.palette) +
  labs(x="", y="",
       title="Gene level") +
  theme(legend.position="bottom") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 10))+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),limits = c(0, 1), position="right")


#pC.known.mid <- Performance_plot_middel(pivoted_gencode_known, main = "Known_transcript level" )
pC.known.mid <- ggplot(pivoted_gencode_known, aes(x=Label, y=value)) +
  geom_segment( aes(x=Label, xend=Label, y=0, yend=value, color=Sample_code), size=0.8) +
  geom_point( size=2, aes( shape=Data_Category, color=Sample_code ))  +
  facet_grid( name ~ Alias, scales = "free", space = "free_x", switch = "y"  ) +
  pub_theme +
  scale_color_manual(values = libplat.palette) +
  labs(x="", y="",
       title="Known transcript level") +
  theme(legend.position="bottom") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 10))+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),limits = c(0, 1), position="right")

#pC.novel <- Performance_plot_right(pivoted_gencode_novel, main = "Novel_transcript level" )
pC.novel <- ggplot(pivoted_gencode_novel, aes(x=Label, y=value)) +
  geom_segment( aes(x=Label, xend=Label, y=0, yend=value, color=Sample_code), size=0.8) +
  geom_point( size=2, aes( shape=Data_Category, color=Sample_code ))  +
  facet_grid( name ~ Alias, scales = "free", space = "free_x", switch = "y"  ) +
  pub_theme +
  scale_color_manual(values = libplat.palette) +
  labs(x="", y="",
       title="Novel transcript level") +
  theme(legend.position="bottom") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 10))+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),limits = c(0, 1), position="right")

#pC.TP <- Performance_plot_TP(pivoted_gencode_TP , main = "Number detected transcripts" )
pC.TP <- ggplot(pivoted_gencode_TP, aes(x=Label, y=value)) +
  geom_segment( aes(x=Label, xend=Label, y=0, yend=value, color=Sample_code), size=0.8) +
  geom_point( size=2, aes( shape=Data_Category, color=Sample_code ))  +
  facet_grid( name ~ Alias, scales = "free", space = "free_x", switch = "y"  ) +
  pub_theme +
  scale_color_manual(values = libplat.palette) +
  labs(x="", y="",
       title="Detected transcripts of manual curation") +
  theme(legend.position="bottom") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 10))+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),limits = c(0, NA), position="right")

#figure2d <- ggarrange(pC.gene.left, pC.known.mid, pC.novel,
#                      labels = c("", "", ""),
#                      ncol = 3, nrow = 1, common.legend = TRUE)

figure2d <- ggarrange(pC.gene.left, pC.known.mid,
                      pC.novel, pC.TP,
                      ncol = 2, nrow = 2, common.legend = TRUE)

#figure2d2 <- ggarrange(  pC.TP,  pC.TP, labels = c("", ""),
#                         ncol = 2, nrow = 1, common.legend = TRUE)


ggsave(file=paste0(outdir, "/Fig2d_1.svg"), plot=figure2d, width=20, height=7)
#ggsave(file=paste0(outdir, "/Fig2d_2.svg"), plot=figure2d2, width=10, height=7)


###############
### Summary figure

realData_simulation_metrics <- get_realData_simulation_metrics(data_sample = "WTC11",
                                                               outdir = outdir, species = "human",
                                                               sim_directory = "Challenge1_Figures_Data/Simulations/")

summary_GENCODE_eval <- genocode_eval_WTC11 %>% select(Row.names, Sensitivity_known, Precision_known, 
                                                       Sensitivity_novel, Precision_novel, Total_detections)

df_summary_metrics <- merge(realData_simulation_metrics, summary_GENCODE_eval, by="Row.names")

summary_SIRVs <- spliced_SIRV_metrics %>% select(Row.names, Sensitivity, Precision)

df_summary_metrics <- merge(df_summary_metrics, summary_SIRVs, by="Row.names")

df_summary_metrics$perc_3illumina <- df_summary_metrics$perc_3illumina/100
df_summary_metrics$perc_5illumina <- df_summary_metrics$perc_5illumina/100
df_summary_metrics$perc_SRTM <- df_summary_metrics$perc_SRTM/100
df_summary_metrics$perc_SNTM <- df_summary_metrics$perc_SNTM/100
df_summary_metrics$perc_cov <- df_summary_metrics$perc_cov/100

df_pivoted_summary <- pivot_longer(df_summary_metrics, 
                                   cols=c("perc_SRTM", "perc_SNTM", "perc_5illumina", "perc_3illumina","perc_cov",
                                          "Sen_kn","Pre_kn","Sen_no", "Pre_no",
                                          "Sensitivity_known", "Precision_known",
                                          "Sensitivity_novel", "Precision_novel",
                                          "Sensitivity", "Precision"))

df_pivoted_summary <- df_pivoted_summary %>%
  mutate(value=ifelse(value=="NaN",
                      0,
                      value)) %>% 
  group_by(name, Lib_Plat) %>%
  mutate(Ranking = rank(value, na.last = "keep")) %>% 
  mutate(max_rank=max(Ranking, na.rm = T))

df_pivoted_summary$Ranking_adj <- df_pivoted_summary$Ranking*11/df_pivoted_summary$max_rank


df_pivoted_summary$name <-df_pivoted_summary$name %>% 
  factor(levels=rev(c("perc_SRTM","perc_SNTM","perc_5illumina","perc_3illumina","perc_cov",
                     "Sensitivity","Precision",
                     "Sensitivity_known","Precision_known",
                     "Sensitivity_novel","Precision_novel",
                     "Sen_kn", "Pre_kn","Sen_no","Pre_no")),
            labels=rev(c("% SRTM", "% SNTM", "% CAGE-Seq", "% Quant-Seq","% SJ cov.",
                     "SIRV Sensitivity", "SIRV Precision",
                     "GENCODE Sensit. (known)", "GENCODE Prec. (known)",
                     "GENCODE Sensit. (novel)", "GENCODE Prec. (novel)",
                     "SIMULATION Sensit. (known)", "SIMULATION Prec. (known)",
                     "SIMULATION Sensit. (novel)", "SIMULATION Prec. (novel)")))

df_pivoted_summary <- df_pivoted_summary %>% mutate(type_metric=ifelse(
  name %in% c("% SRTM", "% SNTM", "% CAGE-Seq", "% Quant-Seq", "% SJ cov."),
  "Real data",
  ifelse(
    name %in% c("SIRV Sensitivity", "SIRV Precision"),
    "SIRVs",
    ifelse(
      name %in% c("GENCODE Sensit. (known)", "GENCODE Prec. (known)",
                  "GENCODE Sensit. (novel)", "GENCODE Prec. (novel)"),
      "GENCODE manual annot",
      "Simulation"
    )
  )
)) 


top_labels <- c("Rest", "Bronze", "Silver", "Gold")
top_shapes <- c("Rest","Top3","Top2", "Top1")
df_pivoted_summary$Rank_top <- cut(df_pivoted_summary$Ranking_adj,
                                   breaks = c(0,8,9,10,11),
                                   labels = top_shapes,
                                   include.lowest = T)


custom_palette <- colorRampPalette(c("#F5962A", "#ECF52A"))
top <- colorRampPalette(c("#B82D2D", "#F12626"))(3)
bottom <- colorRampPalette(c( "#8AE826", "#5CAE05"))(3)

paleta_final <- c(top, custom_palette(10), bottom)

palette_top <- c("Gold"="#FDDA3A",
                 "Silver"="#DAD8D2",
                 "Bronze"="#CA8B2B",
                 "Rest"="white")
shapes_top <- c("Top1"=24,
                 "Top2"=23,
                 "Top3"=22,
                 "Rest"=21)

pdf(paste0(outdir,"/summary_challenge1_metrics.pdf"), width=9, height = 6)
for (i in c("cDNA-PacBio", "cDNA-ONT",
           "CapTrap-PacBio", "CapTrap-ONT",
           "R2C2-ONT", "dRNA-ONT")){
  
  if (startsWith(i, "cDNA")){
    sections <- c(4.5, 8.5, 10.5)
  }else{
    sections <- c(4.5, 6.5)
  }
  
  p.traffic_Rank <- ggplot(df_pivoted_summary  %>% filter(Lib_Plat==i) %>% na.omit(),
         aes(y=name, x=Data_Category))+
    geom_point(size = 4, stroke=0.4, aes(shape=Rank_top, fill=Ranking_adj)) +  
    geom_hline(yintercept = sections,
               linetype="dashed", color="#0578C2") + 
    facet_grid(Lib_Plat~Alias, drop = T, scales="free")+
    scale_fill_viridis(breaks = c(2, 10), labels = c("Bottom", "Top"))+
    #scale_fill_gradientn(colors = paleta_final, values = rescale(c(1:12)) ,
    #                      na.value = "grey50",
    #                      breaks = c(2, 10), labels = c("Bottom", "Top"))  +
    scale_shape_manual(values = shapes_top)+
    pub_theme+
    theme(axis.text.x = element_text(angle=0, size=8),
          axis.text.y = element_text(size=8),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())+
    labs(colour="Ranking\n")
  
  
  p.traffic_Value <- ggplot(df_pivoted_summary  %>% filter(Lib_Plat==i)%>% na.omit(),
         aes(y=name, x=Data_Category))+
    geom_point(size = 4, stroke=0.4, aes(shape=Rank_top, fill=value)) +  
    geom_hline(yintercept = sections,
               linetype="dashed", color="#0578C2") +
    facet_grid(Lib_Plat~Alias, drop = T, scales="free")+
    scale_fill_viridis(breaks = c(0.1, 0.5, 0.9))  +
    #scale_fill_gradientn(colors = paleta_final, values = rescale(c(0:100)) ,
    #                      na.value = "grey50",
    #                      breaks = c(0.1, 0.5, 0.9))  +
    scale_shape_manual(values=shapes_top)+
    pub_theme+
    theme(axis.text.x = element_text(angle=0, size=8),
          axis.text.y = element_text(size=8),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())+
    labs(colour="Value\n")
  
  p.bar_summary <- ggplot(df_summary_metrics  %>% filter(Lib_Plat==i),
         aes(y=total, x=Data_Category, fill=Tool))+
    geom_bar(stat="identity")+
    geom_text(aes(label=paste(round(total*0.001,digits = 0), "K")),
              vjust=0, size=3) +
    facet_grid(Lib_Plat~Alias, drop = T, scales="free")+
    scale_color_gradientn(colors = paleta_final, values = rescale(c(0:100)) ,
                          na.value = "grey50",
                          breaks = c(0.1, 0.5, 0.9))  +
    scale_fill_met_d("Cross")+
    pub_theme+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none",
          strip.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = -20) ) +
    scale_y_continuous(label = unit_format(unit = "K", scale = 0.001, accuracy = 1),
                       expand = expansion(mult=c(0,0.2)))+
    labs(y="Num. transcripts")
  
  sum1 <- p.bar_summary / p.traffic_Value + 
    plot_layout(heights = c(1, 3))
  sum2 <- p.bar_summary / p.traffic_Rank + 
    plot_layout(heights = c(1, 3))
  
  #ggsave(file=paste0(outdir, "/summary_figure.rank.",i,".svg"), plot=sum2, width=9, height=6)
  #ggsave(file=paste0(outdir, "/summary_figure.value.",i,".svg"), plot=sum1, width=9, height=6)
  
  
  print(sum1)
  print(sum2)
  
}
dev.off()
