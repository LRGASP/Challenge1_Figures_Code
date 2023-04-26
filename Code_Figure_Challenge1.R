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

# PATTERNS
#################

pub_theme <- theme_pubclean(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=13),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=13) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=10), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15.5)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  theme(legend.position = "bottom")

old.libplat.palette = c( "cDNA-PacBio"="#e66adb", "CapTrap-PacBio"="#ab0202", "cDNA-Illumina"="#FFCF71",  "Freestyle-Freestyle"="#75b562",
                         "cDNA-ONT"="#005C75", "CapTrap-ONT"="#7482F0", "R2C2-ONT"="#74CDF0", "dRNA-ONT"="#1b36d1"
)

libplat.palette = c( "cDNA-PacBio"="#c06636", "CapTrap-PacBio"="#802417", "cDNA-Illumina"="#e8b960",  "Freestyle-Freestyle"="#ce9344",
                     "cDNA-ONT"="#646e3b", "CapTrap-ONT"="#17486f", "R2C2-ONT"="#508ea2", "dRNA-ONT"="#2b5851"
)

cat.palette = c( "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "GenicGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                 "Intergenic" = "darksalmon", "GenicIntron"="#41B6C4")

# FUNCTIONS
#################
source("Functions_Supplementary_Figures_Challenge1_v3.R")

setwd("./DataFigures")
data_sample = "WTC11_results/WTC11"

# FIGURES
#################

# Figure 2a. Detections
#########################
FirstPanelChl1 (data_sample = "WTC11", ylims = c(160000, 250000), xlims = c(12000, 28000))
# Plots are generated separately and composed in PowerPoint.

# Figure 2b. Overlap
#########################
agreement.pipelines(data_sample = data_sample)

# Figure 2c. SIRVs
#########################
spliced_SIRV_metrics <- read.csv(paste0(data_sample,".splicedSIRVS_metrics.csv"), sep=",", header=T) %>% t() %>% as.data.frame()
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
figure_name <- paste0(data_sample, ".splicedSIRVs.svg")
ggsave(file=paste0("../Figures_Challenge1/",figure_name), plot=pSIRVspliced, width=9, height=4)

#### same for unspliced

unspliced_SIRV_metrics <- read.csv(paste0(data_sample,".unsplicedSIRVS_metrics.csv"), sep=",", header=T) %>% t()
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

figure_name <- paste0(data_sample, ".unsplicedSIRVs.svg")
ggsave(file=paste0("../Figures_Challenge1/",figure_name), plot=pSIRVunspliced, width=9, height=4)


# Figure 2d. Radarplot
#########################
radar.simulation (species = "human", directory = "Simulations/", pdf = "figures/figure2c")

## Figure 2e. Evaluation against GENCODE
########################################
genocode_eval_WTC11 <- peformance.genecode (gencode.pa = pa_GENCODE, ID_UIC = ID_UIC,
                                            pa = pa.WTC11,  code = code, 
                                            selection = FALSE, evaluation = gencode_eval_results,
                                            mypattern = "SQ3_human",
                                            directory = "GENCODE_manualAnnot/classifications/human/")

pivoted_gencode_gene <- pivot_longer(genocode_eval_WTC11, cols = c("Sensitivity.Genes", "Precision.Genes", "F1_score.Genes"))
pivoted_gencode_known <- pivot_longer(genocode_eval_WTC11 , cols = c("Sensitivity_known", "Precision_known", "F1_known"))
pivoted_gencode_novel <- pivot_longer(genocode_eval_WTC11 , cols = c("Sensitivity_novel", "Precision_novel", "F1_novel"))
genocode_eval_WTC11$False_known <- genocode_eval_WTC11$Transcript_models_known - genocode_eval_WTC11$TRUE_known
genocode_eval_WTC11$False_novel <- genocode_eval_WTC11$Transcript_models_novel - genocode_eval_WTC11$TRUE_novel
pivoted_gencode_TP <- pivot_longer(genocode_eval_WTC11 , cols = c("TRUE_known", "TRUE_novel", "FALSE_known", "FALSE_novel"))
pivoted_gencode_TP <- pivoted_gencode_TP[pivoted_gencode_TP$value > 0,]

pC.gene.left <- Performance_plot_left(pivoted_gencode_gene, main = "Gene level" )
pC.known.mid <- Performance_plot_middel(pivoted_gencode_known, main = "Known_transcript level" )
pC.novel <- Performance_plot_right(pivoted_gencode_novel, main = "Novel_transcript level" )
pC.TP <- Performance_plot_TP(pivoted_gencode_TP , main = "Number detected transcripts" )

figure2d <- ggarrange(pC.gene.left, pC.known.mid,  pC.novel, labels = c("", "", ""),
                      ncol = 3, nrow = 1, common.legend = TRUE)
figure2d2 <- ggarrange(  pC.TP,  pC.TP,labels = c("", ""),
                         ncol = 2, nrow = 1, common.legend = TRUE)

pdf(file= "Manual_Curation_WTC11.pdf", width=10, height=7)
figure2d 
dev.off()

pdf(file= "Manual_Curation_WTC11_2.pdf", width=10, height=7)
figure2d2 
dev.off()

# Local Variables:
# ess-indent-offset: 2
# End:
