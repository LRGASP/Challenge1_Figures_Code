##### Code for Supplementary figures for Challenge 1   ####
###########################################################
traceback_and_stop <- function() {
  traceback(2)
  quit("no", status=1, runLast=FALSE) 
}
if (!interactive()) {
  options(error=traceback_and_stop)
}

## Libraries ##
###############
library(RColorConesa)
library(ggplot2)
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
library(ggrepel)
library(grid)
library(ggbreak) 
library(patchwork)
library(gg.gap)
library(ggcorrplot)
library(viridis)
library(fmsb)

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

## Colors ##
############

palette <- c(cDNA_PacBio = "#c06636", CapTrap_PacBio = "#802417", cDNA_Illumina = "#e8b960", Freestyle_Freestyle = "#ce9344",
             cDNA_ONT = "#646e3b", CapTrap_ONT = "#17486f", R2C2_ONT = "#508Ea2", dRNA_ONT = "#2b5851", Illumina = "#ffcf71",
             ONT = "#2d7Ac0", PacBio = "#d6234a", cDNA = "#9932cc", CapTrap = "#deb887", Freestyle =  "#66cd00", R2C2 = "#ff8c00", dRNA = "#53868b" )

palette1 <- sort(c( CapTrap_ONT = "#17486f", CapTrap_PacBio = "#802417", cDNA_ONT = "#646e3b",cDNA_PacBio = "#c06636", dRNA_ONT = "#2b5851", Freestyle_Freestyle = "#ce9344",
                    R2C2_ONT = "#508Ea2"))

palette2 <- c(ONT = "#2d7Ac0", PacBio = "#d6234a", cDNA = "#9932cc", CapTrap = "#deb887", Freestyle =  "#66cd00", R2C2 = "#ff8c00", dRNA = "#53868b" )

libplat.palette = c( "cDNA-PacBio"="#c06636", "CapTrap-PacBio"="#802417", "cDNA-Illumina"="#e8b960",  "Freestyle-Freestyle"="#ce9344",
                     "cDNA-ONT"="#646e3b", "CapTrap-ONT"="#17486f", "R2C2-ONT"="#508ea2", "dRNA-ONT"="#2b5851"
)

cat.palette = c("FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679",
                "NNC"="#EE6A50", "GenicGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                "Intergenic" = "darksalmon", "GenicIntron"="#41B6C4")

tool.abbreviations = c(Bambu = "Ba", FLAMES = "FM", FLAIR = "FL", IsoQuant = "IQ", 
                       IsoTools = "IT", Iso_IB = "IB", LyRic = "Ly", Mandalorion = "Ma", 
                       TALON_LAPA = "TL", Spectra = "Sp", StringTie2 = "ST")

SC.abbreviations = c(FSM ="FSM", ISM="ISM", NIC="NIC",
                     NNC="NNC", GenicGenomic="Genic", Antisense="AntiS", Fusion="Fusion",
                     Intergenic = "InterG", GenicIntron="Intron")


## Functions ##
###############

replace_elements <- function(A, C_to_B) {
  # get the elements in A that are also present in the first column of C_to_B
  common_elements <- intersect(A, C_to_B$A)
  # replace the common elements with their corresponding B association
  for (i in seq_along(common_elements)) {
    A[A == common_elements[i]] <- C_to_B$B[C_to_B$A == common_elements[i]]
  }
  return(as.vector(A))
}

# Function First panel Challenge 1 main figure.
#############################################

FirstPanelChl1 <- function (data_sample = "H1_mix", ylims = c(160000, 250000), xlims = c(18000, 22000)) {

  working_dir <-paste0("./", data_sample, "_results")
  mydir <- getwd()
  setwd(working_dir)
  
  genes_file <- paste0(data_sample,".genes_SJ_table.csv")
  genes_SJ <- read.csv(genes_file, header = T, sep = ",")%>% t()
  num_trx_file <- paste0(data_sample, ".summary_table_SC.csv") 
  num_trx <- read.csv(num_trx_file, header = T, sep = ",")
  code_file <- paste0(data_sample, ".code_updated.txt")
  code=read.csv(code_file, header = T, sep=",")
  
  code$Lib_Plat <- apply(code, 1, function(x){
    paste(x["Library_Preps"], x["Platform"], sep = "-")
  })
  code$Lib_DC=apply(cbind(code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  code$Label <-apply(cbind(code[,c("Platform","Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  
  
  ##### challenge 1 submission efforts
  ggplot(code, aes(x=Tool, fill=Library_Preps)) +
    geom_bar(stat = "count", position = "stack", color="black") +
    #geom_text(aes(label = ..count.. , group=Lib_Plat), stat = "count", position = position_stack(vjust = 0.5) ,size = 4) +
    pub_theme + 
    scale_fill_conesa(palette = "main", reverse = F) +
    theme(axis.text.x = element_text(angle=30, vjust = 0.7)) +
    ylab("# of submissions") 
  
  ggplot(code, aes(x=Tool, fill=Platform)) +
    geom_bar(stat = "count", position = "stack", color="black") +
    #geom_text(aes(label = ..count.. , group=Lib_Plat), stat = "count", position = position_stack(vjust = 0.5) ,size = 4) +
    pub_theme + 
    scale_fill_conesa(palette = "complete", reverse = F) +
    theme(axis.text.x = element_text(angle=30, vjust = 0.7)) +
    ylab("# of submissions") 
  
  ### try donut plot
  count_LibPlat <- table(code$Lib_Plat) %>% as.data.frame()
  colnames(count_LibPlat) <- c("Lib_Plat", "count")
  count_LibPlat$fraction = count_LibPlat$count / sum(count_LibPlat$count)
  count_LibPlat$ymax = cumsum(count_LibPlat$fraction)
  count_LibPlat$ymin = c(0, head(count_LibPlat$ymax, n=-1))
  count_LibPlat$labelPosition <- (count_LibPlat$ymax + count_LibPlat$ymin) / 2
  count_LibPlat$label <- paste0(count_LibPlat$count)
  
  count_LibPlat <- count_LibPlat %>% 
    arrange(desc(Lib_Plat)) %>%
    mutate(prop = count / sum(count_LibPlat$count) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  ggplot(count_LibPlat, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Lib_Plat)) +
    geom_rect() +
    geom_label( x=3.5, aes(y=labelPosition, label=label), size=6, label.size = NA, show.legend = F) +
    geom_label(aes(x = 2, y = 0, label = "Total\nsubmissions: 47"), size=7, label.size = NA, inherit.aes = FALSE)+
    scale_fill_manual(values = libplat.palette, limits=force, name="") +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "right")
  
  big_df <- merge(num_trx, genes_SJ, by.x="ID",by.y=0)
  
  big_df$FSM_perc <- apply(big_df,1, function(x){
    r <- as.numeric(x["FSM"])*100/as.numeric(x["total"])
    r <- round(r, digits=2)
    r
  })
  
  big_df$ISM_perc <- apply(big_df,1, function(x){
    r <- as.numeric(x["ISM"])*100/as.numeric(x["total"])
    r <- round(r, digits=2)
    r
  })
  
  big_df$NIC_perc <- apply(big_df,1, function(x){
    r <- as.numeric(x["NIC"])*100/as.numeric(x["total"])
    r <- round(r, digits=2)
    r
  })
  
  big_df$NNC_perc <- apply(big_df,1, function(x){
    r <- as.numeric(x["NNC"])*100/as.numeric(x["total"])
    r <- round(r, digits=2)
    r
  })
  
  big_df$knowTrx <- apply(big_df, 1, function(x){
    as.numeric(x["FSM"])+as.numeric(x["ISM"])+as.numeric(x["NIC"])+as.numeric(x["NNC"])
  })
  
  big_df$avgTrxGene <- apply(big_df, 1, function(x){
    as.numeric(x["knowTrx"])/as.numeric(x["Num.KnGenes"])
  })
  
  
  big_df <- merge(big_df, code, by.x="ID", by.y="pipelineCode")
  
  p1.0 <- ggplot(big_df, aes(x=Num.KnGenes, y=total, color=Lib_Plat))+
    geom_point(aes(shape=Data_Category), size=3)+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
    scale_x_continuous(expand = expansion(mult=c(0,0.1)))+
    scale_color_manual(values = libplat.palette) +
    pub_theme  +
    xlab("No. Known Genes") + ylab("No. Transcripts") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_numGenes-numTrx.0.svg")
  ggsave(file=figure_name, plot=p1.0, width=8, height=8)
  
  p1.0_zoom <- p1.0 + 
    scale_y_continuous(label = unit_format(unit = "K", scale = 0.001, accuracy = 1), expand = expansion(mult=c(0,0.1)), limits = ylims)+
    scale_x_continuous(label = unit_format(unit = "K", scale = 0.001, accuracy = 1), expand = expansion(mult=c(0,0.1)), limits = xlims) +
    xlab("") + ylab("")
  
  figure_name <- paste0(data_sample, ".plot_numGenes-numTrx.0_zoom.svg")
  ggsave(file=figure_name, plot=p1.0_zoom, width=4, height=4)
  
  
  p1 <- ggplot(big_df, aes(x=Num.KnGenes, y=total, color=Lib_Plat))+
    geom_point(aes(shape=Data_Category), size=3)+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label = unit_format(unit = "K", scale = 0.001, accuracy = 1), expand = expansion(mult=c(0,0.1)), limits = c(0,100000))+
    scale_x_continuous(label = unit_format(unit = "K", scale = 0.001, accuracy = 1), expand = expansion(mult=c(0,0.1)), limits = c(0,20000))+
    scale_color_manual(values = libplat.palette) +
    pub_theme  +
    xlab("No. Known Genes") + ylab("No. Transcripts") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_numGenes-numTrx.svg")
  ggsave(file=figure_name, plot=p1, width=8, height=8)
  
  
  #+
  #  geom_smooth(method="lm", inherit.aes = F, aes(x=Num.KnGenes, y=knowTrx ), formula = "y ~ x", color="black", alpha=0.7, se=F, size=0.3)   
  
  ## limits added to discard Spectra and IB
  p2 <- ggplot(big_df, aes(x=Num.KnGenes, y=knowTrx, color=Lib_Plat))+
    geom_point(aes(shape=Data_Category), size=3)+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label = unit_format(unit = "K", scale = 0.001, accuracy = 1),expand = expansion(mult=c(0,0.1)), limits = c(0,100000))+
    scale_x_continuous(label = unit_format(unit = "K", scale = 0.001, accuracy = 1),expand = expansion(mult=c(0,0.1)), limits = c(0,20000))+
    scale_color_manual(values = libplat.palette) +
    pub_theme  +
    xlab("No. Known Genes") + ylab("No. Transcripts from Known Genes") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_numGenes-numKnTrx.svg")
  ggsave(file=figure_name, plot=p2, width=8, height=8)
  
  p2_tool <- ggplot(big_df, aes(x=Num.KnGenes, y=knowTrx, color=Tool))+
    geom_point(aes(shape=Data_Category), size=3)+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label = unit_format(unit = "K", scale = 0.001, accuracy = 1),expand = expansion(mult=c(0,0.1)), limits = c(0,100000))+
    scale_x_continuous(label = unit_format(unit = "K", scale = 0.001, accuracy = 1),expand = expansion(mult=c(0,0.1)), limits = c(0,20000))+
    scale_color_conesa(palette = "complete") +
    pub_theme  +
    xlab("No. Known Genes") + ylab("No. Transcripts from Known Genes") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_numGenes-numKnTrx_tool.svg")
  ggsave(file=figure_name, plot=p2_tool, width=8, height=8)
  
  #########plot % FSM vs %ISM
  
  p3 <- ggplot(big_df, aes(x=FSM_perc, y=ISM_perc, color=Lib_Plat))+
    geom_point(aes(size=total, shape=Data_Category))+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label=label_percent(scale = 1), expand = expansion(mult=c(0.1,0.1)))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_color_manual(values = libplat.palette) +
    pub_theme +
    xlab("FSM") + ylab("ISM") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_percFSM-ISM.svg")
  ggsave(file=figure_name, plot=p3, width=8, height=8)
  
  p3_tool <- ggplot(big_df, aes(x=FSM_perc, y=ISM_perc, color=Tool))+
    geom_point(aes(size=total, shape=Data_Category))+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label=label_percent(scale = 1), expand = expansion(mult=c(0.1,0.1)))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_color_conesa(palette = "complete") +
    pub_theme +
    xlab("FSM") + ylab("ISM") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_percFSM-ISM_tool.svg")
  ggsave(file=figure_name, plot=p3_tool, width=8, height=8)
  
  #########plot % NIC vs % NNC
  
  p4 <- ggplot(big_df, aes(x=NIC_perc, y=NNC_perc, color=Lib_Plat))+
    geom_point(aes(size=total, shape=Data_Category))+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_color_manual(values = libplat.palette) +
    pub_theme +
    xlab("NIC") + ylab("NNC") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_percNIC-NNC.svg")
  ggsave(file=,figure_name, plot=p4, width=8, height=8)
  
  p4_tool <- ggplot(big_df, aes(x=NIC_perc, y=NNC_perc, color=Tool))+
    geom_point(aes(size=total, shape=Data_Category))+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_color_conesa(palette = "complete") +
    pub_theme +
    xlab("NIC") + ylab("NNC") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_percNIC-NNC_tool.svg")
  ggsave(file=figure_name, plot=p4_tool, width=8, height=8)
  
  ##############################################
  ###### plots support at the TSS / SJ / TTS
  ###### It will only include info about FSM, ISM, NIC and NNC
  
  ## First we need to read all the data. 
  fsm_metrics <- read.csv(paste0(data_sample,".FSM_metrics.csv"), sep=",", header = T) %>% t()
  ism_metrics <- read.csv(paste0(data_sample,".ISM_metrics.csv"), sep=",", header = T) %>% t()
  nic_metrics <- read.csv(paste0(data_sample,".NIC_metrics.csv"), sep=",", header=T) %>% t()
  nnc_metrics <- read.csv(paste0(data_sample,".NNC_metrics.csv"), sep=",", header=T) %>% t()
  
  support_df <- data.frame(FSM_5ref=fsm_metrics[,"5' reference supported (gene)"], FSM_5illumina=fsm_metrics[,"5' CAGE supported"], 
                           FSM_3ref=fsm_metrics[,"3' reference supported (gene)"], FSM_3illumina=fsm_metrics[,"3' QuantSeq supported"],
                           ISM_5ref=ism_metrics[,"5' reference supported (gene)"], ISM_5illumina=ism_metrics[,"5' CAGE supported"], 
                           ISM_3ref=ism_metrics[,"3' reference supported (gene)"], ISM_3illumina=ism_metrics[,"3' QuantSeq supported"],
                           NIC_5ref=nic_metrics[,"5' reference supported (gene)"], NIC_5illumina=nic_metrics[,"5' CAGE supported"], 
                           NIC_3ref=nic_metrics[,"3' reference supported (gene)"], NIC_3illumina=nic_metrics[,"3' QuantSeq supported"],
                           NNC_5ref=nnc_metrics[,"5' reference supported (gene)"], NNC_5illumina=nnc_metrics[,"5' CAGE supported"], 
                           NNC_3ref=nnc_metrics[,"3' reference supported (gene)"], NNC_3illumina=nnc_metrics[,"3' QuantSeq supported"])
  support_df$total_5ref <- apply(support_df,1, function(x){
    as.numeric(x["FSM_5ref"])+as.numeric(x["ISM_5ref"])+as.numeric(x["NIC_5ref"])+as.numeric(x["NNC_5ref"])
  })
  
  support_df$total_3ref <- apply(support_df,1, function(x){
    as.numeric(x["FSM_3ref"])+as.numeric(x["ISM_3ref"])+as.numeric(x["NIC_3ref"])+as.numeric(x["NNC_3ref"])
  })
  
  support_df$total_5illumina <- apply(support_df,1, function(x){
    as.numeric(x["FSM_5illumina"])+as.numeric(x["ISM_5illumina"])+as.numeric(x["NIC_5illumina"])+as.numeric(x["NNC_5illumina"])
  })
  
  support_df$total_3illumina <- apply(support_df,1, function(x){
    as.numeric(x["FSM_3illumina"])+as.numeric(x["ISM_3illumina"])+as.numeric(x["NIC_3illumina"])+as.numeric(x["NNC_3illumina"])
  })
  
  support_df <- merge(support_df, big_df[,c("ID", "knowTrx")], by.x=0, by.y="ID")
  
  support_df$perc_5ref <- apply(support_df,1,function(x){
    as.numeric(x["total_5ref"])*100/as.numeric(x["knowTrx"])
  })
  
  support_df$perc_5illumina <- apply(support_df,1,function(x){
    as.numeric(x["total_5illumina"])*100/as.numeric(x["knowTrx"])
  })
  
  support_df$perc_3ref <- apply(support_df,1,function(x){
    as.numeric(x["total_3ref"])*100/as.numeric(x["knowTrx"])
  })
  
  support_df$perc_3illumina <- apply(support_df,1,function(x){
    as.numeric(x["total_3illumina"])*100/as.numeric(x["knowTrx"])
  })
  
  support_df <- merge(support_df,code, by.x="Row.names", by.y="pipelineCode")
  
  
  p5 <- ggplot(support_df, aes(x=perc_5ref, y=perc_5illumina, color=Lib_Plat))+
    geom_point(aes(size=knowTrx, shape=Data_Category))+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_color_manual(values = libplat.palette) +
    pub_theme +
    xlab("5' end supported by reference") + ylab("5' end supported by CAGE") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_5endSupport.svg")
  ggsave(file=figure_name, plot=p5, width=8, height=8)
  
  
  p6 <- ggplot(support_df, aes(x=perc_3ref, y=perc_3illumina, color=Lib_Plat))+
    geom_point(aes(size=knowTrx, shape=Data_Category) )+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_color_manual(values = libplat.palette) +
    pub_theme +
    xlab("3' end supported by reference") + ylab("3' end supported by QuantSeq") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_3endSupport.svg")
  ggsave(file=figure_name, plot=p6, width=8, height=8)
  
  #### SJ plots
  
  sj_info <- read.csv(paste0(data_sample,".SJ_info.csv"), sep=",", header=T)
  
  sj_info <- merge(sj_info,code, by="pipelineCode")
  sj_info$perc_known <- apply(sj_info,1,function(x){
    as.numeric(x["known"])*100/as.numeric(x["total"])
  })
  sj_info$perc_canonical <- apply(sj_info,1,function(x){
    as.numeric(x["canonical"])*100/as.numeric(x["total"])
  })
  sj_info$perc_cov <- apply(sj_info,1,function(x){
    as.numeric(x["with_Coverage"])*100/as.numeric(x["total"])
  })
  
  sj_info$perc_RTS <- apply(sj_info,1,function(x){
    as.numeric(x["RTS"])*100/as.numeric(x["total"])
  })
  
  p7 <- ggplot(sj_info, aes(x=perc_known, y=perc_cov, color=Lib_Plat))+
    geom_point(aes(size=total, shape=Data_Category))+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_color_manual(values = libplat.palette) +
    pub_theme +
    xlab("known SJ") + ylab("SJ with Short-Read coverage") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_SJ_known-coverage.svg")
  ggsave(file=figure_name, plot=p7, width=8, height=8)
  
  p8 <- ggplot(sj_info, aes(x=perc_canonical, y=perc_cov, color=Lib_Plat))+
    geom_point(aes(size=total, shape=Data_Category))+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_color_manual(values = libplat.palette) +
    pub_theme +
    xlab("% canonical SJ detected") + ylab("% SJ with Short-Read coverage") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_SJ_canonical-coverage.svg")
  ggsave(file=figure_name, plot=p8, width=8, height=8)
  
  
  #### zoom p8
  p8_zoom <- p8 + xlab("") + ylab("") +
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)), limits = c(92,100))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)), limits = c(98.5,100))
  
  figure_name <- paste0(data_sample, ".plot_SJ_canonical-coverage.zoom.svg")
  ggsave(file=figure_name, plot=p8_zoom, width=8, height=8)
  
  #### plot SRTM vs SNTM
  STM_df <- data.frame(FSM=fsm_metrics[,"Number of isoforms"],
                       FSM_SRTM=fsm_metrics[,"Supported Reference Transcript Model (SRTM)"],
                       ISM=ism_metrics[,"Number of isoforms"],
                       ISM_SRTM=ism_metrics[,"Supported Reference Transcript Model (SRTM)"],
                       NIC=nic_metrics[,"Number of isoforms"],
                       NIC_SNTM=nic_metrics[,"Supported Novel Transcript Model (SNTM)"],
                       NNC=nnc_metrics[,"Number of isoforms"],
                       NNC_SNTM=nnc_metrics[,"Supported Novel Transcript Model (SNTM)"])
  
  STM_df$perc_SRTM <- apply(STM_df,1, function(x){
    r <- (as.numeric(x["FSM_SRTM"])+as.numeric(x["ISM_SRTM"]))*100/(as.numeric(x["FSM"])+as.numeric(x["ISM"]))
    round(r, digits=2)
  })
  STM_df$perc_SNTM <- apply(STM_df,1, function(x){
    r <- (as.numeric(x["NIC_SNTM"])+as.numeric(x["NNC_SNTM"]))*100/(as.numeric(x["NIC"])+as.numeric(x["NNC"]))
    round(r, digits=2)
  })
  
  STM_df$total <- apply(STM_df,1, function(x){
    as.numeric(x["NIC"])+as.numeric(x["NNC"])+as.numeric(x["FSM"])+as.numeric(x["ISM"])
  })
  
  
  STM_df <- merge(STM_df, code, by.x=0, by.y="pipelineCode")
  
  p9 <- ggplot(STM_df, aes(x=perc_SRTM, y=perc_SNTM, color=Lib_Plat))+
    geom_point(aes(size=total, shape=Data_Category))+
    geom_text_repel(aes(label=Alias), size=7, max.overlaps = 100) +
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_x_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_color_manual(values = libplat.palette) +
    pub_theme +
    xlab("% SRTM (only FSM and ISM)") + ylab("% SNTM (only NIC and NNC)") +
    theme( axis.text.y  = element_text( size=18),
           axis.text.x = element_text(size=18),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20)) +
    theme(legend.position = "none")
  
  figure_name <- paste0(data_sample, ".plot_SRTM_SNTM.svg")
  ggsave(file=figure_name, plot=p9, width=8, height=8)
  
  #### get legend image
  legend_plot <- ggplot(support_df, aes(x=Tool, y=perc_3illumina, color=Lib_Plat))+
    geom_point(aes(size=knowTrx, shape=Data_Category))+
    scale_y_continuous(label=label_percent(scale=1), expand = expansion(mult=c(0.1,0.1)))+
    scale_color_manual(values = libplat.palette , limits=c("cDNA-PacBio", "CapTrap-PacBio", "Freestyle-Freestyle", 
                                                           "cDNA-ONT", "CapTrap-ONT", "R2C2-ONT", "dRNA-ONT")) +
    pub_theme +
    theme( legend.key = element_blank(),
           legend.text = element_text(size=18),
           legend.title = element_text(size=20),
           legend.key.size = unit(1, 'cm')) +
    theme(legend.position = "right") + guides(colour = guide_legend(override.aes = list(size=10), 
                                                                    title = "Library Prep - Platform", ncol = 1, order = 1),
                                              size = guide_legend(title="No. Transcripts detected" ),
                                              shape = guide_legend(override.aes = list(size=8), 
                                                                   title = "Data Used"))
  leg <- get_legend(legend_plot)
  
  p_leg <- as_ggplot(leg)
  
  figure_name <- paste0(data_sample, ".plot_legend.svg")
  ggsave(file=figure_name, plot=p_leg, width=5, height=8)
  

  setwd(mydir)
}


##### Functions for general boxplots
number.genes.transcripts <- function (data, code, pa, uic, stats = stats) {
  
  #identify known genes
  pa.uic <- merge (pa, uic, by.x = 1, by.y = 1)
  known.genes <- grep(pattern = "novelGene", x = as.vector(pa.uic$gene), invert = TRUE)
  pa.uic.know <- pa.uic[known.genes,]
  
  # compute number of genes
  data$Number_of_genes <- 0
  for (pipe in data$ID) {
    sel <- pa.uic.know[,pipe] == 1
    data[data$ID == pipe,"Number_of_genes"] <- length(unique(pa.uic.know[sel,"gene"]))
  }
  
  #data2 <- data[,c(2:11,16)]
  #data3 <- apply(data2,2,scale)
  colores <- colorConesa(length(unique(as.factor(data$Tool))), palette = "complete")
  names(colores) <- unique(as.factor(data$Tool))
  data$colores <- colores[data$Tool]
  data$sample <- paste(data$Library_Preps, data$Platform, data$Data_Category, sep = "_")
  data$Sample<- paste(data$Library_Preps, data$Platform, sep = "_")
  data$colorSample <- palette1[data$Sample]
  colnames(data)[colnames(data) == "total"] <- "Number_of_transcripts"
  stats$Sample <- paste(stats$Method, stats$Tech, sep = "_") 
  colnames(stats)[colnames(stats) == "Number_of_supplied_reads"] <- "Number_of_reads"
  
  # add stats data to data
  data <- merge(data, stats, by.x = "Sample", by.y = "Sample")
  
  # compute sd in the number of genes and transcripts
  mad.genes.Tool <- tapply(data$Number_of_genes/1000, data$Tool, mad)
  mad.genes.Platform <- tapply(data$Number_of_genes/1000, data$Platform, mad)
  mad.genes.Library <- tapply(data$Number_of_genes/1000, data$Library_Preps, mad)
  mad.transcripts.Tool <- tapply(data$Number_of_transcripts/1000, data$Tool, mad)
  mad.transcripts.Platform <- tapply(data$Number_of_transcripts/1000, data$Platform, mad)
  mad.transcripts.Library <- tapply(data$Number_of_transcripts/1000, data$Library_Preps, mad)
  
  mad.data.frame <- data.frame(mad_genes = c(mad.genes.Tool, mad.genes.Platform , mad.genes.Library),
                               mad_transcripts = c(mad.transcripts.Tool, mad.transcripts.Platform, mad.transcripts.Library),
                               Grouping = c(rep("Tool", length(mad.genes.Tool)), rep("Platform", length(mad.genes.Platform)), rep ("Library_Prep", length(mad.genes.Library))),
                               Names  = c(names(mad.genes.Tool), names(mad.genes.Platform), names(mad.genes.Library)))
  mad.data.frame$Color  <- c(colores[as.vector(mad.data.frame$Names[1:11])], palette[as.vector(mad.data.frame$Names[12:nrow(mad.data.frame)])]
  )
  result <- list(data.genes = data, mad.data.frame = mad.data.frame)
  result
}

fig.number.genes.transcripts <- function (data, sample, var = "Number_of_reads") {
  A <- ggplot(data, aes_string(y = var, x = "Number_of_transcripts", color = "Sample")) + 
    geom_point()  + 
    scale_color_manual(values=palette1 ) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    ggtitle(sample) +
    theme(plot.title = element_text(hjust = 0.5))
  
  B <- ggplot(data, aes_string(y = var, x = "Number_of_genes", color = "Sample")) + 
    geom_point()  + 
    scale_color_manual(values=palette1 ) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    ggtitle(sample) +
    theme(plot.title = element_text(hjust = 0.5))
  result <- list(A = A, B = B)
  result
  
}

fig.boxplots <- function (mad_data, sample, var.x = "mad_genes", my.xlim = NULL, title.size = 12, xlabel.size = 12,
                          var.y = "Grouping", xlabel = "mad_genes (k)", jitter.size = 2,
                          mycolor = mad_data$Color, jitter.color = "Names", rescale = FALSE) {
  A <- ggplot(mad_data, aes_string(x = var.x , y = var.y)) + 
    geom_boxplot(coef=1e30) +
    geom_jitter(shape=16, position=position_jitter(0.2), size = jitter.size, aes_string(color = jitter.color)) +
    scale_color_manual(values = mycolor) +
    xlab(label = xlabel) +
    ggtitle(sample) +
    theme(axis.title.y=element_blank(),
          axis.title.x = element_text(size = xlabel.size),
          title = element_text(size=title.size),
          plot.title = element_text(hjust = 0.5))
  if (rescale){
    A <- A  + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                            labels = trans_format("log10", math_format(10^.x)))
  }
  if (rescale & !is.null(my.xlim)) {
    A <- A  + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                            labels = trans_format("log10", math_format(10^.x)), limits = my.xlim)
  }
  if (!rescale & !is.null(my.xlim)) {
    A <- A  + xlim(my.xlim)
  }
  A
}


##### Long reads coverage
LRC <- function (directory, sample, main) {
  mydir = getwd()
  setwd(directory)
  files <-  dir( pattern = "*.txt")
  data <- NULL
  for ( i in files) {
    lab <- strsplit(i, split = "_")[[1]][3]
    lab <- strsplit(lab, split = "[.]")[[1]][1]
    lrc <- read.delim (i, as.is = TRUE)
    lrc <- data.frame (lrc, lab = rep(lab, nrow(lrc)))
    data <-rbind(data, lrc)
  }
  setwd(mydir)
  colnames(data) <- c("library", "value", "lab")
  data <- as.data.frame(data)
  data2 <- tapply(data[,2], paste(data[,1], data[,3]), collapse = "_", summary)
  
  data2 <- melt(tapply(as.numeric(data[,2]), list(data[,3],data[,1]), function (x) round(length(which(x > 0.98))/length(x),2)))
  data3 <- melt(tapply(as.numeric(data[,2]), list(data[,3],data[,1]), function (x) round(length(which(x < 0.98 & x > 0.75))/length(x),2)))
  data4 <- melt(tapply(as.numeric(data[,2]), list(data[,3],data[,1]), function (x) round(length(which(x < 0.70))/length(x),2)))
  data5 <- data.frame(Tool = rep(data2[,1],3), pipeline = rep(data2[,2],3),
                      LRC_class = c(rep("a.Over_98%", nrow(data2)), rep("b.From_98%_to_75%", nrow(data2)), rep("c.Under_975%", nrow(data2))),
                      LRC_percentage = c(data2[,3], data3[,3], data4[,3])*100)
  data5 <- data5[!is.na(data5$LRC_percentage),]
  data5[,2] <- sapply(data5[,2], function (x) sub("H1_mix", "H1mix", x))
  data5[,2] <- sapply(data5[,2], function (x) sub("ls", "LS", x))
  data5[,2] <- sapply(data5[,2], function (x) sub("long", "LO", x))
  data5[,2] <- sapply(data5[,2], function (x) sub("ES", "ESmouse", x))
  data5[,1] <- sapply(data5[,1], function (x) sub("FLAIR", "FL", x))
  data5[,1] <- sapply(data5[,1], function (x) sub("LyRic", "Ly", x))
  data5[,1] <- sapply(data5[,1], function (x) sub("Spectra", "Sp", x))
  data5[,1] <- sapply(data5[,1], function (x) sub("TALON", "TL", x))
  data5[,1] <- sapply(data5[,1], function (x) sub("Mandalorion", "Ma", x))
  data5[,1] <- sapply(data5[,1], function (x) sub("JHU", "ST", x))
  data5[,1] <- sapply(data5[,1], function (x) sub("Bambu", "Ba", x))
  data5[,1] <- sapply(data5[,1], function (x) sub("FLAMES", "FM", x))
  data5[,1] <- sapply(data5[,1], function (x) sub("StringTie2", "ST", x))
  data5[,1] <- sapply(data5[,1], function (x) sub("Iso", "IB", x, fixed = TRUE))
  data5[,1] <- sapply(data5[,1], function (x) sub("IBTools", "IT", x, fixed = TRUE))
  data5[,1] <- sapply(data5[,1], function (x) sub("IBQuant", "IQ", x, fixed = TRUE))
  library <- lapply(as.character(data5[,2]), function (x) strsplit(x, split = "_")[[1]])
  data6 <- matrix("NA", ncol = 4 , nrow = length(library))
  data6[,4] <- rep("LO", nrow(data6))
  for (i in 1:nrow(data6)){
    a <- length(library[[i]])
    data6[i,1:a] <- library[[i]]
  }
  data7 <- data.frame(data5, Library = data6[,1], Platform = data6[,2], Sample2 = data6[,3], Data_Category = data6[,4], Lib_DC = paste(data6[,1], data6[,4], sep = "_"))
  data7$Sample <- paste(data7$Library, data7$Platform, sep = "_")
  data7_W <- data7[data7$Sample2 == sample,]
  
  p <- ggplot(data7_W, aes(x= LRC_class, y = LRC_percentage))+
    geom_point() +
    facet_grid(.~ Tool, drop = TRUE) +
    geom_jitter(shape=16, position=position_jitter(0.2), size = 2, aes(color = Sample)) +
    scale_color_manual(values = palette1) +
    ylab("% TM") +
    xlab("LRC_class") +
    ggtitle(main) +
    scale_x_discrete(labels=c("a.Over_98%" = ">98%", "b.From_98%_to_75%" = "98-75%",
                              "c.Under_975%" = "<75%")) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 8),
          axis.title.y = element_text(size = 8),
          title = element_text(size = 8)
    )
  
  p
  
}

##### Biotype
# Creates a plot that identifies the biotype distribution for each pipeline and creates a compartive plot.
# pa: presence_absence file for one sample
# info: file with info for UIC identified within one sample
# code: pipeline codification file
# gtf:biotype annotation file run after prepare script that simplifies the biotype classes
# name: name of the sample
biotype.plot <- function (pa, info, code, gtf, name) {
  code$label <- paste(code$Lab, code$Library_Preps, code$Platform, code$Data_Category, sep = "_")
  sample <- merge(pa, info[,c(1,9)], by.x = 1, by.y = 1)
  ENS <- sapply(sample$gene, function (x)  strsplit(x, split = "_")[[1]])
  sample$ENS <- sapply(ENS, function (x) grep("ENS", x, value = TRUE))
  sample$biotype <- gtf$V3[match(sample$ENS, gtf$V1)]
  novels <-grep ("novelGene", sample$gene)
  sample.known <- sample[setdiff(c(1:nrow(sample)), novels),]
  
  
  pipes <- grep("pipeline", colnames(sample.known), value = TRUE)
  biotypes_by_pipe <-as.data.frame(matrix(0, nrow = length(unique(gtf$V3)), ncol = length(pipes)))
  rownames(biotypes_by_pipe) <- sort(unique(gtf$V3))[c(6,7,2,3,1,4,9,8,5)]
  colnames(biotypes_by_pipe) = code$label
  for (i in 1: length(pipes)) {
    j = code$label[code$pipelineCode == pipes[i]]
    a <- sample.known[,c(pipes[i], "ENS", "biotype")]
    a <- a[a[,1] > 0,]
    b <- unique(a[,c(2,3)])
    c <- table(b[,2])
    biotypes_by_pipe[names(c), j] <- c
  }
  biotypes_by_pipe2<- melt(as.matrix(biotypes_by_pipe))
  colnames(biotypes_by_pipe2) <- c("Biotype", "Pipeline", "Number_Genes")
  
  # Stacked
  p <- ggplot(biotypes_by_pipe2, aes(fill=Biotype, y=Number_Genes, x=Pipeline)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_conesa(palette = "complete") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
    ggtitle (name) +
    theme(plot.title = element_text(hjust = 0.5))
  print (p)
  
  biotypes_by_pipe_rel <- apply(as.matrix(biotypes_by_pipe), 2, function (x) round(x/sum(x),6))
  biotypes_by_pipe_rel2<- melt(biotypes_by_pipe_rel)
  colnames(biotypes_by_pipe_rel2) <- c("Biotype", "Pipeline", "Percentage_Genes")
  
  # Stacked
  q <- ggplot(biotypes_by_pipe_rel2, aes(fill=Biotype, y=Percentage_Genes, x=Pipeline)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_conesa(palette = "complete") +
    ggtitle (name) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
          axis.title.y = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8))
  
  # Boxplot
  h <- ggplot(biotypes_by_pipe_rel2, aes(x=Biotype, y=Percentage_Genes, fill = Biotype)) + 
    geom_boxplot(outlier.colour="black", outlier.shape=16,
                 outlier.size=2, notch=FALSE) +
    scale_fill_conesa(palette = "complete") +
    ggtitle (name) +
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  a <- tapply(biotypes_by_pipe_rel2$Percentage_Genes, biotypes_by_pipe_rel2$Biotype, function (x) round(mean(x)*100,2))
  b <- tapply(biotypes_by_pipe_rel2$Percentage_Genes, biotypes_by_pipe_rel2$Biotype, function (x) round(sd(x)*100,2))
  
  mytable <- as.data.frame(t(data.frame(mean = c("mean",a), sd = c("sd", b))))
  g <- gt(mytable, rowname_col = "V1") %>% tab_header (title = "WTC11 biotype detection stats")
  
  plots <- list(p = p, q = q, h = h, g = g)
  return(plots)
}

##### Simplify Biotypes in GTF
# Function to reduce the diversity of biotypes
simplify.biotypes <- function (gtf) {
  gtf$V3 <- gtf$V2
  gtf$V3[grep("IG_", gtf$V2)] <- "IG_gene"
  gtf$V3[grep("TC_", gtf$V2)] <- "TR_gene"
  gtf$V3[grep("TR_", gtf$V2)] <- "TR_gene"
  gtf$V3[grep("Mt_", gtf$V2)] <- "Mt_gene"
  gtf$V3[grep("pseudogene", gtf$V2)] <- "pseudogene"
  gtf$V3[grep("snRNA", gtf$V2)] <- "other_RNA"
  gtf$V3[grep("misc_RNA", gtf$V2)] <- "other_RNA"
  gtf$V3[grep("TEC", gtf$V2)] <- "other_RNA"
  gtf$V3[grep("snoRNA", gtf$V2)] <- "other_RNA"
  gtf$V3[grep("scaRNA", gtf$V2)] <- "other_RNA"
  gtf$V3[grep("scRNA", gtf$V2)] <- "other_RNA"
  gtf$V3[grep("vault_RNA", gtf$V2)] <- "other_RNA"
  gtf$V3[grep("sRNA", gtf$V2)] <- "other_RNA"
  gtf$V3[grep("ribozyme", gtf$V2)] <- "rRNA"
  unique(gtf$V3)
  gtf
}

# me_and_others function evaluates the number of pipelies that agree with
# the transcript models of one tool for a given sample. The function also returns the agreement 
# of different samples analyzed with the same tool. Uses to.table function
#
to.table <- function (SC.distr, col1 = "Number_of_tools" ) {
  cols <- unique(as.vector(unlist(sapply(SC.distr, names))))
  tab <- matrix(0, nrow = length(SC.distr), ncol = length(cols))
  rownames(tab) <- names(SC.distr)
  colnames(tab) <- cols
  for ( i in 1: nrow(tab)) {
    tab[i,] <- SC.distr[[i]][colnames(tab)]
  }
  tab[is.na(tab)] <- 0
  tab <- melt(tab)
  colnames(tab) <- c(col1, "Structural_category", "Number_of_UIC")
  tab
}
me_and_others <- function (pa, code, subset= c("PacBio", "cDNA"), name = "WTC11",remove = c("Iso_IB","Spectra"), replace = FALSE) {
  if (replace) {
    code$Lab <- tool.abbreviations[match(code$Lab, names(tool.abbreviations))]
    remove <- tool.abbreviations[match(remove, names(tool.abbreviations))]
  }
  pa.pipe <- pa[,code[,1]] 
  all.labs <- unique(code$Lab)
  Lab.colors = colorConesa(n = length(all.labs)) ; names (Lab.colors) = all.labs
  if (!is.null(subset)) { # opens bracket for across tools analysis
    sel.subset <- code$Library_Preps == subset[2] & code$Platform == subset[1] 
    pa.pipe <- pa.pipe[,sel.subset]
    code2 <- code[sel.subset,]
    pa.pipe2 <- t(apply(pa.pipe, 1, function (x) tapply (x, code2$Lab, max)))
    rownames(pa.pipe2) <- pa[,1]
    Labs <- unique(code2$Lab)
    if (!is.null(remove) & any(is.element(remove, Labs))) {
      Labs <- setdiff(Labs,remove)
      pa.pipe2 <- pa.pipe2[,setdiff(colnames(pa.pipe2),remove)] 
    }
    SC.distr2 <- NULL
    for (i in c(1:length(Labs))) {
      sel.lab <- which(pa.pipe2[,Labs[i]] == 1) # UIC positions detected by the Lab pipeline
      which.col = which(colnames(pa.pipe2) == Labs[i])
      rS <- rowSums(pa.pipe2[sel.lab,- which.col])
      SC.distr <- tapply(pa$structural_category[sel.lab], rS , table)
      table.SC.distr <- to.table(SC.distr)
      SC.distr2 <- rbind(SC.distr2, data.frame(table.SC.distr, Lab = rep(Labs[i],nrow(table.SC.distr) )))
    }
    h <- ggplot(data=SC.distr2 , aes(x = Number_of_tools, y = Number_of_UIC, group = Structural_category, colour = Structural_category)) +
      geom_line() +
      geom_point(size = 0.5) +
      scale_y_continuous(trans='log10') +
      ggtitle(paste(c(subset,name), collapse = "_")) + 
      theme(axis.text.x = element_text(size = 8),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(hjust = 0.5, size = 10),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8)) +
      scale_color_manual(name = "Structural_category", values = cat.palette , 
                         limits = names(cat.palette)) +
      facet_grid(.~ Lab, drop = TRUE)
  } # opens option for within tool analysis
  code2 <- code
  pa.pipe <- pa[,code[,1]] 
  subset <- "Within_tool"
  index <- paste(code2$Lab, code2$Library_Preps, code2$Platform, sep = "_")
  pa.pipe2 <-  t(apply(pa.pipe, 1, function (x) tapply (x, index, max)))
  annot <- sapply(colnames(pa.pipe2), strsplit, split = "_")
  annot.lab <- sapply(annot, function (x) x[1])
  Labs <- unique(annot.lab)
  if (!is.null(remove) & any(is.element(remove, Labs))) {
    Labs <- setdiff(Labs,remove)
    pa.pipe3 <- pa.pipe2[,is.element(annot.lab,Labs)]
    annot.lab <- annot.lab[is.element(annot.lab,Labs)]
  }
  SC.distr3 = NULL
  for (i in c(1:length(Labs))) {
    if(length(which(annot.lab == Labs[i])) > 1) {
      sel.lab <- pa.pipe3[,annot.lab == Labs[i]] # selects the columns of the tool
      sel.UIC <- rowSums(sel.lab) > 0 # selects the UIC detected by the tool in at least one pipeline
      rS <- rowSums(sel.lab[sel.UIC,]) # how many times the UIC was detected by the tool
      SC.distr <- tapply(pa$structural_category[sel.UIC], rS , table) 
      SC.distr_b <- to.table(SC.distr, col1 = "Number_of_samples")
      SC.distr3 <- rbind(SC.distr3, data.frame(SC.distr_b, Tool = rep(Labs[i],nrow(SC.distr_b) )))
    }
  }
  j <- ggplot(data=SC.distr3 , aes(x = Number_of_samples, y = Number_of_UIC, group = Structural_category, colour = Structural_category)) +
    geom_line() +
    geom_point(size = 0.5) +
    scale_y_continuous(trans='log10') +
    ggtitle(paste(c(subset,name), collapse = "_")) + 
    theme(axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8)) +
    scale_color_manual(name = "Structural_category", values = cat.palette , 
                       limits = names(cat.palette)) +
    facet_grid(.~ Tool, drop = TRUE)
  
  return(list(Overlap = SC.distr2, h = h, j = j))
}

#####  Analysis of UIC that were found by multiple experimental datasets (num.samples) and analysis tools (num.tools). Uses enrichment function
enrichment <- function (a, b, round = 1) {
  per.a <- a/sum(a)
  per.b <- b/sum(b)
  per.c <- data.frame (c = rep(0, length(per.b)), b = as.vector(per.b))
  rownames(per.c) <- names(per.b)
  per.c[names(per.a),1] <-per.a
  per.c$Enrichment <- round(per.c$c/ per.c$b,round)
  per.c
  
}
highly.detected <- function (pa, code, num.samples = 3, num.tools = 3, name = "WTC11", replace = TRUE){
  if (replace) {
    code$Lab <- tool.abbreviations[match(code$Lab, names(tool.abbreviations))]
    pa$structural_category <- SC.abbreviations[match(pa$structural_category, names(SC.abbreviations))]
    names(cat.palette)<- SC.abbreviations[match(names(cat.palette), names(SC.abbreviations))]
    
  }
  pa.pipe <- pa[,code[,1]] 
  code$Pipeline <- paste(code$Library_Preps, code$Platform, code$Lab, sep = "_")
  pa.pipe2 <-  t(apply(pa.pipe, 1, function (x) tapply (x, code$Pipeline, max)))
  annot <- sapply(colnames(pa.pipe2), strsplit, split = "_")
  annot.lab <- sapply(annot, function (x) x[3])
  annot.prep <- sapply(annot, function (x) x[1])
  annot.plat <- sapply(annot, function (x) x[2])
  sample <- paste(annot.prep, annot.plat, sep = "_")
  count.lab <- t(apply(pa.pipe2, 1, function (x) tapply(x, annot.lab, sum)))
  count.sample <- t(apply(pa.pipe2, 1, function (x) tapply(x, sample, sum)))
  condition.lab <- apply(count.lab,1, function (x) any(x >= num.tools))
  condition.sample <- apply(count.sample,1, function (x) any(x >= num.samples))
  sel.UIC <- pa[(condition.lab & condition.sample) ,]
  pa.pipe3 <- pa.pipe2[(condition.lab & condition.sample),]
  hc.sc <- table(sel.UIC$structural_category)
  sc.base <- table(pa$structural_category)
  enrichment.sc <- enrichment(hc.sc, sc.base)
  enrichment.sc$Structural_category <- rownames(enrichment.sc)
  enrichment.sc$count <- hc.sc[rownames(enrichment.sc)]
  hc.sc.df <- data.frame(Structural_category = names(hc.sc), Count = as.vector(hc.sc))
  
  hc.lab <- table(unlist(apply(pa.pipe3 ,1, function (x) annot.lab[x!=0])))
  enrichment.lab <- enrichment(hc.lab , table(unlist(apply(pa.pipe2 ,1, function (x) annot.lab[x!=0]))))
  enrichment.lab$Tool <- rownames(enrichment.lab)
  enrichment.lab$Count <- hc.lab[rownames(enrichment.lab)] 
  
  hc.plat <- table(unlist(apply(pa.pipe3 ,1, function (x) annot.plat[x!=0])))
  enrichment.plat <- enrichment(hc.plat, table(unlist(apply(pa.pipe2 ,1, function (x) annot.plat[x!=0]))))
  hc.prep <- table(unlist(apply(pa.pipe3 ,1, function (x) annot.prep[x!=0])))
  enrichment.prep <- enrichment(hc.prep, table(unlist(apply(pa.pipe2 ,1, function (x) annot.prep[x!=0]))))
  
  a <- ggplot(hc.sc.df, aes(x=Structural_category, y=Count, fill=Structural_category)) +
    geom_bar(stat="identity", width=1, color="white", show.legend = FALSE) +
    theme_bw() +
    ggtitle(paste("Frequently detected UICs \n ", name, sep = "")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_discrete(limits = names(cat.palette)) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
          axis.title.x = element_text(color = "white")) +
    scale_fill_manual(values = cat.palette) +
    annotate(geom = "table",
             x = 9,
             y = 40000,
             label = list(enrichment.sc[names(cat.palette),c(4,3)]))
  
  
  
  b <- ggplot(enrichment.sc, aes(x = Structural_category, y = Enrichment, fill = Structural_category)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_hline(yintercept=1) +
    ggtitle(paste("Structural category enrichment within frequently detected UICs \n Sample ", name, sep = "")) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_discrete(limits = names(cat.palette)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
          axis.title.x = element_text(color = "white")) +
    scale_fill_manual(values = cat.palette)
  
  tool.colors = colorConesa(n = nrow(enrichment.lab)) ; names (tool.colors) = rownames(enrichment.lab)
  c <- ggplot(enrichment.lab, aes(x=Tool, y=Count, fill=Tool)) +
    geom_bar(stat="identity", width=1, color="white", show.legend = FALSE) +
    theme_bw() +
    ggtitle(paste("Tools count for frequently detected UICs in sample ", name, sep = "")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
          axis.title.x = element_text(color = "white")) +
    scale_fill_manual(values = tool.colors)
  
  
  mytable <- as.data.frame(t(enrichment.lab[,c(4,5)]))[2,]
  d <- ggplot(enrichment.lab, aes(x = Tool, y = Enrichment, fill = Tool)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_hline(yintercept=1) +
    ggtitle(paste("Tool enrichment within frequently detected UICs \n ", name, sep = "")) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    scale_fill_manual(values = tool.colors) +
    annotate(geom = "table",
             x = 6,
             y = 2,
             cex= 3,
             label = list(mytable))
  
  result <- list(a = a, b = b, c = c, d = d, 
                 enrichment.lab = enrichment.lab, hc.sc.df = hc.sc.df, enrichment.sc = enrichment.sc,
                 pa.pipe3 = pa.pipe3, pa.pipe2 = pa.pipe2, annot.plat = annot.plat)
  result
}

##### Function to calculate the coverage of simulated long reads of the targeted transcripts
line.plot <- function (df , title) {
  p <- ggplot(data=df, aes(x=transcript_position, y=coverage, fill=type_of_data)) +
    geom_line(aes(color=type_of_data))+
    scale_color_manual(values=c("#17486f", "#c06636")) +
    theme_classic() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    labs(title= title,
         x ="Transcript position, %") +
    theme(plot.title = element_text(hjust = 0.5))
  p
}

##### Function to evaluate the differences in TM length, number of exons and CPM between LRGASP experimental methods
trend.analysis <- function (pa, code, info) {
  pa <- merge (pa, info[,c("LRGASP_id", "exons",  "length", "FL_cpm")], by.x = 1, by.y = 1)
  ONT <- code$pipelineCode[code$Platform == "ONT"]
  PacBio <- code$pipelineCode[code$Platform == "PacBio"]
  cDNA <- code$pipelineCode[code$Library_Preps == "cDNA"]
  CapTrap <- code$pipelineCode[code$Library_Preps == "CapTrap"]
  dRNA <- code$pipelineCode[code$Library_Preps == "dRNA"]
  R2C2<- code$pipelineCode[code$Library_Preps == "R2C2"]
  cDNA.PacBio <- intersect(PacBio, cDNA)
  cDNA.ONT <- intersect(ONT, cDNA)
  CapTrap.PacBio <- intersect(PacBio, CapTrap)
  CapTrap.ONT <- intersect(ONT, CapTrap)
  
  ONT.rs <- rowSums(pa[,ONT])
  PacBio.rs <- rowSums(pa[,PacBio])
  cDNA.rs <- rowSums(pa[,cDNA])
  CapTrap.rs <- rowSums(pa[,CapTrap])
  dRNA.rs <- rowSums(pa[,dRNA])
  R2C2.rs <- rowSums(pa[,R2C2])
  cDNA.PacBio.rs <- rowSums(pa[,cDNA.PacBio])
  cDNA.ONT.rs <- rowSums(pa[,cDNA.ONT])
  CapTrap.PacBio.rs <- rowSums(pa[,CapTrap.PacBio])
  CapTrap.ONT.rs <- rowSums(pa[,CapTrap.ONT])
  
  onlyONT <- ONT.rs > 0 & PacBio.rs == 0 ; a <- length(which (onlyONT))
  onlyPacBio <- ONT.rs == 0 & PacBio.rs >= 0 ; b <- length(which (onlyPacBio))
  both <- ONT.rs >0 & PacBio.rs > 0 ; c<- length(which(both))
  
  onlycDNA <- cDNA.rs > 0 & CapTrap.rs == 0  & dRNA.rs == 0 & R2C2.rs ==0 ; only.c <- length(which (onlycDNA))
  onlyCapTrap <- CapTrap.rs > 0 & cDNA.rs == 0  & dRNA.rs == 0 &R2C2.rs ==0 ; only.t <- length(which (onlyCapTrap))
  onlydRNA <- dRNA.rs > 0 & cDNA.rs == 0  & CapTrap.rs == 0 & R2C2.rs ==0 ; only.d <- length(which (onlydRNA))
  onlyR2C2 <- R2C2.rs  > 0 & cDNA.rs == 0  & CapTrap.rs == 0 & dRNA.rs ==0 ; only.r <- length(which (onlyR2C2))
  
  onlycDNA.PacBio <- cDNA.PacBio.rs  > 0 & cDNA.ONT.rs== 0   ; only.cb <- length(which (onlycDNA.PacBio))
  onlycDNA.ONT <- cDNA.ONT.rs > 0 & cDNA.PacBio.rs == 0 ; only.cn <- length(which (onlycDNA.ONT))
  onlyCapTrap.PacBio <- CapTrap.PacBio.rs > 0 & CapTrap.ONT.rs == 0 ;  only.db <- length(which(onlyCapTrap.PacBio))
  onlyCapTrap.ONT <- CapTrap.ONT.rs  > 0 & CapTrap.PacBio.rs == 0 ; only.rt <- length(which(onlyCapTrap.ONT ))
  
  
  bp.data <- data.frame(Length = c(pa[onlyONT,"length"],pa[onlyPacBio,"length"],pa[both,"length"]),
                        Exons = c(pa[onlyONT,"exons"],pa[onlyPacBio,"exons"],pa[both,"exons"]),
                        CPM = c(pa[onlyONT,"FL_cpm"],pa[onlyPacBio,"FL_cpm"],pa[both,"FL_cpm"]),
                        Detected_by = c(rep("onlyONT",a), rep("onlyPacBio",b), rep ("both", c)))
  
  
  bp.data2 <- data.frame(Length = c(pa[onlycDNA,"length"],pa[onlyCapTrap,"length"],pa[onlydRNA,"length"],pa[onlyR2C2,"length"]),
                         Exons = c(pa[onlycDNA,"exons"],pa[onlyCapTrap,"exons"],pa[onlydRNA,"exons"],pa[onlyR2C2,"exons"]),
                         CPM = c(pa[onlycDNA,"FL_cpm"],pa[onlyCapTrap,"FL_cpm"],pa[onlydRNA,"FL_cpm"],pa[onlyR2C2,"FL_cpm"]),
                         Detected_by = c(rep("onlycDNA",only.c), rep("onlyCapTrap",only.t), rep ("onlydRNA", only.d), rep ("onlyR2C2", only.r)))
  
  
  bp.data3 <- data.frame(Length = c(pa[onlycDNA.PacBio,"length"],pa[onlycDNA.ONT,"length"],pa[onlyCapTrap.PacBio,"length"],pa[onlyCapTrap.ONT,"length"]),
                         Exons = c(pa[onlycDNA.PacBio,"exons"],pa[onlycDNA.ONT,"exons"],pa[onlyCapTrap.PacBio,"exons"],pa[onlyCapTrap.ONT,"exons"]),
                         CPM = c(pa[onlycDNA.PacBio,"FL_cpm"],pa[onlycDNA.ONT,"FL_cpm"],pa[onlyCapTrap.PacBio,"FL_cpm"],pa[onlyCapTrap.ONT,"FL_cpm"]),
                         Detected_by = c(rep("onlycDNA.PacBio",only.cb), rep("onlycDNA.ONT",only.cn), rep ("onlyCapTrap.PacBio", only.db), rep ("onlyCapTrap.ONT", only.rt)))
  
  result <- list(bp.data = bp.data, bp.data2 = bp.data2, bp.data3 = bp.data3)
  result
}

##### Function to evaluate the differences in TM length across analysis methods
length_pipelines <- function(pa, uic, name ){
  data.41 <- merge(pa, uic[,c("LRGASP_id", "length")], by.x = "TAGS", by.y = "LRGASP_id")
  lengths  = NULL ;    pipes = NULL
  for (h in 3: (ncol(data.41)-1) ){
    l <- data.41[data.41[,h] == 1,"length"]
    lengths = c(lengths, l)
    pipes = c(pipes, rep(colnames(data.41)[h], length(l)))
  }
  code$Lib_Plat <- paste(code$Library_Preps, code$Platform, sep = "_")
  data_lengths <- data.frame(length = lengths, pipeline = pipes)
  data_lengths <- merge(data_lengths , code, by.x = "pipeline", by.y = "pipelineCode" )
  
  bp <- ggplot(data_lengths, aes(x=pipeline, y=length, group=pipeline, fill = Lib_Plat )) + 
    geom_boxplot(outlier.size=0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(name) +
    scale_y_continuous(trans='log10') +
    scale_fill_manual(values= palette1) +
    facet_grid( .~ Alias, scales = "free_x", space="free")
  bp
}
comparisons.plot <- function (bp.data, my_comparisons, name) {
  P <- vector(mode = "list", length = 3)
  features <- c("Length" ,"Exons", "CPM")
  for ( i in 1:length(features)) {
    p <- ggplot(bp.data, aes_string(x = "Detected_by", y = features[i],  fill = "Detected_by")) +
      geom_boxplot() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(name) +
      scale_y_continuous(trans='log10') +
      scale_fill_manual(values= colorConesa(length(unique(bp.data$Detected_by)), palette = "main"))
    if (!is.null (my_comparisons)) {
      p <- p + stat_compare_means(comparisons = my_comparisons)
    }
    P[[i]] <- p
  }
  P
}

##### Function to display the agreement in transcript detections across pipelines
agreement.pipelines <- function (data_sample) {
  pa_table <- read.csv(paste0("./",data_sample,"_comparison.pa.csv"), sep=",", header = T)
  pa_table <- pa_table[!grepl("SIRV", pa_table$TAGS),]
  pa_table <- pa_table[!grepl("ERCC", pa_table$TAGS),]
  
  pa_table$found_by <- apply(pa_table,1,function(x){
    sum(as.numeric(x[3:49]))
  })
  
  p10.1 <- ggplot(pa_table, aes(x=found_by)) +
    geom_bar() + 
    pub_theme +
    labs(y="Number UIC", x="") +
    scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001, accuracy = 1),  expand = expansion(mult=c(0,0.1)))+
    scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
    theme(plot.margin=unit(c(0.5,0.5,-0.5,0.5), "cm"),
          axis.title.y = element_text(size = 10))
  
  p10.1_2 <- ggplot(pa_table, aes(x=found_by)) +
    geom_bar() + 
    pub_theme +
    labs(y="Number of UIC, log10", x="") +
    scale_y_continuous(expand = expansion(mult=c(0,0.1)),
                       trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))+
    scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
    theme(plot.margin=unit(c(0.5,0.5,-0.5,0.5), "cm"),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 10))
  
  p10.2 <- ggplot(pa_table, aes(x=found_by)) +
    geom_bar(position = "fill", aes(fill=structural_category)) + 
    pub_theme +
    labs(x="Number of submissions that found a UIC", y="Structural Category\ndistribution")+
    scale_fill_manual(values=cat.palette, name="" ) +
    scale_y_continuous(labels = label_percent(), expand = expansion(mult=c(0,0))) +
    scale_x_continuous(expand = expansion(mult=c(0.01,0.01)))+
    theme(plot.margin=unit(c(-0.5,0.5,0,0.5), "cm"), 
          axis.title.y = element_text(size = 10)) +
    guides(fill = guide_legend( nrow = 1, order = 1))
  
  p10 <- p10.1_2 / p10.2 + 
    plot_layout(heights = c(2, 1))
  
  p10
}
##### Functions to compute the pairwise overlap between pipelies
get_overlap_matrix <- function(pa, pipe){
  only_pipelines <- pa[, pipe]
  pip <- colnames(only_pipelines)
  id <- expand.grid(pipe,pipe)
  pairwise_intersection_matrix <- matrix( colSums( only_pipelines[ , id[,1] ] == only_pipelines[ , id[,2] ] & only_pipelines[ , id[,1] ]!=0) , ncol = length(pip) )
  pairwise_ov1_matrix <-  matrix( colSums( (only_pipelines[ , id[,1] ]==1) ),
                                  ncol = length(pip) )
  
  p2 <- pairwise_intersection_matrix/pairwise_ov1_matrix
  overlap_matrix <- p2
  colnames(overlap_matrix) <- pipe
  rownames(overlap_matrix) <- pipe
  
  return(overlap_matrix)
}
overlap.plot <- function (code, data_sample){
  
  code$pipeline_name <- apply(code,1,function(x){paste0(c(x["Library_Preps"], x["Platform"], x["Alias"], x["Data_Category"]), collapse="-")})
  
  ### read pa file and 
  pa_table = read.csv(paste0(data_sample,".pa.csv"), sep = ",", header = T )
  pips <- colnames(pa_table)[colnames(pa_table) %in% code$pipelineCode]
  pips<-merge(as.data.frame(pips), code, by.y="pipelineCode", by.x="pips", sort=F)
  
  colnames(pa_table)[colnames(pa_table) %in% code$pipelineCode] <- pips$pipeline_name
  
  overlap_mat_UJC <- get_overlap_matrix(pa_table, as.character(code$pipeline_name))
  melted_overlap <- melt(overlap_mat_UJC)
  pOverlap <- ggplot(melted_overlap, aes(Var1, Var2, fill= value)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE, begin = 0, end=1, limit = c(0,1), 
                       name="Overlap\nindex") +
    pub_theme +
    theme(axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(angle=0, size=10), 
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          legend.position = "right") +
    coord_flip()
  pOverlap
}

#####  Radar plots showing performance metrics of LRGASP tools on simlulated data

radarchart2 <- function (df, axistype = 0, seg = 4, pty = 16, pcol = 1:8, plty = 1:6, 
                         plwd = 1, pdensity = NULL, pangle = 45, pfcol = NA, cglty = 3, 
                         cglwd = 1, cglcol = "navy", axislabcol = "blue", vlabcol = "black", title = "", 
                         maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlabels = NULL, 
                         vlcex = NULL, caxislabels = NULL, calcex = NULL, paxislabels = NULL, 
                         palcex = NULL, ...) 
{
  if (!is.data.frame(df)) {
    cat("The data must be given as dataframe.\n")
    return()
  }
  if ((n <- length(df)) < 3) {
    cat("The number of variables must be 3 or more.\n")
    return()
  }
  if (maxmin == FALSE) {
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  plot(c(-1.2, 1.2), c(-1.2, 1.2), type = "n", frame.plot = FALSE, 
       axes = FALSE, xlab = "", ylab = "", main = title, asp = 1, 
       ...)
  theta <- seq(90, 450, length = n + 1) * pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) {
    polygon(xx * (i + CGap)/(seg + CGap), yy * (i + CGap)/(seg + 
                                                             CGap), lty = cglty, lwd = cglwd, border = cglcol)
    if (axistype == 1 | axistype == 3) 
      CAXISLABELS <- paste(i/seg * 100, "(%)")
    if (axistype == 4 | axistype == 5) 
      CAXISLABELS <- sprintf("%3.2f", i/seg)
    if (!is.null(caxislabels) & (i < length(caxislabels))) 
      CAXISLABELS <- caxislabels[i + 1]
    if (axistype == 1 | axistype == 3 | axistype == 4 | 
        axistype == 5) {
      if (is.null(calcex)) 
        text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS, 
             col = axislabcol)
      else text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS, 
                col = axislabcol, cex = calcex)
    }
  }
  if (centerzero) {
    arrows(0, 0, xx * 1, yy * 1, lwd = cglwd, lty = cglty, 
           length = 0, col = cglcol)
  }
  else {
    arrows(xx/(seg + CGap), yy/(seg + CGap), xx * 1, yy * 
             1, lwd = cglwd, lty = cglty, length = 0, col = cglcol)
  }
  PAXISLABELS <- df[1, 1:n]
  if (!is.null(paxislabels)) 
    PAXISLABELS <- paxislabels
  if (axistype == 2 | axistype == 3 | axistype == 5) {
    if (is.null(palcex)) 
      text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol)
    else text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol, 
              cex = palcex)
  }
  VLABELS <- colnames(df)
  if (!is.null(vlabels)) 
    VLABELS <- vlabels
  if (is.null(vlcex)) 
    text(xx * 1.2, yy * 1.2, VLABELS, col = vlabcol)
  else text(xx * 1.2, yy * 1.2, VLABELS, cex = vlcex, col = vlabcol)
  series <- length(df[[1]])
  SX <- series - 2
  if (length(pty) < SX) {
    ptys <- rep(pty, SX)
  }
  else {
    ptys <- pty
  }
  if (length(pcol) < SX) {
    pcols <- rep(pcol, SX)
  }
  else {
    pcols <- pcol
  }
  if (length(plty) < SX) {
    pltys <- rep(plty, SX)
  }
  else {
    pltys <- plty
  }
  if (length(plwd) < SX) {
    plwds <- rep(plwd, SX)
  }
  else {
    plwds <- plwd
  }
  if (length(pdensity) < SX) {
    pdensities <- rep(pdensity, SX)
  }
  else {
    pdensities <- pdensity
  }
  if (length(pangle) < SX) {
    pangles <- rep(pangle, SX)
  }
  else {
    pangles <- pangle
  }
  if (length(pfcol) < SX) {
    pfcols <- rep(pfcol, SX)
  }
  else {
    pfcols <- pfcol
  }
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    scale <- CGap/(seg + CGap) + (df[i, ] - df[2, ])/(df[1, 
                                                         ] - df[2, ]) * seg/(seg + CGap)
    if (sum(!is.na(df[i, ])) < 3) {
      cat(sprintf("[DATA NOT ENOUGH] at %d\n%g\n", i, 
                  df[i, ]))
    }
    else {
      for (j in 1:n) {
        if (is.na(df[i, j])) {
          if (na.itp) {
            left <- ifelse(j > 1, j - 1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left > 1, left - 1, n)
            }
            right <- ifelse(j < n, j + 1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right < n, right + 1, 
                              1)
            }
            xxleft <- xx[left] * CGap/(seg + CGap) + 
              xx[left] * (df[i, left] - df[2, left])/(df[1, 
                                                         left] - df[2, left]) * seg/(seg + CGap)
            yyleft <- yy[left] * CGap/(seg + CGap) + 
              yy[left] * (df[i, left] - df[2, left])/(df[1, 
                                                         left] - df[2, left]) * seg/(seg + CGap)
            xxright <- xx[right] * CGap/(seg + CGap) + 
              xx[right] * (df[i, right] - df[2, right])/(df[1, 
                                                            right] - df[2, right]) * seg/(seg + 
                                                                                            CGap)
            yyright <- yy[right] * CGap/(seg + CGap) + 
              yy[right] * (df[i, right] - df[2, right])/(df[1, 
                                                            right] - df[2, right]) * seg/(seg + 
                                                                                            CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft
              yytmp <- yyleft
              xxleft <- xxright
              yyleft <- yyright
              xxright <- xxtmp
              yyright <- yytmp
            }
            xxs[j] <- xx[j] * (yyleft * xxright - yyright * 
                                 xxleft)/(yy[j] * (xxright - xxleft) - 
                                            xx[j] * (yyright - yyleft))
            yys[j] <- (yy[j]/xx[j]) * xxs[j]
          }
          else {
            xxs[j] <- 0
            yys[j] <- 0
          }
        }
        else {
          xxs[j] <- xx[j] * CGap/(seg + CGap) + xx[j] * 
            (df[i, j] - df[2, j])/(df[1, j] - df[2, 
                                                 j]) * seg/(seg + CGap)
          yys[j] <- yy[j] * CGap/(seg + CGap) + yy[j] * 
            (df[i, j] - df[2, j])/(df[1, j] - df[2, 
                                                 j]) * seg/(seg + CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i - 
                                                            2], border = pcols[i - 2], col = pfcols[i - 
                                                                                                      2])
      }
      else {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i - 
                                                            2], border = pcols[i - 2], density = pdensities[i - 
                                                                                                              2], angle = pangles[i - 2], col = pfcols[i - 
                                                                                                                                                         2])
      }
      points(xx * scale, yy * scale, pch = ptys[i - 2], 
             col = pcols[i - 2])
    }
  }
}

  radar.simulation <- function (species = "human", directory = "simulations/", pdf, text = NULL, a = 2) {
    
    # select species
    sim_data <- species
    code <- read.csv(paste0(directory, sim_data, "_simulation.code.txt"), header = T, sep=",")
    
    # create labels
    code$Lib_Plat <- apply(code, 1, function(x){
      paste(x["Library_Preps"], x["Platform"], sep = "_")
    })
    code$Lib_DC=apply(cbind(code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="_")
    code$Label <-apply(cbind(code[,c("Library_Preps", "Platform", "Data_Category")]), 1, paste, collapse="_")
    
    #set colors
    mycolors = c( "cDNA_PacBio_LO"="#d8527c", "cDNA_PacBio_LS"="#9a133d", 
                  "cDNA_ONT_LO"="#6996e3", "cDNA_ONT_LS"="#1a318b" )
    if(species == "mouse") { 
      names(mycolors) = sort(unique(code$Label))
    }
    
    # read data
    sim_file_prefix <- paste0("simulations/", sim_data, "_simulation_comparison.")
    sim_metrics <- read.csv(paste0(sim_file_prefix, "totsim_metrics.csv"), header = T, as.is = T) %>% t()
    sim0_metrics <- read.csv(paste0(sim_file_prefix, "totsim0_metrics.csv"), header = T, as.is = T) %>% t()
    sim1_metrics <- read.csv(paste0(sim_file_prefix, "totsim1_metrics.csv"), header = T, as.is = T) %>% t()
    sim5_metrics <- read.csv(paste0(sim_file_prefix, "totsim5_metrics.csv"), header = T, as.is = T) %>% t()
    simN_metrics <- read.csv(paste0(sim_file_prefix, "novelsim_metrics.csv"), header = T, as.is = T) %>% t()
    
    # calculate performance
    sim_SC <- paste0(directory, sim_data, "_simulation_comparison.summary_table_SC.csv")
    sim_SC <- read.csv(sim_SC, header = TRUE, as.is = TRUE) 
    sim_SC$total.novel <- sim_SC$total-sim_SC$FSM
    sim_metrics <- as.data.frame(sim_metrics)
    simN_metrics <- as.data.frame(simN_metrics)
    simN_metrics$FP <- sim_SC$total.novel
    sim_metrics$Precision_k <- sim_metrics[,"True Positive detections (TP)"] /sim_SC$FSM
    simN_metrics$Precision <- simN_metrics[,"Number of transcripts associated to TP (Reference Match)"]/(simN_metrics[,"Number of transcripts associated to TP (Reference Match)"]+ simN_metrics$FP)
    simN_metrics$nrPrecision <- simN_metrics[,"True Positive detections (TP)"] /(simN_metrics[,"True Positive detections (TP)"]+ simN_metrics$FP)
    Inv.Redundancy = round(1/sim_metrics[,"Redundancy"],2)
    metrics2 <- cbind(sim0_metrics[,c("Sensitivity")],
                      sim5_metrics[,c("Sensitivity")],
                      sim_metrics[,c("Precision_k")],
                      simN_metrics[,c("Sensitivity")],
                      simN_metrics[,c("Precision")],
                      Inv.Redundancy)
    
    colnames(metrics2)<- c("Sen_kn", "   Sen_kn \n >5TPM", "Pre_kn", "Sen_no", "Pre_no", "1/Red")
    
    # merge data
    metrics2 <- merge(metrics2, code, by.x=0, by.y="pipelineCode")
    metrics2$Tool <- gsub(pattern = "_", replacement = " ", x = metrics2$Tool)
    tool = sort(unique(metrics2$Tool))
    plat_lib <- unique(metrics2$Label)
    
    figure_name <- paste0(pdf, ".pdf" )
    pdf(figure_name ,width = 7, height = 8)
    #pdf(figure_name ,width = 10, height = 8) # size used for main paper panel
    
    par(mar=c(0,0.5,1.5,0) + 0.1)
    #layout(matrix(c(1:12,13,13,13,13),  byrow = TRUE, ncol = 4,nrow = 4)) # for main figure
    layout(matrix(c(1:12,13,13,13),  byrow = TRUE, ncol = 3,nrow = 5))
    
    sel <- metrics2$Tool == tool[11]
    subms <- metrics2[sel,c(2:7)]
    radarchart2(subms , axistype=1 ,vlabcol = c("darkgoldenrod4", "darkgoldenrod4", "darkgoldenrod4", "darkorchid4", "darkorchid4", "black"),
                #custom polygon
                pcol= mycolors[metrics2[sel,38]], plwd=2 , plty=1.3,
                #custom the grid
                cglcol="grey", cglty=1, axislabcol="grey", cglwd = 0.8, caxislabels=seq(0,0.25,1), 
                #custom labels
                vlabels = colnames(subms), vlcex=1.5, palcex = 1.5, cex.main=1.5)
    
    for ( i in 1:11) {
      sel <- metrics2$Tool == tool[i]
      subms <- metrics2[sel,c(2:7)]
      subms <- rbind(rep(1,ncol(subms)) , rep(0,ncol(subms)), subms)
      radarchart(subms, axistype=1 , title= tool[i],
                 #custom polygon
                 pcol= mycolors[metrics2[sel,16]], plwd=2 , plty=1.3,
                 #custom the grid
                 cglcol = "grey", cglty=1, axislabcol="white", cglwd = 0.8, caxislabels=seq(0,0.25,1), 
                 #custom labels
                 vlabels = rep("",6), palcex = 0, cex.main=1.5, paxislabels = 0
                 #vlabels = short.labls , palcex = 1, cex.main=1.7, paxislabels = 0, col = c(2,2,2,3,3,1)
      )
    }
    
    
    # Add a legend
    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    plot_colors <- c("blue","black", "green", "orange", "pink")
    l <- legend(x = "top",inset = 0, legend = names(mycolors), horiz = TRUE, bty = "n", pch=20 , col=mycolors , text.col = "black", cex=1.2, pt.cex=3)
    if (!is.null(text)) {
      text(x=l$text$x-a , y=l$text$y-l$rect$h, c(text,"","",""), adj=c(0,1))
    }
    
    dev.off()
  }

#### Gencode manual annotation analysis
# function to evaluate the SC of Genconde manually curated transcripts 
gencode.analysis <- function (gencode.self, uic, name) {
  gencode.data <- matrix(as.numeric(as.matrix(gencode.self[c(3:nrow(gencode.self)),c(2:7)])), nrow = nrow(gencode.self)-2, ncol = 6)
  colnames(gencode.data) <- colnames(gencode.self)[2:7]
  gencode.data <- as.data.frame(gencode.data)
  gencode.data$SC <- gencode.self[-c(1,2),"structural_category"]
  
  # expression distribution
  gencode.UIC <- sapply(gencode.self$Row.names[-c(1,2)], function (a) {paste(strsplit(a, "_")[[1]][-1], collapse = "_")}) # find UIC genecode
  UIC <- sapply(uic$LRGASP_id, function (a) {paste(strsplit(a, "_")[[1]][-1], collapse = "_")}) # find UIC
  genes.gencode <- unique(uic[(is.element(UIC, gencode.UIC)),"gene"])
  expression.gencode.genes <- uic[is.element(uic$gene, genes.gencode),]
  expression.gencode.genes <- tapply(expression.gencode.genes$FL_cpm, expression.gencode.genes$gene, sum)
  dens.data <- data.frame(CPM = c(uic[,"FL_cpm"],expression.gencode.genes),
                          Data = c(rep("All genes", nrow(uic)), rep("GENCODE manual", length(expression.gencode.genes))))
  
  A <- ggplot(dens.data , aes(x=CPM, color = Data, fill = Data)) +
    geom_density(alpha=0.6) + 
    theme_classic() +
    scale_x_log10() +
    ggtitle(paste0("Gene expression distribution ",name)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.title.y=element_blank(),
          legend.position = c(0.8, 0.4))
  
  #pie.chart
  pie.data <- as.data.frame(t(table(gencode.data$SC)))
  colnames(pie.data)[2] <- "Structural_Category"
  pie.data$Structural_Category <- as.vector(pie.data$Structural_Category)
  pie.data$Structural_Category[2] <- "GenicG"
  pie.data2 <- pie.data %>% 
    mutate(csum = rev(cumsum(rev(Freq))), 
           pos = Freq/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Freq/2, pos))
  
  B <- ggplot(pie.data2, aes(x="", y=Freq, fill=fct_inorder(Structural_Category))) +
    geom_col(width = 1, color = 1) +
    geom_text(aes(label = Freq),
              size=3,
              position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    scale_fill_manual(values=cat.palette[pie.data[,2]], name = "" ) +
    geom_label_repel(data = pie.data2,
                     aes(y = pos, label = Structural_Category),
                     segment.color = NA,
                     size = 3, nudge_x = 0.5, show.legend = FALSE)+
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 11, hjust=2), 
          legend.position = "none", # Removes the legend
          panel.background = element_rect(fill = "white"))
  
  #upsetplot
  preps = colnames(gencode.data)[1:6]
  C <- upset(
    gencode.data,
    preps,
    base_annotations=list(
      '# UIC'=intersection_size(
        counts=FALSE,
        mapping=aes(fill=SC)
      ) + scale_fill_manual(values=cat.palette, name="" )
      +  theme(legend.position = c(0.8,0.8), 
               legend.key.size = unit(0.3, 'cm'),
               panel.grid = element_blank())
    ),
    width_ratio=0.1 ,
    set_sizes=(
      upset_set_size(
        geom=geom_bar(
          aes(fill=SC),
          width=0.8
        ),
        position ='left'
      ) + scale_fill_manual(values=cat.palette, name="" ) 
      + theme(legend.position="none")
    ) +
      scale_y_continuous(position = "left") +
      theme(axis.text.y.left = element_text(angle = 270,vjust = -10))
  )
  return <- list(A = A, B = B, C= C)
  return
  
}

# function to read Genconde manually curated transcript IDs
read.gencode <- function (mypattern, directory, cols = 1, header = FALSE) {
  files <- dir(directory)
  gencode <- list()
  files.known <- grep(pattern = mypattern, files, value = TRUE)
  exclude <- grep("all", files.known) 
  if (length(exclude) > 0 ) {
    files.known <- files.known[-c(exclude)]
  }  
  for (i in 1: length(files.known)) {
    gencode [[i]] <- read.delim(paste0(directory, files.known[i]), as.is = TRUE, header = header) [,cols]
  }
  names(gencode) <- sapply(files.known , function (x) strsplit(x, split = "[.]")[[1]][2])
  gencode
}

# function to compute performance metrics of pipelines agains the GENCODE manual annotation
peformance.genecode <- function (gencode.pa, pa, code, directory, mypattern, evaluation, ID_UIC = ID_UIC, selection = NULL, gencode.self = NULL) {
  
  #read gencode files with the type of transcript (known or novel) each transcript is
  SQ3.data <- read.gencode(mypattern= mypattern, directory = directory, cols =c(1,6,7,8), header = TRUE)
  
  if(!is.null(selection)){
    isoform_selection <- ID_UIC[is.element(ID_UIC$LRGASP_id, selection),"isoform"]  # from UIC to transcript ID in the selection
    isoform_exclusion <- setdiff(ID_UIC$isoform, isoform_selection)
    SQ3.data <- lapply(SQ3.data, function (x) x[x$isoform %in% isoform_selection,])
    logical <- is.element(pa_GENCODE$associated_transcript,isoform_exclusion)
    b <- length(which(logical))
    gencode.pa$associated_transcript[logical] <- paste0("novel", c(1001:(1001+b-1)))
    gencode.pa$structural_category[logical] <- "NOVEL"
    }
  number.gencode <- data.frame (genes = sapply(SQ3.data, function (x) length(unique(x$associated_gene))),
                                transcripts =  sapply(SQ3.data , nrow), 
                                known = sapply(SQ3.data, function (x) length(which(x$structural_category == "full-splice_match"))),
                                novel = sapply(SQ3.data, function (x) length(which(x$structural_category != "full-splice_match"))))
  gencode.manual <- data.frame(transcript_id = unlist(sapply(SQ3.data, function (x) x$isoform)), 
                                library = rep(names(SQ3.data), sapply(SQ3.data , nrow)),
                                manual_vs_existing = unlist(sapply(SQ3.data, function (x) x$structural_category)))

  #obtain the UICs without the SC
  pa_SCs <- pa[,1:2] # this is the PA of the pipelines against the NORMAL annotation, with all genes
  colnames(pa_SCs)[2] <- "SC_normal_annotation"
  pa_SCs[,1] <- sapply(pa_SCs[,1], function (x) paste(strsplit(x, split = "_")[[1]][-1], collapse= "_"))
  
  # Additional variables to code
  code$Sample_code <- paste(code$Library_Preps, code$Platform, sep = "-")
  code$Label<- paste(code$Library_Preps, code$Platform, code$Data_Category, sep = "-")
  code$Lib_Plat <- paste(code$Library_Preps, code$Platform, sep = "-")

  # Merge in one table the SC of transcripts in the 47 pipelines when evaluated against normal annotation and against manual annotation
  pa.GENCODE.sample <- merge(gencode.pa, pa_SCs, by.x = "Row.names", by.y = "TAGS") # This contains the UIC of the pipelines associated to the 50 loci
  
  FSM_analysis <- matrix(0, ncol = 47, nrow = 10)
  colnames(FSM_analysis) = colnames(pa.GENCODE.sample)[2:48]
  for ( i in 1:47) {
    pipeline <- colnames(FSM_analysis)[i]
    library_pipeline <- code[code$pipelineCode == pipeline , "Sample_code"]
    # manual annotations on the prep
    
    h <- pa.GENCODE.sample[,c(pipeline ,"structural_category","SC_normal_annotation", "associated_transcript")]
    hh <- h[h[,1] == 1, ]
    hh_novel <- hh[grep("novel", hh$associated_transcript), "associated_transcript"] # identifies the novel transcripts of this pipeline
    annotated_genecode_library_specific <- c(SQ3.data[library_pipeline][[1]]$isoform, hh_novel) # list of transcripts specific library plus novel of the pipeline
    
    g <- pa.GENCODE.sample[is.element(pa.GENCODE.sample$associated_transcript, annotated_genecode_library_specific),] # select only the transcripts of that library preparation plus the novel
    d <- g[,c(pipeline ,"structural_category","SC_normal_annotation", "associated_transcript")]
    d <- d[d[,1] == 1, ]
    
    FSM_analysis[1,i] = length(which(d[,3] == "FSM")) # predicted transcript models identified as known against the normal gencode annotation. 
    FSM_analysis[2,i] = length(unique(d[d[,3] == "FSM","associated_transcript"])) # of these, how many unique transcript IDs in the manual annotation
    
    FSM_analysis[3,i] = length (which(d[,3] == "FSM" & d[,2] == "FSM")) # Known confirmed with GENCODE manual annotation
    FSM_analysis[4,i] = length(unique(d[d[,3] == "FSM" & d[,2] == "FSM","associated_transcript"]))  # number of known manually curated transcripts that were identifed by the pipelines
    
    FSM_analysis[5,i] = length(which(d[,3] != "FSM")) # predicted transcript models identified as novel against normal genome annotation
    FSM_analysis[6,i] = length(unique(d[d[,3] != "FSM","associated_transcript"])) # of these, how many unique transcript IDs in the manual annotation
    
    FSM_analysis[7,i] = length (which(d[,3] != "FSM" & d[,2] == "FSM")) # novel confirmed as TP against GENCODE manual annotation
    FSM_analysis[8,i] = length(unique(d[d[,3] != "FSM" & d[,2] == "FSM","associated_transcript"])) # of these, how many unique transcript IDs in the manual annotation
    
    FSM_analysis[9,i] = nrow(d) # total number of transcript models predicted by the pipeline
    FSM_analysis[10,i] = length(unique(d$associated_transcript)) # unique associated transcripts in the manual annotation
    
  }
  rownames(FSM_analysis) = c("Transcript_models_known", "Unique_detected_known", "TRUE_known", "Unique_TRUE_known",
                             "Transcript_models_novel", "Unique_detected_novel", "TRUE_novel", "Unique_TRUE_novel",
                             "Total_transcript_models", "Total_detections")
  
  FSM_analysis <- as.data.frame(t(FSM_analysis))
  
  # Total detections refers to the existing transcripts in the manual annotation that have been found at least once by the pipeline
  # Total transcripts refers to the number of transcripts reported by the pipeline
  # Difference between the two numbers is when the same reference transcript is reported twice, for example at ISM and as FSM
  
  FSM_analysis$Precision_total <- (FSM_analysis$TRUE_known + FSM_analysis$TRUE_novel) / FSM_analysis$Total_transcript_models
  FSM_analysis$nrPrecision_total <- (FSM_analysis$Unique_TRUE_known + FSM_analysis$Unique_TRUE_novel) / FSM_analysis$Total_detections
  FSM_analysis$Precision_known <- FSM_analysis$TRUE_known / FSM_analysis$Transcript_models_known
  FSM_analysis$nrPrecision_known <-  FSM_analysis$Unique_TRUE_known / FSM_analysis$Transcript_models_known
  FSM_analysis$Precision_novel <-  FSM_analysis$TRUE_novel / FSM_analysis$Transcript_models_novel
  FSM_analysis$nrPrecision_novel <-  FSM_analysis$Unique_TRUE_novel / FSM_analysis$Transcript_models_novel
  
  #Add information on ground truth to compute sensitivity
  number.gencode.2  <- merge(number.gencode, code, by.x = 0, by.y = "Lib_Plat" )
  evaluation1 <- merge(FSM_analysis, evaluation, by.x = 0, by.y = "pipelineCode",)
  evaluation2 <- merge(evaluation1, number.gencode.2, by.x = "Row.names", by.y = "pipelineCode",)
  
  evaluation2$Sensitivity_total <- (evaluation2$Unique_TRUE_known + evaluation2$Unique_TRUE_novel) / evaluation2$transcripts
  evaluation2$Sensitivity_known <- evaluation2$Unique_TRUE_known / evaluation2$know
  evaluation2$Sensitivity_novel <- evaluation2$Unique_TRUE_novel / evaluation2$novel
  evaluation2$F1_total <- (2 * evaluation2$Sensitivity_total * evaluation2$Precision_total) / (evaluation2$Sensitivity_total + evaluation2$Precision_total)
  evaluation2$F1_known <- (2 * evaluation2$Sensitivity_known * evaluation2$Precision_known) / (evaluation2$Sensitivity_known + evaluation2$Precision_known)
  evaluation2$F1_novel <- (2 * evaluation2$Sensitivity_novel * evaluation2$Precision_novel) / (evaluation2$Sensitivity_novel + evaluation2$Precision_novel)
  evaluation2$F1.Genes <- (2 * evaluation2$Sensitivity.Genes * evaluation2$Precision.Genes) / (evaluation2$Sensitivity.Genes + evaluation2$Precision.Genes)
  evaluation2$Redundancy_known <- evaluation2$Transcript_models_known / evaluation2$Unique_detected_known
  evaluation2$Redundancy_novel <- evaluation2$Transcript_models_novel / evaluation2$Unique_detected_novel
  evaluation2$Redundancy_total <- evaluation2$Total_transcript_models / evaluation2$Total_detections
  
  evaluation2
}

Performance_plot_left <- function (mydata, main,  axis.text.x = 11) { 
  ggplot(mydata, aes(x=value,y=Label, color=name, shape=name)) +
    geom_point(size=3) +
    scale_shape_manual(name="",
                       labels=c("F1 score", "Precision", "Sensitivity"),
                       values=c(4,19,19)) +
    scale_color_manual(name="",
                       labels=c("F1 score", "Precision","Sensitivity"),
                       values=c("#EE446F", "#F58A53","#15918A" )) +
    facet_grid(Alias ~., scales = "free_y", space="free") +
    pub_theme + 
    theme(
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size=7),
      axis.text.x = element_text(size=axis.text.x),
      title = element_text(size=8),
      legend.text = element_text(size=20),
      legend.key = element_rect(colour = NA, fill = NA)
    ) +
    ggtitle(main)+
    theme(plot.title = element_text(hjust = 0.5))
}
Performance_plot_right <- function (mydata, main,  axis.text.x = 11) { ggplot(mydata, aes(x=value,y=Label, color=name, shape=name))+
    geom_point(size=3)+
    scale_shape_manual(name="",
                       labels=c("F1 score", "Precision", "Sensitivity"),
                       values=c(4,19,19))+
    scale_color_manual(name="",
                       labels=c("F1 score", "Precision","Sensitivity"),
                       values=c("#EE446F", "#F58A53","#15918A" )) +
    facet_grid(Alias ~., scales = "free_y", space="free")+
    pub_theme+ 
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size= axis.text.x),
      title = element_text(size = 8),
      legend.position = "none",
      strip.text = element_text(size=8)
    ) +
    ggtitle(main)+
    theme(plot.title = element_text(hjust = 0.5))
}
Performance_plot_middel <- function (mydata, main,  axis.text.x = 11) { ggplot(mydata, aes(x=value,y=Label, color=name, shape=name))+
    geom_point(size=3)+
    scale_shape_manual(name="",
                       labels=c("F1 score", "Precision", "Sensitivity"),
                       values=c(4,19,19))+
    scale_color_manual(name="",
                       labels=c("F1 score", "Precision","Sensitivity"),
                       values=c("#EE446F", "#F58A53","#15918A" )) +
    facet_grid(Alias ~., scales = "free_y", space="free")+
    pub_theme+ 
    theme(
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      title =element_text(size = 8),
      axis.text.x = element_text(size = axis.text.x),
      legend.position = "none",
      #strip.text = element_text(size=8),
      legend.key = element_rect(colour = NA, fill = NA),
    ) +
    ggtitle(main)+
    theme(plot.title = element_text(hjust = 0.5))
}
Performance_plot_TP <- function (mydata, main,  axis.text.x = 11) { 
  ggplot(mydata, aes(x=value,y=Label, color=name, shape=name))+
    geom_point(size=3)+
    scale_shape_manual(name="",
                       labels=c("FALSE_known", "FALSE_novel", "TRUE_known", "TRUE_novel"),
                       values=c(19,4,19,4))+
    scale_color_manual(name="",
                       labels=c("FALSE_known", "FALSE_novel", "TRUE_known", "TRUE_novel"),
                       values=c("#EE446F", "#EE446F","#2d7Ac0", "#2d7Ac0" )) +
    facet_grid(Alias ~., scales = "free_y", space="free")+
    pub_theme+ 
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size= 7),
      axis.text.x = element_text(size= axis.text.x),
      title = element_text(size = 8),
      strip.text = element_text(size=8)
    ) +
    ggtitle(main)+
    theme(plot.title = element_text(hjust = 0.5))
}

# function to compute performance metrics of pipelines agains the GENCODE manual annotation imposing a selection. DOES NOT WORK YET. UNDER DEVELOPMENT
performance.genecode2 <- function (pa = pa.WTC11, pa_GENCODE = pa_GENCODE, code = code, gencode.self = gencode.self, selection = selection, directory) {
  
  # Calcualations on the manual curation data (number of known and number of novel transcripts)
  #############################################################################################
  # Selection of UIC in more than one sample
  gencode.self2 <- gencode.self[-c(1:2),c(1:9)]  
  gencode.self2[,1] <- sapply(gencode.self2[,1], function (x) paste(strsplit(x, split = "_")[[1]][-1], collapse= "_"))
  gencode.self3 <- matrix(as.numeric(as.matrix(gencode.self2[,c(2:7)])), nrow = nrow(gencode.self2), ncol = 6)
  colnames(gencode.self3) <- colnames(gencode.self2)[c(2:7)]
  rownames(gencode.self3) <- gencode.self2[,1] # matrix with UIC of the gencode manual transcripts indicating in which library they were found.
  
  gencode.k <- read.gencode (mypattern = "known_UIC", directory = directory)
  longK <- sapply(gencode.k , length)
  gencode.k  <- data.frame ( ID = unlist(gencode.k), sample = rep(names(longK), longK), status = "known")
  
  gencode.n <- read.gencode (mypattern = "novel_UIC", directory = directory)
  longN <- sapply(gencode.n , length)
  gencode.n  <- data.frame ( ID = unlist(gencode.n), sample = rep(names(longN), longN), status = "novel")
  
  gencode <- as.data.frame(rbind(gencode.k , gencode.n)) #annotations of the genome manual transcripts indicating if known or novel
  
  if(!is.null(selection)) {
    selection.uic <- selection[,1] # selected UIC
    selection2 <- selection[,c(2:7)] # in which libraries they are
    selection2 <- matrix(as.numeric(as.matrix(selection2)), nrow= nrow(selection2), ncol = ncol(selection2)) # the same as a matrix
    colnames(selection2) <- colnames(selection)[c(2:7)]
    rownames(selection2) <- selection.uic
    
    gencode.self2$selection <- rep("NO" ,nrow(gencode.self2))
    gencode.self2$selection[is.element(gencode.self2[,1], selection.uic)] <- "YES"  # add a YES if the UIC is in the selection
    table(gencode.self2$selection)
    gencode.self2 <- gencode.self2[gencode.self2$selection == "YES",] # select the selection
    gencode.self3 <- gencode.self3[gencode.self2$selection == "YES",] # select the selection
    
    selection_ID <- unique(pa_GENCODE[pa_GENCODE[,1] %in% selection.uic,"associated_transcript" ]) 
    pa_GENCODE[setdiff(pa_GENCODE[,1], selection_ID), "structural_category"] <-"NOVEL"
    gencode <- gencode[gencode$ID %in% selection_ID , ]
  }
  
  
  # Calcuating the KNOWN and NOVEL transcripts of the manual curation with respect to the normal annotation
  #calc.SC <- apply(gencode.self3, 2, function (x) tapply(x, gencode.self2$structural_category, sum))
  #Known_curated <- calc.SC[1,]
  #Novel_curated <- colSums(calc.SC[c(2:5),])
  
  #calc.SC.selec <- apply(selection2, 2, function (x) tapply(x, selection$structural_category, sum))
  #Known_curated_selection <- calc.SC.selec[1,]
  #Novel_curated_selection <-  colSums(calc.SC.selec[c(2:5),])
  
  # Calcuated the knonw and novel transcripts of the manual curation with respect to the normal annotation
  
  ###############################################################################################################
  
  # Merging evaluation 50 loci versus normal annotation and versus genecode manual annotation. Adding gencode.self data
  ########################################################################################################################
  
  # In pa_GENCODE evaluation doc (indicates PA of all transcripts in manual loci using GENCODE manual annotation as reference) add selection column
  pipelines <- colnames(pa_GENCODE) [2:48]
  
  #We should have here all the info to compute performance metrics
  Performance_analysis <- matrix(0, ncol = 47, nrow = 18)
  colnames(Performance_analysis) = pipelines
  for ( i in 1:length(pipelines)) {
    pipeline <- pipelines[i]
    library_pipeline <- code[code$pipelineCode == pipeline , "Sample_code"]
    #select the library in pa_GENCODE
    lib.GENCODE <- pa_GENCODE[,c("Row.names","structural_category", "associated_transcript", pipeline )]
    lib.GENCODE <- lib.GENCODE[lib.GENCODE[,pipeline] ==1,]
    lib.detections <- unique(lib.GENCODE [,1])# total number of transcripts of the pipeline in the 50 loc
    
    #select the library in pa
    lib.normal <- pa[,c("TAGS","structural_category", pipeline )]  # these are ALL transcripts of the pipeline
    lib.normal <- lib.normal[lib.normal[,pipeline] ==1,] ; nrow(lib.normal) # these are ALL transcripts of the pipeline
    
    sc.normal.loci <- lib.normal[lib.normal[,1] %in% lib.detections, "structural_category"]
    
    # filter for the transcripts that are part of the library preparation according to GENCODE
    lib.specific <- lib.GENCODE[lib.GENCODE$associated_transcript %in% gencode$ID[gencode$sample ==library_pipeline], ]
    nrow(lib.specific)
    
    #Now, merge the two
    merged.data <- merge(lib.normal , lib.specific, by.x = 1, by.y = 1)
    colnames(merged.data)[colnames(merged.data) =="structural_category.x"] <- "SC_against_normal"
    colnames(merged.data)[colnames(merged.data) =="structural_category.y"] <- "SC_against_manual"
    colnames(merged.data)[colnames(merged.data) =="associated_transcript"] <- "AT_manual"
    merged.data.u<- unique(merged.data[,c(1:4)])
    merged.data.u$"AT_manual" <- merged.data$AT_manual[match(merged.data.u$TAGS, merged.data$TAGS)]
    
    p <- merged.data.u
    
    Performance_analysis[1,i] = length(which(p[,"SC_against_normal"] == "FSM" & p[,"SC_against_manual"] == "FSM")) # True Known
    Performance_analysis[2,i] = length(unique(p[p[,"SC_against_normal"] == "FSM" & p[,"SC_against_manual"] == "FSM","AT_manual"])) # True known detections
    
    Performance_analysis[3,i] = length(which(p[,"SC_against_normal"] == "FSM" & p[,"SC_against_manual"] != "FSM")) # False known
    Performance_analysis[4,i] = length(unique(p[p[,"SC_against_normal"] == "FSM" & p[,"SC_against_manual"] != "FSM","AT_manual"])) # False known detections
    
    Performance_analysis[5,i] = length(which(p[,"SC_against_normal"] != "FSM" & p[,"SC_against_manual"] == "FSM")) # True novel
    Performance_analysis[6,i] = length(which(p[,"SC_against_normal"] != "FSM" & p[,"SC_against_manual"] != "FSM")) # False novel
    
    Performance_analysis[7,i]   <- length(gencode$ID[gencode$sample ==library_pipeline & gencode$status == "known"]) #ground.truth known
    Performance_analysis[8,i]   <- length(gencode$ID[gencode$sample ==library_pipeline & gencode$status == "novel"]) #ground.truth novel
    
    Performance_analysis[9,i]   <- length(which(sc.normal.loci == "FSM")) # detected as known
    Performance_analysis[10,i]   <- length(which(sc.normal.loci != "FSM"))# detected as novel
    
    Performance_analysis[11,i] <- Performance_analysis[1,i]/ Performance_analysis[7,i] #sensivity for known
    Performance_analysis[12,i] <- Performance_analysis[2,i]/ Performance_analysis[7,i] #non_redundant sensivity for known
    Performance_analysis[13,i] <- Performance_analysis[1,i]/ Performance_analysis[9,i] #precision for known
    Performance_analysis[14,i] <- Performance_analysis[1,i]/ Performance_analysis[9,i] #non_redundant precision for known
    
    Performance_analysis[15,i] <- Performance_analysis[5,i]/ Performance_analysis[8,i] #sensivity for novel
    Performance_analysis[16,i] <- Performance_analysis[6,i]/ Performance_analysis[8,i] #non_redundant sensivity for novel
    Performance_analysis[17,i] <- Performance_analysis[5,i]/ Performance_analysis[10,i] #precision for novel
    Performance_analysis[18,i] <- Performance_analysis[6,i]/ Performance_analysis[10,i] #non_redundant precision for novel
  }
  
  rownames(Performance_analysis) = c("True Known", "True known detections", "False known", "False known detections", 
                                     "True novel", "False novel", "Ground.truth known",  "Ground.truth novel", "Detected_known", 
                                     "Detected_novel", "Sensitivity_Known", "NR_Sensitivity_Known", "Precision_Known", "NR_Precision_Known", 
                                     "Sensitivity_Novel", "NR_Sensitivity_Novel", "Precision_Novel", "NR_Precision_Novel")
  Performance_analysis <- as.data.frame(t(Performance_analysis))
  Performance_analysis$F1_Known <- (2 * Performance_analysis$Sensitivity_Known * Performance_analysis$Precision_Known) / ( Performance_analysis$Sensitivity_Known + Performance_analysis$Precision_Known)
  Performance_analysis$F1_Novel <- (2 * Performance_analysis$Sensitivity_Novel * Performance_analysis$Precision_Novel) / ( Performance_analysis$Sensitivity_Novel + Performance_analysis$Precision_Novel)
  
  Performance_analysis <- merge(Performance_analysis, code , by.x = 0, by.y = "pipelineCode")
  Performance_analysis
  
}

performance.genecode3 <- function (pa = pa.WTC11, pa_GENCODE = pa_GENCODE, code = code, gencode.self = gencode.self, selection = selection, directory) {

  #read gencode files with the type of transcript (known or novel) each transcript is
gencode.FSM <- read.gencode(mypattern = "known_UIC", directory = directory); n1 = sapply(gencode.FSM , length); n1 # gencode manual annotation that are FSM of existing transcripts 
libraries.fsm <- names(gencode.FMS)
all.FSM <- unlist(gencode.FSM)
gencode.novel <- read.gencode(mypattern = "novel_UIC", directory = directory); n2 = sapply(gencode.novel , length); n2 # gencode manual annotation that are novel of existing transcripts 
all.novel <- unlist(gencode.novel)
libraries.novel <- names(gencode.novel)
gencode.manual <- data.frame(transcript_id = c(all.FMS, all.novel), 
                             library = c(rep(libraries.fsm, n1), rep(libraries.novel, n2)),
                             manual_vs_existing = c(rep("FSM", length(all.FSM)), rep("Novel", length(all.novel)))
)
head(gencode.manual)

if(!is.null(selection)) {
  not.selection <-  setdiff(gencode.manual$transcript_id, selection) # transcripts to exclude
  gencode.manual <- gencode.manual[gencode.manual$transcript_id %in% selection,]
  new.novel.positions <- which(pa_GENCODE$associated_transcript %in% not.selection)
  pa_GENCODE$structural_category[new.novel.positions] <- paste0("novel", 1:length(new.novel.positions))
}

ground.truth.numbers <- table(gencode.manual$library,gencode.manual$manual_vs_existing )

#obtain the UICs without the SC
pa_SCs <- pa[,1:2] # this is the PA of the pipelines against the NORMAL annotation
colnames(pa_SCs)[2] <- "SC_normal_annotation"
pa_SCs[,1] <- sapply(pa_SCs[,1], function (x) paste(strsplit(x, split = "_")[[1]][-1], collapse= "_"))

# Additional variables to code
code$Sample_code <- paste(code$Library_Preps, code$Platform, sep = "_")
code$Label<- paste(code$Library_Preps, code$Platform, code$Data_Category, sep = "_")
code$Lib_Plat <- paste(code$Library_Preps, code$Platform, sep = "_")

#merge the evaluation results against GENCODE and against the normal annotation
pa.GENCODE_2 <- merge(pa_GENCODE, pa_SCs, by.x = "Row.names", by.y = "TAGS")
}




# Local Variables:
# ess-indent-offset: 2
# End:
