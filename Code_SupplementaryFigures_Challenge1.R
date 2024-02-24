#####################################################
######### SUPPLEMENTARY FIGURES CHALLENGE 1 #########
#####################################################
## author: Francisco J. Pardo-Palacios, f.pardo.palacios@gmail.com
## author:Ana Conesa, ana.conesa@csic.es
## Last modified: March 16th 2023
#####################################################

outdir <- "output/extended"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

###############
source("Functions_Supplementary_Figures_Challenge1_v4.R")

### Data
code<- read.csv("Challenge1_Figures_Data/code.csv", as.is = TRUE,)
stats <- read.csv("Challenge1_Figures_Data/Stats_data.csv", as.is = TRUE)

Data.WTC11 <- read.csv("Challenge1_Figures_Data/WTC11_results/WTC11_comparison.summary_table_SC.csv", as.is = TRUE)
pa.WTC11 <- read.csv("Challenge1_Figures_Data/WTC11_results/WTC11_comparison.pa.csv", as.is = TRUE)
uic.WTC11 <- read.csv("Challenge1_Figures_Data/WTC11_results/WTC11_comparison.UIC_info.csv", as.is = TRUE)
bed.WTC11 <- read.delim("Challenge1_Figures_Data/WTC11_results/WTC11_comparison.bed", as.is = TRUE, header = FALSE)

Data.ES <- read.csv("Challenge1_Figures_Data/ES_results/ES_comparison.summary_table_SC.csv", as.is = TRUE)
pa.ES <- read.csv("Challenge1_Figures_Data/ES_results/ES_comparison.pa.csv", as.is = TRUE)
uic.ES <- read.csv("Challenge1_Figures_Data/ES_results/ES_comparison.UIC_info.csv", as.is = TRUE)

Data.H1mix <- read.csv("Challenge1_Figures_Data/H1_mix_results/H1_mix_comparison.summary_table_SC.csv", as.is = TRUE)
pa.H1mix <- read.csv("Challenge1_Figures_Data/H1_mix_results/H1_mix_comparison.pa.csv", as.is = TRUE)
uic.H1mix <- read.csv("Challenge1_Figures_Data/H1_mix_results/H1_mix_comparison.UIC_info.csv", as.is = TRUE)

gtf_human <- read.delim ("Challenge1_Figures_Data/human.gene_biotype.tsv", header = FALSE, as.is = TRUE)
gtf_human <- simplify.biotypes(gtf_human)
gtf_mouse <- read.delim ("Challenge1_Figures_Data/mouse.gene_biotype.tsv", header = FALSE, as.is = TRUE)
gtf_mouse <- simplify.biotypes(gtf_mouse)

gencode.self <- read.csv("Challenge1_Figures_Data/GENCODE_manualAnnot/presence_absence.agreement_within_GENCODE_manualAnnot.csv", as.is = TRUE)
gencode.self_mouse <- read.csv("Challenge1_Figures_Data/GENCODE_manualAnnot/GENCODE_mouse/GENCODE_manualAnnot.presence_absence.GENCODE_loci.csv", as.is = TRUE)

pa_GENCODE <- read.csv("Challenge1_Figures_Data/GENCODE_manualAnnot/presence_absence.GENCODE_loci_2.csv", sep=",", header = T) [,1:51] # Presence absence analysis of all transcripits of the 50 loci in pipelines evaluated against manual annotation. 
gencode_eval_results <- read.csv("Challenge1_Figures_Data/GENCODE_manualAnnot/new_GENCODE_manualAnnot_evaluation.csv", header = T) # evaluation result
pa_GENCODE_mouse <- read.csv("Challenge1_Figures_Data/GENCODE_manualAnnot/GENCODE_mouse/presence_absence.GENCODE_loci.csv", sep=",", header = T) [,1:51] # pa de pipelines evaluated against manual annotation
gencode_eval_results_mouse <- read.csv("Challenge1_Figures_Data/GENCODE_manualAnnot/GENCODE_mouse/GENCODE_manualAnnot_evaluation.csv", header = T) # evaluation result

gencode_human <- readxl::read_xlsx("Challenge1_Figures_Data/GENCODE_manualAnnot/LRGASP_human_read_support_percent.xlsx")
gencode_mouse <- readxl::read_xlsx("Challenge1_Figures_Data/GENCODE_manualAnnot/LRGASP_mouse_read_support_percent.xlsx")

ID_UIC <- read.delim("Challenge1_Figures_Data/GENCODE_manualAnnot/ID_UJC_manual_human/ID_UJC.txt", as.is = TRUE, header = TRUE )
gencode_human_SQ3 <- read.gencode (mypattern = "SQ3", directory = "Challenge1_Figures_Data/GENCODE_manualAnnot/classifications/human/") #SQ3 analysis of Gencode manual annotations human
gencode_mouse_SQ3 <- read.gencode (mypattern = "SQ3", directory = "Challenge1_Figures_Data/GENCODE_manualAnnot/classifications/mouse/") #SQ3 analysis of Gencode manual annotations human

usedreads<- read.delim("Challenge1_Figures_Data/NumberReadsInTranscripts.txt", as.is = TRUE, header = TRUE)

sim.data <- read.delim ("Challenge1_Figures_Data/Simulations/Coverage_simulated_data.txt",as.is = TRUE, header = FALSE )


## Extended Data Fig.1
##########################
usedreads$Sample <- paste(usedreads$BioSample, usedreads$Library_Preps, usedreads$Platform, sep = "_")
usedreads$Library_Platform <- paste(usedreads$Library_Preps, usedreads$Platform, sep = "_")
usedreads2 <- usedreads[usedreads$BioSample == "WTC11" | usedreads$BioSample == "H1-mix" | usedreads$BioSample == "Mouse ES",]
stats$Sample <- paste(stats$Bio_Sample2, stats$Method, stats$Tech, sep = "_")
usedreads2$Number_of_supplied_reads <- stats$Number_of_supplied_reads[match( usedreads2$Sample,stats$Sample)]
usedreads2$PercentageUsed <- (as.numeric(usedreads2$UsedReads) / usedreads2$Number_of_supplied_reads)*100
S1A <- fig.boxplots (usedreads2[usedreads2$BioSample == "WTC11",], sample = "WTC11",  var.y = "Lab", var.x = "PercentageUsed", 
                     xlabel = "PRU", jitter.color = "Library_Platform", mycolor = palette1) +
  scale_x_break(c(310, 1600), scales= 0.4) 

S1B <- fig.boxplots (usedreads2[usedreads2$BioSample == "H1-mix",], sample = "H1-mix",  var.y = "Lab", var.x = "PercentageUsed", 
                     xlabel = "PRU", jitter.color = "Library_Platform", mycolor = palette1) +
  scale_x_break(c(310, 900), scales= 0.4) 

S1C <- fig.boxplots (usedreads2[usedreads2$BioSample == "Mouse ES",], sample = "Mouse ES",  var.y = "Lab", var.x = "PercentageUsed", 
                     xlabel = "PRU", jitter.color = "Library_Platform", mycolor = palette1)+
  scale_x_break(c(310, 1375), scales= 0.4) 

trsPread <- read.delim("Challenge1_Figures_Data/reads2trans/TranscriptsPerReads.txt", as.is = TRUE, header = TRUE, row.names = 1)
trsPread2 <-round(trsPread / trsPread$Total,4)
trsPread2$Total <- rownames(trsPread2)
trsPread3 <- melt(trsPread2)
colnames(trsPread3) <-c("Tool", "Transcripts_per_read", "Fraction")
S1D <- ggplot(trsPread3, aes(fill=Transcripts_per_read, y=Fraction, x=Tool)) +
  geom_bar(position="stack", stat="identity") +
  ggtitle("Distribution number of transcripts per long-read") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_conesa(palette = "complete") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
  theme(legend.title = element_text(size=9), legend.text = element_text(size=9))

# Figures were copied to the clipboard, size adjusted, saved and used to compose Figure S1 in PowerPoint


## Extended Data Fig. 2 and 3  Detection plots for H1-mix and Mouse ES
#########################################################################

FirstPanelChl1 (data_sample = "H1_mix", outdir = outdir, ylims = c(160000, 250000), xlims = c(12000, 28000))
FirstPanelChl1 (data_sample = "ES", outdir = outdir, ylims = c(100000, 300000), xlims = c(17000, 20000))

# Plots are generated separately.Supplementary figures 2 and 3 are composed in PowerPoint.

## Extended Data Fig.4
#########################
data.WTC11 <-  number.genes.transcripts(data = Data.WTC11, code = code, pa = pa.WTC11, uic = uic.WTC11, stats = stats[stats$Bio_Sample == "WTC11",])
W <- fig.number.genes.transcripts(data.WTC11$data.genes, sample = "WTC11", var = "Number_of_reads")

data.ES <-  number.genes.transcripts(data = Data.ES, code = code, pa = pa.ES, uic = uic.ES, stats = stats[stats$Bio_Sample == "ES_mouse",])
E <- fig.number.genes.transcripts(data.ES$data.genes, sample = "Mouse ES", var = "Number_of_reads")

data.H1 <-  number.genes.transcripts(data = Data.H1mix, code = code, pa = pa.H1mix, uic = uic.H1mix, stats = stats[stats$Bio_Sample == "H1_mix",])
H <- fig.number.genes.transcripts(data.H1$data.genes, sample = "H1-mix", var = "Number_of_reads")

suppl = "4"
figureS4 <- ggarrange(W$A,H$A,E$A, W$B,H$B,E$B,
                    labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                    ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Relationship between sequencing depth and number of detected features. a-c) Transcripts, \n     d-f) Genes.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_seq-depth-by-detected-features",".pdf"))
annotate_figure(figureS4,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.5
#########################
W1 <- fig.number.genes.transcripts(data.WTC11$data.genes, sample = "WTC11", var = "Median_Read_length")
E1 <- fig.number.genes.transcripts(data.ES$data.genes, sample = "Mouse ES", var = "Median_Read_length")
H1 <- fig.number.genes.transcripts(data.H1$data.genes, sample = "H1-mix", var = "Median_Read_length")

suppl = "5"
figureS5 <- ggarrange(W1$A,H1$A,E1$A, W1$B,H1$B,E1$B,
                    labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                    ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Relationship between read length and number of detected features. a-c) Transcripts, \n     d-f) Genes.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_read-length-by-detected-features",".pdf"))
annotate_figure(figureS5,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.6
#########################
W2 <- fig.number.genes.transcripts(data.WTC11$data.genes, sample = "WTC11", var = "Median_Identity")
E2 <- fig.number.genes.transcripts(data.ES$data.genes, sample = "Mouse ES", var = "Median_Identity")
H2 <- fig.number.genes.transcripts(data.H1$data.genes, sample = "H1-mix", var = "Median_Identity")

suppl = "6"
figureS6 <- ggarrange(W2$A,H2$A,E2$A, W2$B,H2$B,E2$B,
                    labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                    ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Relationship between read quality and number of detected features. a-c) Transcripts, \n     d-f) Genes.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_read-qual-by-detected-features",".pdf"))
annotate_figure(figureS6,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.7
#########################
S7W1 <- fig.boxplots (data.WTC11$mad.data.frame, sample = "WTC11", my.xlim = c(0,32), var.x = "mad_transcripts", xlabel = "mad_transcripts (k)")
S7W2 <- fig.boxplots (data.WTC11$mad.data.frame, sample = "WTC11", my.xlim = c(0,9), var.x = "mad_genes", xlabel = "mad_genes (k)")
S7H1 <- fig.boxplots (data.H1$mad.data.frame, sample = "H1-mix", my.xlim = c(0,32), var.x = "mad_transcripts", xlabel = "mad_transcripts (k)")
S7H2 <- fig.boxplots (data.H1$mad.data.frame, sample = "H1-mix", my.xlim = c(0,9), var.x = "mad_genes" ,xlabel = "mad_genes (k)")
S7E1 <- fig.boxplots (data.ES$mad.data.frame, sample = "Mouse ES", my.xlim = c(0,32), var.x = "mad_transcripts", xlabel = "mad_transcripts (k)")
S7E2 <- fig.boxplots (data.ES$mad.data.frame, sample = "Mouse ES", my.xlim = c(0,9), var.x = "mad_genes", xlabel = "mad_genes (k)")

suppl = "7"
figureS7 <- ggarrange(S7W1,S7H1,S7E1, S7W2,S7H2,S7E2,
                     labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                     ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Median Absolute Deviance of detected features by experimental factor. a-c) Transcripts, \n     d-f) Genes.")
pdf(paste0(outdir, "/Extended_Fig._",suppl, "_deviance-dectected-features-by-expr-factors",".pdf"))
annotate_figure(figureS7,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.8
#########################
S8W1 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11", var.x = "Number_of_genes", my.xlim = c(250, 50000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Sample", mycolor = palette1)
S8W2 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11", var.x = "Number_of_transcripts", my.xlim = c(250, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1)
S8H1 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix", var.x = "Number_of_genes", my.xlim = c(250, 50000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Sample", mycolor = palette1)
S8H2 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix", var.x = "Number_of_transcripts" , my.xlim = c(250, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1)
S8E1 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES", var.x = "Number_of_genes", my.xlim = c(250, 50000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Sample", mycolor = palette1)
S8E2 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES", var.x = "Number_of_transcripts", my.xlim = c(250, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1)

suppl = "8"
figureS8 <- ggarrange(S8W2,S8H2,S8E2, S8W1,S8H1,S8E1,
                      labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                      ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of detected transcripts and genes per analysis tool.  a-c) Transcripts, d-f) Genes.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_detected-trans-gene-by-tool", ".pdf"))
annotate_figure(figureS8,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.9
#########################
S9W1 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11", var.x = "Number_of_genes", my.xlim = c(500, 32000), rescale = TRUE, var.y = "Platform", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Library_Preps", mycolor = palette2)
S9W2 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11", var.x = "Number_of_genes", my.xlim = c(500, 32000), rescale = TRUE, var.y = "Library_Preps", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Platform", mycolor = palette2)
S9H1 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix", var.x = "Number_of_genes", my.xlim = c(500, 32000), rescale = TRUE, var.y = "Platform", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Library_Preps", mycolor = palette2)
S9H2 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix", var.x = "Number_of_genes" , my.xlim = c(500, 32000), rescale = TRUE, var.y = "Library_Preps", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Platform", mycolor = palette2)
S9E1 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES", var.x = "Number_of_genes", my.xlim = c(500, 32000), rescale = TRUE, var.y = "Platform", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Library_Preps", mycolor = palette2)
S9E2 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES", var.x = "Number_of_genes", my.xlim = c(500, 32000), rescale = TRUE, var.y = "Library_Preps", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Platform", mycolor = palette2)

legend1 <- get_legend(S9W1)
legend2 <- get_legend(S9E2)

suppl = "9"

figureS91 <- ggarrange(S9W1,S9H1,S9E1,
                       labels = c( "a)", "b)", "c)"),
                       ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom") 
figureS92 <- ggarrange( S9W2,S9H2,S9E2,
                       labels = c(  "d)", "e)", "f)"),
                       ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom") 

figureS9 <- ggarrange ( figureS91, figureS92, ncol = 1, nrow = 2) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 

mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of detected genes per Platform and Library Preparation. a-c) Platform, \n     d-f) Library Preparation.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_detected-genes-by-plat-prep",".pdf"))
annotate_figure(figureS9,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.10
##########################
S10W1 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11", var.x = "Number_of_transcripts", my.xlim = c(500, 200000), rescale = TRUE, var.y = "Platform", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2)
S10W2 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11", var.x = "Number_of_transcripts", my.xlim = c(500, 200000), rescale = TRUE, var.y = "Library_Preps", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)
S10H1 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix", var.x = "Number_of_transcripts", my.xlim = c(500, 200000), rescale = TRUE, var.y = "Platform", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2)
S10H2 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix", var.x = "Number_of_transcripts" , my.xlim = c(500, 200000), rescale = TRUE, var.y = "Library_Preps", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)
S10E1 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES", var.x = "Number_of_transcripts", my.xlim = c(500, 200000), rescale = TRUE, var.y = "Platform", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2)
S10E2 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES", var.x = "Number_of_transcripts", my.xlim = c(500, 200000), rescale = TRUE, var.y = "Library_Preps", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)

suppl = "10"
figureS101 <- ggarrange(S10W1,S10H1,S10E1,
                       labels = c( "a)", "b)", "c)"),
                       ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom") 
figureS102 <- ggarrange( S10W2,S10H2,S10E2,
                        labels = c(  "d)", "e)", "f)"),
                        ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom") 

figureS10 <- ggarrange ( figureS101, figureS102, ncol = 1, nrow = 2) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 

mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of detected transcripts per Platform and Library Preparation. a-c) Platform, \n     d-f) Library Preparation.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_detected-genes-by-plat-prep",".pdf"))
annotate_figure(figureS10,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.11
##########################
S11W1 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Library_Preps"] == "cDNA",], sample = "WTC11_cDNA", var.x = "Number_of_transcripts", my.xlim = c(3000, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 9)
S11H1 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Library_Preps"] == "cDNA",], sample = "H1-mix_cDNA", var.x = "Number_of_transcripts", my.xlim = c(3000, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 9)
S11E1 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Library_Preps"] == "cDNA",], sample = "Mouse ES_cDNA", var.x = "Number_of_transcripts", my.xlim = c(3000, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 9)

S11W2 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Library_Preps"] == "CapTrap",], sample = "WTC11_CapTrap", var.x = "Number_of_transcripts", my.xlim = c(3000, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 9)
S11H2 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Library_Preps"] == "CapTrap",], sample = "H1-mix_CapTrap", var.x = "Number_of_transcripts", my.xlim = c(3000, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 9)
S11E2 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Library_Preps"] == "CapTrap",], sample = "Mouse ES_CapTrap", var.x = "Number_of_transcripts", my.xlim = c(3000, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 9)

suppl = "11"

figureS11 <- ggarrange(S11W1,S11H1,S11E1, S11W2,S11H2,S11E2, 
                       labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                      ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of detected transcripts in cDNA and CapTrap libraries. a-c) cDNA, d-f) CapTrap.") 
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_detected-trans-by-cDNA-CapTrap",".pdf"))
annotate_figure(figureS11,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.12
##########################
S11W3 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Platform"] == "PacBio",], sample = "WTC11_PacBio", var.x = "Number_of_transcripts", my.xlim = c(3000, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 9)
S11H3 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Platform"] == "PacBio",], sample = "H1-mix_PacBio", var.x = "Number_of_transcripts", my.xlim = c(3000, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 9)
S11E3 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Platform"] == "PacBio",], sample = "Mouse ES_PacBio", var.x = "Number_of_transcripts", my.xlim = c(3000, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 9)

S11W4 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Platform"] == "ONT",], sample = "WTC11_ONT", var.x = "Number_of_transcripts", my.xlim = c(3000, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 9)
S11H4 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Platform"] == "ONT",], sample = "H1-mix_ONT", var.x = "Number_of_transcripts", my.xlim = c(3000, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 9)
S11E4 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Platform"] == "ONT",], sample = "Mouse ES_ONT", var.x = "Number_of_transcripts", my.xlim = c(3000, 200000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 9)

suppl = "12"

figureS121 <- ggarrange(S11W3 ,S11H3,S11E3,
                        labels = c( "a)", "b)", "c)"),
                        ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom") 
figureS122 <- ggarrange( S11W4,S11H4,S11E4,
                         labels = c(  "d)", "e)", "f)"),
                         ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom") 

figureS12 <- ggarrange ( figureS121, figureS122, ncol = 1, nrow = 2) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 

mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of detected transcripts in PacBio and Nanopore platforms. a-c) PacBio, d-f) Nanopore.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_detected-trans-by-PacBio-ONT",".pdf"))
annotate_figure(figureS12,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.13
##########################
S13W1 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Library_Preps"] == "cDNA",], sample = "WTC11_cDNA", var.x = "Number_of_genes", my.xlim = c(3000, 32000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Platform", mycolor = palette2, title.size = 9)
S13H1 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Library_Preps"] == "cDNA",], sample = "H1-mix_cDNA", var.x = "Number_of_genes", my.xlim = c(3000, 32000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Platform", mycolor = palette2, title.size = 9)
S13E1 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Library_Preps"] == "cDNA",], sample = "Mouse ES_cDNA", var.x = "Number_of_genes", my.xlim = c(3000, 32000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Platform", mycolor = palette2, title.size = 9)

S13W2 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Library_Preps"] == "CapTrap",], sample = "WTC11_CapTrap", var.x = "Number_of_genes", my.xlim = c(3000, 32000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Platform", mycolor = palette2, title.size = 9)
S13H2 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Library_Preps"] == "CapTrap",], sample = "H1-mix_CapTrap", var.x = "Number_of_genes", my.xlim = c(3000, 32000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Platform", mycolor = palette2, title.size = 9)
S13E2 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Library_Preps"] == "CapTrap",], sample = "Mouse ES_CapTrap", var.x = "Number_of_genes", my.xlim = c(3000, 32000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Platform", mycolor = palette2, title.size = 9)

suppl = "13"
figureS13 <- ggarrange(S13W1,S13H1,S13E1, S13W2,S13H2,S13E2, 
                       labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                       ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of detected genes in cDNA and CapTrap libraries. a-c) cDNA, d-f) CapTrap.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_detected-genes-by-PacBio-ONT",".pdf"))
annotate_figure(figureS13,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.14
##########################
S13W3 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Platform"] == "PacBio",], sample = "WTC11_PacBio", var.x = "Number_of_genes", my.xlim = c(500, 32000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Library_Preps", mycolor = palette2, title.size = 9)
S13H3 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Platform"] == "PacBio",], sample = "H1-mix_PacBio", var.x = "Number_of_genes", my.xlim = c(500, 32000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Library_Preps", mycolor = palette2, title.size = 9)
S13E3 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Platform"] == "PacBio",], sample = "Mouse ES_PacBio", var.x = "Number_of_genes", my.xlim = c(500, 32000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Library_Preps", mycolor = palette2, title.size = 9)

S13W4 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Platform"] == "ONT",], sample = "WTC11_ONT", var.x = "Number_of_genes", my.xlim = c(500, 32000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Library_Preps", mycolor = palette2, title.size = 9)
S13H4 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Platform"] == "ONT",], sample = "H1-mix_ONT", var.x = "Number_of_genes", my.xlim = c(500, 32000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Library_Preps", mycolor = palette2, title.size = 9)
S13E4 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Platform"] == "ONT",], sample = "Mouse ES_ONT", var.x = "Number_of_genes", my.xlim = c(500, 32000), rescale = TRUE, var.y = "Tool", jitter.size= 1.5, xlabel = "# genes", jitter.color = "Library_Preps", mycolor = palette2, title.size = 9)

suppl = "14"

figureS141 <- ggarrange(S13W3 ,S13H3,S13E3,
                        labels = c( "a)", "b)", "c)"),
                        ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom") 
figureS142 <- ggarrange( S13W4,S13H4,S13E4,
                         labels = c(  "d)", "e)", "f)"),
                         ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom") 

figureS14 <- ggarrange ( figureS141, figureS142, ncol = 1, nrow = 2) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 

mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of detected genes in PacBio and Nanopore platforms. a-c) PacBio, d-f) Nanopore.")

pdf(paste0(outdir, "/Extended_Fig._",suppl,"_detected-genes-by-PacBio-ONT",".pdf"))
annotate_figure(figureS14,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.15
##########################
S15W1 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11_FSM", var.x = "FSM", var.y = "Library_Preps", my.xlim = c(500, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)
S15H1 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix_FSM", var.x = "FSM", var.y = "Library_Preps", my.xlim = c(500, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)
S15E1 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES_FSM", var.x = "FSM", var.y = "Library_Preps", my.xlim = c(500, 100000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)

S15W2 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11_ISM", var.x = "ISM", var.y = "Library_Preps", my.xlim = c(500, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)
S15H2 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix_ISM", var.x = "ISM", var.y = "Library_Preps", my.xlim = c(500, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)
S15E2 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES_ISM", var.x = "ISM", var.y = "Library_Preps", my.xlim = c(500, 100000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)

suppl = "15"
figureS15 <- ggarrange(S15W1,S15H1,S15E1, S15W2,S15H2,S15E2, 
                       labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                       ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of FSM and ISM by sequencing platform and library preparation. a-c) FSM, d-f) ISM.") 
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_FSM-ISM-by-plat-lib",".pdf"))
annotate_figure(figureS15,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.16
##########################
S16W1 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11_NIC", var.x = "NIC", var.y = "Library_Preps", my.xlim = c(500, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)
S16H1 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix_NIC", var.x = "NIC", var.y = "Library_Preps", my.xlim = c(500, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)
S16E1 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES_NIC", var.x = "NIC", var.y = "Library_Preps", my.xlim = c(500, 100000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)

S16W2 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11_NNC", var.x = "NNC", var.y = "Library_Preps", my.xlim = c(500, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)
S16H2 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix_NNC", var.x = "NNC", var.y = "Library_Preps", my.xlim = c(500, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)
S16E2 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES_NNC", var.x = "NNC", var.y = "Library_Preps", my.xlim = c(500, 100000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2)

suppl = "16"
figureS16 <- ggarrange(S16W1,S16H1,S16E1, S16W2,S16H2,S16E2, 
                       labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                       ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of NIC and NNC by sequencing platform and library preparation. a-c) NIC, d-f) NNC.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_NIC-NCC-by-plat-prep",".pdf"))
annotate_figure(figureS16,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.17
##########################
S17W1 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Library_Preps"] == "cDNA",], sample = "WTC11_cDNA", var.x = "FSM", var.y = "Tool", my.xlim = c(500, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 7, xlabel.size = 8)
S17H1 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Library_Preps"] == "cDNA",], sample = "H1-mix_cDNA", var.x = "FSM", var.y = "Tool", my.xlim = c(500, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 7, xlabel.size = 8)
S17E1 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Library_Preps"] == "cDNA",], sample = "Mouse ES_cDNA", var.x = "FSM", var.y = "Tool", my.xlim = c(500, 100000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 7, xlabel.size = 8)

S17W2 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Library_Preps"] == "CapTrap",], sample = "WTC11_CapTrap", var.x = "FSM", var.y = "Tool", my.xlim = c(500, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 7, xlabel.size = 8)
S17H2 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Library_Preps"] == "CapTrap",], sample = "H1-mix_CapTrap", var.x = "FSM", var.y = "Tool", my.xlim = c(500, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 7, xlabel.size = 8)
S17E2 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Library_Preps"] == "CapTrap",], sample = "Mouse ES_CapTrap", var.x = "FSM", var.y = "Tool", my.xlim = c(500, 100000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 7, xlabel.size = 8)

suppl = "17"
figureS17 <- ggarrange(S17W1,S17H1,S17E1, S17W2,S17H2,S17E2, 
                       labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                       ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of FSM transcripts by library preparation and analysis tool. a-c) cDNA. d-f) CapTrap.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_FSM-by-lib-tool",".pdf"))
annotate_figure(figureS17,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.18
##########################
S18W1 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Platform"] == "PacBio",], sample = "WTC11_PacBio", var.x = "FSM", var.y = "Tool", my.xlim = c(900, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 7, xlabel.size = 8)
S18H1 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Platform"] == "PacBio",], sample = "H1-mix_PacBio", var.x = "FSM", var.y = "Tool", my.xlim = c(900, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 7, xlabel.size = 8)
S18E1 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Platform"] == "PacBio",], sample = "Mouse ES_PacBio", var.x = "FSM", var.y = "Tool", my.xlim = c(900, 100000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 7, xlabel.size = 8)

S18W2 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Platform"] == "ONT",], sample = "WTC11_ONT", var.x = "FSM", var.y = "Tool", my.xlim = c(900, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 7, xlabel.size = 8)
S18H2 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Platform"] == "ONT",], sample = "H1-mix_ONT", var.x = "FSM", var.y = "Tool", my.xlim = c(900, 100000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 7, xlabel.size = 8)
S18E2 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Platform"] == "ONT",], sample = "Mouse ES_ONT", var.x = "FSM", var.y = "Tool", my.xlim = c(900, 100000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 7, xlabel.size = 8)

suppl = "18"
figureS18 <- ggarrange(S18W1,S18H1,S18E1, S18W2,S18H2,S18E2, 
                       labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                       ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of FSM transcripts by sequencing platform and analysis tool. a-c) PacBio, \n     d-f) Nanopore.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_FSM-by-plat-tool",".pdf"))
annotate_figure(figureS18,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.19
##########################
S19W1 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Library_Preps"] == "cDNA",], sample = "WTC11_cDNA", var.x = "ISM", var.y = "Tool", my.xlim = c(500, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 7, xlabel.size = 8)
S19H1 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Library_Preps"] == "cDNA",], sample = "H1-mix_cDNA", var.x = "ISM", var.y = "Tool", my.xlim = c(500, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 7, xlabel.size = 8)
S19E1 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Library_Preps"] == "cDNA",], sample = "Mouse ES_cDNA", var.x = "ISM", var.y = "Tool", my.xlim = c(500, 50000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 7, xlabel.size = 8)

S19W2 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Library_Preps"] == "CapTrap",], sample = "WTC11_CapTrap", var.x = "ISM", var.y = "Tool", my.xlim = c(500, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 7, xlabel.size = 8)
S19H2 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Library_Preps"] == "CapTrap",], sample = "H1-mix_CapTrap", var.x = "ISM", var.y = "Tool", my.xlim = c(500, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 7, xlabel.size = 8)
S19E2 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Library_Preps"] == "CapTrap",], sample = "Mouse ES_CapTrap", var.x = "ISM", var.y = "Tool", my.xlim = c(500, 50000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Platform", mycolor = palette2, title.size = 7, xlabel.size = 8)

suppl = "19"
figureS19 <- ggarrange(S19W1,S19H1,S19E1,S19W2,S19H2,S19E2, 
                       labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                       ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of ISM transcripts by library preparation and analysis tool. a-c) cDNA. d-f) CapTrap.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_ISM-by-prep-tool",".pdf"))
annotate_figure(figureS19,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.20
##########################
S20W1 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Platform"] == "PacBio",], sample = "WTC11_PacBio", var.x = "ISM", var.y = "Tool", my.xlim = c(900, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 7, xlabel.size = 8)
S20H1 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Platform"] == "PacBio",], sample = "H1-mix_PacBio", var.x = "ISM", var.y = "Tool", my.xlim = c(900, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 7, xlabel.size = 8)
S20E1 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Platform"] == "PacBio",], sample = "Mouse ES_PacBio", var.x = "ISM", var.y = "Tool", my.xlim = c(900, 50000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 7, xlabel.size = 8)

S20W2 <- fig.boxplots (data.WTC11$data.genes[data.WTC11$data.genes[,"Platform"] == "ONT",], sample = "WTC11_ONT", var.x = "ISM", var.y = "Tool", my.xlim = c(900, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 7, xlabel.size = 8)
S20H2 <- fig.boxplots (data.H1$data.genes[data.H1$data.genes[,"Platform"] == "ONT",], sample = "H1-mix_ONT", var.x = "ISM", var.y = "Tool", my.xlim = c(900, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 7, xlabel.size = 8)
S20E2 <- fig.boxplots (data.ES$data.genes[data.ES$data.genes[,"Platform"] == "ONT",], sample = "Mouse ES_ONT", var.x = "ISM", var.y = "Tool", my.xlim = c(900, 50000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Library_Preps", mycolor = palette2, title.size = 7, xlabel.size = 8)

suppl = "20"
figureS20 <- ggarrange(S19W1,S19H1,S19E1,S19W2,S19H2,S19E2, 
                       labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                       ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of ISM transcripts by sequencing platform and analysis tool. a-c) Intergenic. \n     d-f) GenicGenomic.") 
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_ISM-by-plat-tool",".pdf"))
annotate_figure(figureS20,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.21
##########################
S22W1 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11_Intergenic", var.x = "Intergenic", var.y = "Tool", my.xlim = c(10, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1, title.size = 7, xlabel.size = 8)
S22H1 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix_Intergenic", var.x = "Intergenic", var.y = "Tool", my.xlim = c(10, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1, title.size = 7, xlabel.size = 8)
S22E1 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES_Intergenic", var.x = "Intergenic", var.y = "Tool", my.xlim = c(10, 50000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1, title.size = 7, xlabel.size = 8)
  
S22W2 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11_GenicGenomic", var.x = "GenicGenomic", var.y = "Tool", my.xlim = c(10, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1, title.size = 7, xlabel.size = 8)
S22H2 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix_GenicGenomic", var.x = "GenicGenomic", var.y = "Tool", my.xlim = c(10, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1, title.size = 7, xlabel.size = 8)
S22E2 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES_GenicGenomic", var.x = "GenicGenomic", var.y = "Tool", my.xlim = c(10, 50000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1, title.size = 7, xlabel.size = 8)
  
suppl = "21"
figureS21 <- ggarrange(S22W1,S22H1,S22E1, S22W2,S22H2,S22E2, 
                         labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                         ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
   theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of Intergenic and GenicGenomic by sequencing platform and library preparation. \n     a-c) Intergenic, d-f) GenicGenomic.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_genic-by-plat-prep",".pdf"))
annotate_figure(figureS21,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()
  
## Extended Data Fig.22
##########################
S23W1 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11_Fusion", var.x = "Fusion", var.y = "Tool", my.xlim = c(10, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1, title.size = 7, xlabel.size = 8)
S23H1 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix_Fusion", var.x = "Fusion", var.y = "Tool", my.xlim = c(10, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1, title.size = 7, xlabel.size = 8)
S23E1 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES_Fusion", var.x = "Fusion", var.y = "Tool", my.xlim = c(10, 50000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1, title.size = 7, xlabel.size = 8)
  
S23W2 <- fig.boxplots (data.WTC11$data.genes, sample = "WTC11_Antisense", var.x = "Antisense", var.y = "Tool", my.xlim = c(10, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1, title.size = 7, xlabel.size = 8)
S23H2 <- fig.boxplots (data.H1$data.genes, sample = "H1-mix_Antisense", var.x = "Antisense", var.y = "Tool", my.xlim = c(10, 50000), rescale = TRUE,  jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1, title.size = 7, xlabel.size = 8)
S23E2 <- fig.boxplots (data.ES$data.genes, sample = "Mouse ES_Antisense", var.x = "Antisense", var.y = "Tool", my.xlim = c(10, 50000), rescale = TRUE, jitter.size= 1.5, xlabel = "# transcripts", jitter.color = "Sample", mycolor = palette1, title.size = 7, xlabel.size = 8)
  
suppl = "22"
figureS22 <- ggarrange(S23W1,S23H1,S23E1, S23W2,S23H2,S23E2, 
                         labels = c( "a)", "b)", "c)", "d)", "e)", "f)"),
                         ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom") +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of Fusion and Antisense by sequencing platform and library preparation. \n     a-c) Fusion. d-f) Antisense.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_fusion-anti-by-plat-prep",".pdf"))
annotate_figure(figureS22,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.23
##########################
S23A <- LRC(directory  = "Challenge1_Figures_Data/coverage_files", sample = "WTC11", main = "WTC11")
S23B <- LRC(directory  = "Challenge1_Figures_Data/coverage_files", sample = "H1mix", main = "H1-mix")
S23C <- LRC(directory  = "Challenge1_Figures_Data/coverage_files", sample = "ESmouse", main = "Mouse ES")
  
suppl = "23"
figureS23 <- ggarrange(S23A,S23B,S23C,
                         labels = c( "a)", "b)", "c)"),
                         ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.5,0.2,0.5, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Percentage of transcript models (TM) with different ranges of sequence coverage by long reads.\n     a) WTC11. c) H1-mix. c) Mouse ES. Ba: Bambu, FM: FLAMES, FL: FLAIR, IQ: IsoQuant, IT: IsoTools, IB: Iso_IB, Ly: LyRic, \n     Ma: Mandalorion, TL: TALON-LAPA, Sp: Spectra, ST: StringTie2.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_TMs-coverage-by-tool",".pdf"))
annotate_figure(figureS23,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()
  
## Extended Data Fig.24
##########################
S24A  <- biotype.plot(pa = pa.WTC11 , info = uic.WTC11 , code = code, gtf = gtf_human, name = "WTC11")
S24B  <- biotype.plot(pa = pa.H1mix , info = uic.H1mix, code = code, gtf = gtf_human, name = "H1-mix")
S24C  <- biotype.plot(pa = pa.ES , info = uic.ES , code = code, gtf = gtf_mouse, name = "Mouse ES")

suppl = "24"
figureS24 <- ggarrange(S24A$h,S24B$h,S24C$h, as_ggplot(get_legend(S24A$h)),
                       labels = c( "a)", "b)", "c)", ""),
                       ncol = 2, nrow = 2, legend = "none") +
  theme(plot.margin = margin(0.2,0.1,0.2,0.1, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Distribution of Biotypes across pipelines. a) WTC11, c) H1-mix, c) Mouse ES.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_biotype-by-tool",".pdf"))
annotate_figure(figureS24,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.25
########################## 
suppl = "25"
figureS25 <- ggarrange(S24A$q,S24B$q,S24C$q, as_ggplot(get_legend(S24A$q)),
                       labels = c( "a)", "b)", "c)", ""),
                       ncol = 2, nrow = 2, legend = "none") +
  theme(plot.margin = margin(0.2,0.1,0.2,0.1, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Biotypes per pipeline. a) WTC11, c) H1-mix, c) Mouse ES.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_biotype-per-pipeline",".pdf"))
annotate_figure(figureS25,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.26
##########################  Pipelines overalp for H1-mix and WTC11
suppl = "26"
S26A <- agreement.pipelines(data_sample = "H1_mix_results/H1_mix")
S26B <- agreement.pipelines(data_sample = "ES_results/ES")
figureS26 <- ggarrange(S26A ,S26B,
                       labels = c( "a)", "b)"),
                       ncol = 1, nrow = 2, legend = "none") +
  theme(plot.margin = margin(0.2,0.1,0.2,0.1, "cm")) 

mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number and SQANTI category distribution of Unique Intron Chain (UIC) consistently detected by\n     an increasing number of submissions. a) H1-mix sample, b) Mouse ES sample.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_cnt_squanti_dist-UIC",".pdf"))
annotate_figure(figureS26,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig. 27 to 29  Intersection plots
#############################################
suppl = "27"
S27 <- overlap.plot (code = code, data_sample = "Challenge1_Figures_Data/WTC11_results/WTC11")
S28 <- overlap.plot (code = code, data_sample = "Challenge1_Figures_Data/H1_mix_results/H1_mix") 
S29 <- overlap.plot (code = code, data_sample = "Challenge1_Figures_Data/ES_results/ES")

## Extended Data Fig.30
########################## (this is  a time-consuming figure)
subset = c("PacBio", "cDNA")
S30A <- me_and_others(pa = pa.WTC11, code = code, subset = subset, name = "WTC11", remove = c("Iso_IB","Spectra" ), replace = TRUE)
S30B <- me_and_others(pa = pa.H1mix, code = code, subset = subset, name = "H1-mix", remove = c("Iso_IB","Spectra" ), replace = TRUE)
S30C <- me_and_others(pa = pa.ES, code = code, subset = subset, name = "Mouse ES", remove = c("Iso_IB","Spectra" ), replace = TRUE)

suppl = "30"
figureS30 <- ggarrange(S30A$h,S30B$h,S30C$h,
                       labels = c( "a)", "b)", "c)"),
                       ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of UIC detected by a tool and shared with an increasing number of other tools, \n     processing ", paste(subset, collapse = "_"), " data. a) WTC11, c) H1-mix, c) Mouse ES.") 
pdf(paste0(outdir, "/Extended_Fig._",suppl, "_UIC-by-tool-",paste(subset, collapse = "_"), ".pdf"))
annotate_figure(figureS30,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.31
########################## (this is  a time-consuming figure)
subset = c("PacBio", "CapTrap")
S31A <- me_and_others(pa = pa.WTC11, code = code, subset = subset, name = "WTC11", remove = c("Iso_IB","Spectra" ), replace = TRUE)
S31B <- me_and_others(pa = pa.H1mix, code = code, subset = subset, name = "H1-mix", remove = c("Iso_IB","Spectra" ), replace = TRUE)
S31C <- me_and_others(pa = pa.ES, code = code, subset = subset, name = "Mouse ES", remove = c("Iso_IB","Spectra" ), replace = TRUE)

suppl = "31"
figureS31 <- ggarrange(S31A$h,S31B$h,S31C$h,
                       labels = c( "a)", "b)", "c)"),
                       ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of UIC detected by a tool and shared with an increasing number of other tools, \n     processing ", paste(subset, collapse = "_"), " data. a) WTC11, c) H1-mix, c) Mouse ES.") 
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_UIC-by-tool-",paste(subset, collapse = "_"), ".pdf"))
annotate_figure(figureS31,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.32
########################## (this is  a time-consuming figure)
subset = c("ONT", "cDNA")
S32A <- me_and_others(pa = pa.WTC11, code = code, subset = subset, name = "WTC11", remove = c("Iso_IB","Spectra" ), replace = TRUE)
S32B <- me_and_others(pa = pa.H1mix, code = code, subset = subset, name = "H1-mix", remove = c("Iso_IB","Spectra" ), replace = TRUE)
S32C <- me_and_others(pa = pa.ES, code = code, subset = subset, name = "Mouse ES", remove = c("Iso_IB","Spectra" ), replace = TRUE)

suppl = "32"
figureS32 <- ggarrange(S32A$h,S32B$h,S32C$h,
                       labels = c( "a)", "b)", "c)"),
                       ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Biotypes per pipeline. a) WTC11, c) H1-mix, c) Mouse ES.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_UIC-by-tool-",paste(subset, collapse = "_"), ".pdf"))
annotate_figure(figureS32,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.33
########################## (this is  a time-consuming figure)
subset = c("ONT", "CapTrap")
S33A <- me_and_others(pa = pa.WTC11, code = code, subset = subset, name = "WTC11", remove = c("Iso_IB","Spectra" ), replace = TRUE)
S33B <- me_and_others(pa = pa.H1mix, code = code, subset = subset, name = "H1-mix", remove = c("Iso_IB","Spectra" ), replace = TRUE)
S33C <- me_and_others(pa = pa.ES, code = code, subset = subset, name = "Mouse ES", remove = c("Iso_IB","Spectra" ), replace = TRUE)

suppl = "33"
figureS33 <- ggarrange(S33A$h,S33B$h,S33C$h,
                       labels = c( "a)", "b)", "c)"),
                       ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of UIC detected by a tool and shared with an increasing number of other tools, \n     processing ", paste(subset, collapse = "_"), " data. a) WTC11, c) H1-mix, c) Mouse ES.") 
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_UIC-by-tool-",paste(subset, collapse = "_"),".pdf"))
annotate_figure(figureS33,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.34
########################## (this is  a time-consuming figure)
subset = c("ONT", "R2C2")
S34A <- me_and_others(pa = pa.WTC11, code = code, subset = subset, name = "WTC11", remove = c("Iso_IB","Spectra" ), replace = TRUE)
S34B <- me_and_others(pa = pa.H1mix, code = code, subset = subset, name = "H1-mix", remove = c("Iso_IB","Spectra" ), replace = TRUE)
S34C <- me_and_others(pa = pa.ES, code = code, subset = subset, name = "Mouse ES", remove = c("Iso_IB","Spectra" ), replace = TRUE)

suppl = "34"
figureS34 <- ggarrange(S34A$h,S34B$h,S34C$h,
                       labels = c( "a)", "b)", "c)"),
                       ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of UIC detected by a tool and shared with an increasing number of other tools, \n     processing ", paste(subset, collapse = "_"), " data. a) WTC11, c) H1-mix, c) Mouse ES") 
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_UIC-by-tool-",paste(subset, collapse = "_"),".pdf"))
annotate_figure(figureS34,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.35
########################## (this is  a time-consuming figure)
subset = c("ONT", "dRNA")
S35A <- me_and_others(pa = pa.WTC11, code = code, subset = subset, name = "WTC11", remove = c("Iso_IB","Spectra" ), replace = TRUE)
S35B <- me_and_others(pa = pa.H1mix, code = code, subset = subset, name = "H1-mix", remove = c("Iso_IB","Spectra" ), replace = TRUE)
S35C <- me_and_others(pa = pa.ES, code = code, subset = subset, name = "Mouse ES", remove = c("Iso_IB","Spectra" ), replace = TRUE)

suppl = "35"
figureS35 <- ggarrange(S35A$h,S35B$h,S35C$h,
                       labels = c( "a)", "b)", "c)"),
                       ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of UIC detected by a tool and shared with an increasing number of other tools, \n     processing ", paste(subset, collapse = "_"), " data. a) WTC11, c) H1-mix, c) Mouse ES") 
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_UIC-by-tool-",paste(subset, collapse = "_"),".pdf"))
annotate_figure(figureS35,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.36
########################## 
suppl = "36"
figureS36 <- ggarrange(S30A$j,S30B$j,S30C$j,
                       labels = c( "a)", "b)", "c)"),
                       ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Number of UIC consistently detected by a tool across samples.  \n     a) WTC11, c) H1-mix, c) Mouse ES") 
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_UIC-by-tool-sample",".pdf"))
annotate_figure(figureS36,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.37
##########################  # figure is set up in powerpoint from individual figures
S37A <- highly.detected(pa = pa.WTC11, code = code,  name = "WTC11",  replace = TRUE)
S37B <- highly.detected(pa = pa.H1mix, code = code,  name = "H1-mix",  replace = TRUE)
S37C <- highly.detected(pa = pa.ES, code = code,  name = "Mouse ES",  replace = TRUE)

suppl = "37"
figureS37 <- ggarrange(S37A$a,S37A$d, S37B$a,S37B$d, S37C$a,S37C$d,
                       labels = c( "a)", "b)", "c)","d)", "e)", "f)","g)"),
                       ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Characterization of highly detected UICs. /n      
                   a,c,e) Strucutral categories distribution. The table indicates the fold enrichment /n   
                   of each structural category within the frequently detected transcripts respect to their global count. /n
                   b,d,e) Tools identifying highly detected transcripts. The graph shows the enrichment in the number   /n
                       HDT found by a tool with respect to their global number of reported transcripts") 
#pdf(paste0(outdir, "/Extended_Fig._",suppl,".pdf"))
#annotate_figure(figureS37,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
#dev.off()

## Extended Data Fig.38
########################## 
data.38.WTC11 <- trend.analysis (pa = pa.WTC11, info = uic.WTC11, code = code)
data.38.ESmouse <- trend.analysis (pa = pa.ES, info = uic.ES, code = code)
data.38.H1mix <- trend.analysis (pa = pa.H1mix, info = uic.H1mix, code = code)

S38AW <- comparisons.plot(bp.data = data.38.WTC11[[1]],name = "WTC11", my_comparisons = NULL)
S38AH <- comparisons.plot(bp.data = data.38.H1mix[[1]],name = "H1-mix", my_comparisons = NULL)
S38AE <- comparisons.plot(bp.data = data.38.ESmouse[[1]],name = "Mouse ES", my_comparisons = NULL)

suppl = "38"
figureS38 <- ggarrange(S38AW[[1]], S38AW[[2]],S38AW[[3]],
                       S38AH[[1]], S38AH[[2]],S38AH[[3]],
                       S38AE[[1]], S38AE[[2]],S38AE[[3]],
                       labels = c( "a)", "b)", "c)","d)", "e)", "f)","g)", "h)", "i)"),
                       ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Properties of detected transcripts by library preparation. \n     a,d,g) Length distribution. b,e,h) Exon number distribution. c,f,i) Counts per million")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_trans-prop-by-prep",".pdf"))
annotate_figure(figureS38,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()  

## Extended Data Fig.39
##########################
S39AW <- comparisons.plot(bp.data = data.38.WTC11[[2]],name = "WTC11", my_comparisons = NULL)
S39AH <- comparisons.plot(bp.data = data.38.H1mix[[2]],name = "H1-mix", my_comparisons = NULL)
S39AE <- comparisons.plot(bp.data = data.38.ESmouse[[2]],name = "Mouse ES", my_comparisons = NULL)

suppl = "39"
figureS39 <- ggarrange(S39AW[[1]], S39AW[[2]],S39AW[[3]],
                       S39AH[[1]], S39AH[[2]],S39AH[[3]],
                       S39AE[[1]], S39AE[[2]],S39AE[[3]],
                       labels = c( "a)", "b)", "c)","d)", "e)", "f)","g)", "h)", "i)"),
                       ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Properties of detected transcripts by platform. \n     a,d,g) Length distribution. b,e,h) Exon number distribution. c,f,i) Counts per million")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_trans-prop-by-plat",".pdf"))
annotate_figure(figureS39,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.40
##########################
S40AW <- comparisons.plot(bp.data = data.38.WTC11[[3]],name = "WTC11", my_comparisons = NULL)
S40AH <- comparisons.plot(bp.data = data.38.H1mix[[3]],name = "H1-mix", my_comparisons = NULL)
S40AE <- comparisons.plot(bp.data = data.38.ESmouse[[3]],name = "Mouse ES", my_comparisons = NULL)

suppl = "40"
figureS40 <- ggarrange(S40AW[[1]], S40AW[[2]],S40AW[[3]],
                       S40AH[[1]], S40AH[[2]],S40AH[[3]],
                       S40AE[[1]], S40AE[[2]],S40AE[[3]],
                       labels = c( "a)", "b)", "c)","d)", "e)", "f)","g)", "h)", "i)"),
                       ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Properties of detected transcripts by experimental protocol. \n     a,d,g) Length distribution. b,e,h) Exon number distribution. c,f,i) Counts per million")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_trans-by-protocol",".pdf"))
annotate_figure(figureS40,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

#  Extended Data Fig. 41. Transcript_length distribution by analysis pipeline
###############################################################################

S41A <- length_pipelines (pa.WTC11, uic.WTC11, name = "WTC11")
S41B <- length_pipelines (pa.H1mix, uic.H1mix, name = "H1-mix")
S41C <- length_pipelines (pa.ES, uic.ES, name = "Mouse ES")

suppl = "41"
figureS41 <- ggarrange(S41A, S41B, S41C,
                       labels = c( "a)", "b)", "c)"),
                       ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Distribution of transcript length by analysis tool.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_trans-len-by-tool",".pdf"))
annotate_figure(figureS41,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.42
########################### coverage SIRVs

## Extended Data Fig.43
##########################
suppl = "43"
text = "Extended Data Fig.43. Performance metrics on mouse simulated data. Sen_kn: sensitivity known transcripts, \nSen_kn > 5TMP: sensitivity known transcripts with expression > 5 TPM, Pre_kn: precision known transcripts,
Sen_no: sensitivity novel transcripts, Pre_no: precision novel transcripts, 1/Red: inverse of redundancy."
radar.simulation (species = "mouse", directory = "Challenge1_Figures_Data/Simulations/",
                  pdf = paste0(outdir, "/Extended_Fig._43", "_metrics-mouse-simul"), text = text, a = 0.05)

## Extended Data Fig.44. Coverage simulated data
####################################################
df1 <- as.data.frame(t(sim.data[1:3, 3:102])) ; colnames(df1) <- sim.data[1:3,2] ; 
df1 <- data.frame(transcript_position = rep(df1[,1],2), coverage = c(df1[,2], df1[,3]), type_of_data = c(rep("Real", nrow(df1)), rep("Simulated", nrow(df1))))

df2 <- as.data.frame(t(sim.data[c(1,4,5), 3:102])) ; colnames(df2) <- sim.data[c(1,4,5),2] ; 
df2 <- data.frame(transcript_position = rep(df2[,1],2), coverage = c(df2[,2], df2[,3]), type_of_data = c(rep("Real", nrow(df2)), rep("Simulated", nrow(df2))))

df3 <- as.data.frame(t(sim.data[c(1,6,7), 3:102])) ; colnames(df3) <- sim.data[c(1,6,7),2] ; 
df3 <- data.frame(transcript_position = rep(df3[,1],2), coverage = c(df3[,2], df3[,3]), type_of_data = c(rep("Real", nrow(df3)), rep("Simulated", nrow(df3))))

df4 <- as.data.frame(t(sim.data[c(1,8,9), 3:102])) ; colnames(df4) <- sim.data[c(1,8,9),2] ; 
df4 <- data.frame(transcript_position = rep(df4[,1],2), coverage = c(df4[,2], df4[,3]), type_of_data = c(rep("Real", nrow(df4)), rep("Simulated", nrow(df4))))

S44A <- line.plot(df1, title = sim.data[2,1])
S44B <- line.plot(df2, title = sim.data[4,1])
S44C <- line.plot(df3, title = sim.data[6,1])
S44D <- line.plot(df4, title = sim.data[8,1])

suppl = "44"
figureS44 <- ggarrange(S44A, S44B, S44C, S44D,
                       labels = c( "a)", "b)", "c)","d)"),
                       ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) 
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Comparison of long-read transcript coverage between real and simulated datasets.")
pdf(paste0(outdir, "/Extended_Fig._",suppl, "_trans-cover-real-simul",".pdf"))
annotate_figure(figureS44,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.45. Figure GENCODE data for human
#########################################################
S45 <- gencode.analysis ( gencode.self, uic.WTC11, name = "WTC11")

suppl = "45"
figureS45 <- ggarrange(ggarrange(S45$A, S45$B, ncol = 2, labels = c("a)", "b)")),S45$C,
          nrow = 2, 
          labels = c("","c)")                                        
) 

mylegend <- paste0("     Extended Data Fig. ", suppl, ". Properties of GENCODE manually annotated loci for WTC11 sample.a) Distributon of gene \n     expression. b) Distribution of SQANTI categories. c) Intersection of Unique Intron Chains (UIC) among experimental \n     protocols.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_gencode-manual-props-wtc11",".pdf"))
annotate_figure(figureS45,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.46. Figure GENCODE data for mouse
#########################################################
S46 <- gencode.analysis(gencode.self_mouse, uic.ES, name = "mouse ES")

suppl = "46"
figureS46 <- ggarrange(ggarrange(S46$A, S46$B, ncol = 2, labels = c("a)", "b)")),S46$C, # Second row with box and dot plots
                       nrow = 2, 
                       labels = c("","c)")                                        # Labels of the scatter plot
) 

mylegend <- paste0("     Extended Data Fig. ", suppl, ". Properties of GENCODE manually annotated loci for mouse ES sample.a) Distributon of gene \n     expression. b) Distribution of SQANTI categories. c) Intersection of Unique Intron Chains (UIC) among experimental \n     protocols.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_gencode-manual-props-es",".pdf"))
annotate_figure(figureS46,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.47. Evaluation against GENCODE_mouse
############################################################
pa_GENCODE_mouse = read.csv("Challenge1_Figures_Data/GENCODE_manualAnnot/GENCODE_mouse/presence_absence.GENCODE_loci.csv", sep=",", header = T, as.is = TRUE) [,1:51] 
pa_GENCODE_mouse[,1] <- sapply(pa_GENCODE_mouse[,1] , function (x) paste(strsplit(x, split = "_")[[1]][-1], collapse= "_"))

genocode_eval_mouse <- performance.genecode (gencode.pa = pa_GENCODE_mouse, ID_UIC = NULL,
                                             pa = pa.ES,  code = code, selection = NULL,
                                             mypattern = "SQ3_mouse", evaluation = gencode_eval_results_mouse,
                                             directory = "Challenge1_Figures_Data/GENCODE_manualAnnot/classifications/mouse/")

pivoted_gencode_gene_M <- pivot_longer(genocode_eval_mouse, cols = c("Sensitivity.Genes", "Precision.Genes", "F1_score.Genes"))
pivoted_gencode_gene_M$name <- pivoted_gencode_gene_M$name %>% factor(levels = c("Sensitivity.Genes", "Precision.Genes", "F1_score.Genes"),
                                                                  labels = c("Sensitivity", "Precision", "F1-score"))

pivoted_gencode_known_M <- pivot_longer(genocode_eval_mouse , cols = c("Sensitivity_known", "Precision_known", "F1_known"))
pivoted_gencode_known_M$name <- pivoted_gencode_known_M$name %>% factor(levels = c("Sensitivity_known", "Precision_known", "F1_known"),
                                                                    labels = c("Sensitivity", "Precision", "F1-score"))

pivoted_gencode_novel_M <- pivot_longer(genocode_eval_mouse , cols = c("Sensitivity_novel", "Precision_novel", "F1_novel"))
pivoted_gencode_novel_M$name <- pivoted_gencode_novel_M$name %>% factor(levels = c("Sensitivity_novel", "Precision_novel", "F1_novel"),
                                                                    labels = c("Sensitivity", "Precision", "F1-score"))

genocode_eval_mouse$False_known <- genocode_eval_mouse$Transcript_models_known - genocode_eval_mouse$TRUE_known
genocode_eval_mouse$False_novel <- genocode_eval_mouse$Transcript_models_novel - genocode_eval_mouse$TRUE_novel
pivoted_gencode_M_TP <- pivot_longer(genocode_eval_mouse , cols = c("TRUE_known", "TRUE_novel", "False_known", "False_novel"))
pivoted_gencode_M_TP <- pivoted_gencode_M_TP[pivoted_gencode_M_TP$value > 0,]
pivoted_gencode_M_TP$name <- pivoted_gencode_M_TP$name %>% factor(levels = c("TRUE_known","TRUE_novel",
                                                                         "False_known", "False_novel"),
                                                              labels = c("TRUE\nknown", "TRUE\nnovel",
                                                                         "FALSE\nknown", "FALSE\nnovel"))



#pC.gene.left <- Performance_plot_left(pivoted_gencode_gene_M , main = "Gene level", axis.text.x = 7)
pC.gene.left <- ggplot(pivoted_gencode_gene_M, aes(x=Label, y=value)) +
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


#pC.known.mid <- Performance_plot_middel(pivoted_gencode_known_M , main = "Known_transcript level", axis.text.x = 7)
pC.known.mid <- ggplot(pivoted_gencode_known_M, aes(x=Label, y=value)) +
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

#pC.novel <- Performance_plot_right(pivoted_gencode_novel_M , main = "Novel_transcript level", axis.text.x = 7)
pC.novel <- ggplot(pivoted_gencode_novel_M, aes(x=Label, y=value)) +
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

#pC.TP <- Performance_plot_TP(pivoted_gencode_M_TP , main = "Number detected transcripts" )
pC.TP <- ggplot(pivoted_gencode_M_TP, aes(x=Label, y=value)) +
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


suppl = "47"
figureS47 <- ggarrange(pC.gene.left, pC.known.mid,
                       pC.novel, pC.TP,
                       labels = c("a)", "b)", "c)", "d)"),
                       ncol = 2, nrow = 2, common.legend = TRUE)
  
  #ggarrange(pC.gene.left, pC.known.mid,  pC.novel, labels = c("", "", ""),
  #                   ncol = 3, nrow = 1, common.legend = TRUE)

mylegend <- paste0("     Extended Data Fig. ", suppl, ". Performance metrics of LRGASP pipelines evaluate against GENCODE manual annotation \n     of mouse ES sample. Ba: Bambu, FM: Flames, FR: FLAIR, IQ: IsoQuant, IT: IsoTools, IB: Iso_IB, Ly: LyRic, \n     Ma: Mandalorion, TL: TALON-LAPA, Sp: Spectra, ST: StringTie2."  )
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_pipeline-by-gencode-manual",".pdf"), width=20, height = 7)
print(figureS47)
dev.off()

pdf(paste0(outdir, "/Extended_Fig._31_pipeline-by-gencode-manual-Genes",".pdf"), width=10, height = 5)
print(pC.gene.left) 
dev.off()

pdf(paste0(outdir, "/Extended_Fig._31_pipeline-by-gencode-manual-KnownTrx",".pdf"), width=10, height = 5)
print(pC.known.mid)
dev.off()

pdf(paste0(outdir, "/Extended_Fig._31_pipeline-by-gencode-manual-Novel",".pdf"), width=10, height = 5)
print(pC.novel)
dev.off()

pdf(paste0(outdir, "/Extended_Fig._31_pipeline-by-gencode-manual-NumDetected",".pdf"), width=10, height = 5)
print(pC.TP)
dev.off()

## Extended Data Fig.48. Evaluation against GENCODE_mouse TP
###############################################################
pC.TP <- Performance_plot_TP(pivoted_gencode_M_TP , main = "Number of detected transcripts. Mouse ES" )

suppl = "48"
mylegend <- paste0("     Extended Data Fig. ", suppl, ". Detection of Unique Intron Chains (UIC) at GENCODE manual annotation loci. Ba: Bambu, \n     FM: Flames, FL: FLAIR, IQ: IsoQuant, IT: IsoTools, IB: Iso_IB, Ly: LyRic, Ma: Mandalorion, TL: TALON-LAPA,\n     Sp: Spectra, ST: StringTie2."  )
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_uic-gencode-manual-by-tool",".pdf"))
annotate_figure(pC.TP,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()  # for some reason the pdf includes one empty page. Delete from pdf itself.


## Extended Data Fig.S49. Selecting GENCODE transcripts found in more than one sample
########################################################################################
# Identify and select the transcript that are in more than one sample looking at the GENCODE data self
gencode.self2 <- gencode.self[-c(1:2),c(1:9)]
#gencode.self2[,1] <- sapply(gencode.self2[,1], function (x) paste(strsplit(x, split = "_")[[1]][-1], collapse= "_")) 
gencode.self3 <- matrix(as.numeric(as.matrix(gencode.self2[,c(2:7)])), nrow= 271, ncol = 6)
#gencode.self4 <- data.frame(dRNA_ONT = tapply(gencode.self3[,1],gencode.self2[,1], sum ), 
#                            cDNA_ONT = tapply(gencode.self3[,2],gencode.self2[,1], sum ),
#                            CapTrap_ONT = tapply(gencode.self3[,3],gencode.self2[,1], sum ),
#                            R2C2_ONT = tapply(gencode.self3[,4],gencode.self2[,1], sum ),
#                            cDNA_PacBio = tapply(gencode.self3[,5],gencode.self2[,1], sum ),
#                            CapTrap_PacBio = tapply(gencode.self3[,6],gencode.self2[,1], sum ))
colnames(gencode.self3) <- colnames(gencode.self2)[c(2:7)]

selection <- gencode.self2[rowSums(gencode.self3) > 1,1] # UIC of manually annotated transcripts that are in more than one sample

genocode_eval_WTC11.more_samples <- performance.genecode (gencode.pa = pa_GENCODE, ID_UIC = ID_UIC,
                                                          pa = pa.WTC11,  code = code, selection = selection,
                                                          evaluation = gencode_eval_results,
                                                          mypattern = "SQ3_human",
                                                          directory = "Challenge1_Figures_Data/GENCODE_manualAnnot/classifications/human/")

pivoted_gencode_gene_S <- pivot_longer(genocode_eval_WTC11.more_samples, cols = c("Sensitivity.Genes", "Precision.Genes", "F1_score.Genes"))
pivoted_gencode_known_S <- pivot_longer(genocode_eval_WTC11.more_samples, cols = c("Sensitivity_known", "Precision_known", "F1_known"))
pivoted_gencode_novel_S <- pivot_longer(genocode_eval_WTC11.more_samples, cols = c("Sensitivity_novel", "Precision_novel", "F1_novel"))

genocode_eval_WTC11.more_samples$False_known <- genocode_eval_WTC11.more_samples$Transcript_models_known - genocode_eval_WTC11.more_samples$TRUE_known
genocode_eval_WTC11.more_samples$False_novel <- genocode_eval_WTC11.more_samples$Transcript_models_novel - genocode_eval_WTC11.more_samples$TRUE_novel
pivoted_gencode_S_TP <- pivot_longer(genocode_eval_WTC11.more_samples , cols = c("TRUE_known", "TRUE_novel", "False_known", "False_novel"))

pC.known.S <- Performance_plot_left(pivoted_gencode_known_S , main = "Known_transcript level", axis.text.x = 7)
pC.novel.S <- Performance_plot_right(pivoted_gencode_novel_S , main = "Novel_transcript level", axis.text.x = 7)

suppl = "49"
figureS49 <- ggarrange(pC.known.S, pC.novel.S, labels = c("", ""),
                      ncol = 2, nrow = 1, common.legend = TRUE)

mylegend <- paste0("     Extended Data Fig. ", suppl, ". Performance on GENCODE manually curated data. Curated transcripts selected to be present \n     in at least two experimental datasets. Ba: Bambu, FM: Flames, FL: FLAIR, IQ: IsoQuant, IT: IsoTools, IB: Iso_IB,  \n     Ly: LyRic, Ma: Mandalorion, TL: TALON-LAPA, Sp: Spectra, ST: StringTie2.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_gencode-manual-by-tool-two-datasets",".pdf"))
annotate_figure(figureS49,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.S50. Selecting GENCODE transcripts found with a CPM filter (WTC11)
########################################################################################

merged.gencode <- merge(gencode.self2, ID_UIC, by.x = "Row.names", by.y = "LRGASP_id")
merged.gencode2 <- merge (merged.gencode , gencode_human, by.x = "isoform", by.y = "Transcript Id")
merged.gencode3 <- unique(merged.gencode2 [,c("Row.names", "Read count (full-length)", "Library")])
merged.gencode3.count <- tapply(merged.gencode3$'Read count (full-length)', merged.gencode3$Row.names, sum)
select.count <- names(merged.gencode3.count) [merged.gencode3.count > 2]

genocode_eval_WTC11.more_reads <- performance.genecode (gencode.pa = pa_GENCODE, ID_UIC = ID_UIC,
                                                         pa = pa.WTC11,  code = code, selection = select.count ,
                                                         evaluation = gencode_eval_results,
                                                         mypattern = "SQ3_human",
                                                         directory = "Challenge1_Figures_Data/GENCODE_manualAnnot/classifications/human/")

pivoted_gencode_gene_R <- pivot_longer(genocode_eval_WTC11.more_reads, cols = c("Sensitivity.Genes", "Precision.Genes", "F1_score.Genes"))
pivoted_gencode_known_R <- pivot_longer(genocode_eval_WTC11.more_reads, cols = c("Sensitivity_known", "Precision_known", "F1_known"))
pivoted_gencode_novel_R <- pivot_longer(genocode_eval_WTC11.more_reads, cols = c("Sensitivity_novel", "Precision_novel", "F1_novel"))

genocode_eval_WTC11.more_reads$FALSE_known <- genocode_eval_WTC11.more_reads$Transcript_models_known - genocode_eval_WTC11.more_reads$TRUE_known
genocode_eval_WTC11.more_reads$FALSE_novel <- genocode_eval_WTC11.more_reads$Transcript_models_novel - genocode_eval_WTC11.more_reads$TRUE_novel
pivoted_gencode_R_TP <- pivot_longer(genocode_eval_WTC11.more_reads , cols = c("TRUE_known", "TRUE_novel", "FALSE_known", "FALSE_novel"))

pC.known.R <- Performance_plot_left(pivoted_gencode_known_R , main = "Known_transcript level", axis.text.x = 7)
pC.novel.R <- Performance_plot_right(pivoted_gencode_novel_R , main = "Novel_transcript level", axis.text.x = 7)


suppl = "50"
figureS50 <- ggarrange(pC.known.R , pC.novel.R, labels = c("", ""),
                       ncol = 2, nrow = 1, common.legend = TRUE)

mylegend <- paste0("     Extended Data Fig. ", suppl, ". Performance on GENCODE manually curated data. The ground truth is the set of manually\n     annotated transcripts with more than two reads. Ba: Bambu, FM: Flames, FL: FLAIR, IQ: IsoQuant, IT: IsoTools,\n     IB: Iso_IB, Ly: LyRic, Ma: Mandalorion, TL: TALON-LAPA, Sp: Spectra, ST: StringTie2.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_gencode-manual-by-tool-three-reads",".pdf"))
annotate_figure(figureS50,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.51. Boxplox summary performance for Library preparations (WTC11)
#########################################################################################

genocode_eval_WTC11 <- performance.genecode (gencode.pa = pa_GENCODE, ID_UIC = NULL,
                                             pa = pa.WTC11,  code = code, 
                                             selection = NULL, evaluation = gencode_eval_results,
                                             mypattern = "SQ3_human",
                                             directory = "Challenge1_Figures_Data/GENCODE_manualAnnot/classifications/human/")


all.labs <- unique(code$Lab)
Lab.colors = colorConesa(n = length(all.labs)) ; names (Lab.colors) = all.labs

A <- fig.boxplots (genocode_eval_WTC11, sample = "WTC11", var.y = "Library_Preps", var.x = "Precision_known", xlabel = "Precision_known", jitter.color = "Lab", mycolor = Lab.colors)
B <- fig.boxplots (genocode_eval_WTC11, sample = "WTC11", var.y = "Library_Preps", var.x = "Precision_novel", xlabel = "Precision_novel", jitter.color = "Lab", mycolor = Lab.colors)
B <- B + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
C <- fig.boxplots (genocode_eval_WTC11, sample = "WTC11", var.y = "Library_Preps", var.x = "Sensitivity_known", xlabel = "Sensitivity_known", jitter.color = "Lab", mycolor = Lab.colors)
C <- C + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
D <- fig.boxplots (genocode_eval_WTC11, sample = "WTC11", var.y = "Library_Preps", var.x = "Sensitivity_novel", xlabel = "Sensitivity_novel", jitter.color = "Lab", mycolor = Lab.colors)
D <- D + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


A2 <- fig.boxplots (genocode_eval_mouse, sample = "Mouse ES", var.y = "Library_Preps", var.x = "Precision_known", xlabel = "Precision_known", jitter.color = "Lab", mycolor = Lab.colors)
B2 <- fig.boxplots (genocode_eval_mouse, sample = "Mouse ES", var.y = "Library_Preps", var.x = "Precision_novel", xlabel = "Precision_novel", jitter.color = "Lab", mycolor = Lab.colors)
B2 <- B2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
C2 <- fig.boxplots (genocode_eval_mouse, sample = "Mouse ES", var.y = "Library_Preps", var.x = "Sensitivity_known", xlabel = "Sensitivity_known", jitter.color = "Lab", mycolor = Lab.colors)
C2 <- C2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
D2 <- fig.boxplots (genocode_eval_mouse, sample = "Mouse ES", var.y = "Library_Preps", var.x = "Sensitivity_novel", xlabel = "Sensitivity_novel", jitter.color = "Lab", mycolor = Lab.colors)
D2 <- D2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

suppl = "51"
figureS51 <- ggarrange(A,B,C,D,A2,B2,C2,D2, labels = c("a)", "b)", "c)", "d)", "e)", "f)", "g)"),
                       ncol = 4, nrow = 2, common.legend = TRUE)

mylegend <- paste0("     Extended Data Fig. ", suppl, ". Performance on GENCODE manually curated data by Library Preparation.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_gencode-manual-by-prep",".pdf"))
annotate_figure(figureS51,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()

## Extended Data Fig.52. Boxplox summary performance for sequencing platforms (WTC11)

A <- fig.boxplots (genocode_eval_WTC11, sample = "WTC11", var.y = "Platform", var.x = "Precision_known", xlabel = "Precision_known", jitter.color = "Lab", mycolor = Lab.colors)
B <- fig.boxplots (genocode_eval_WTC11, sample = "WTC11", var.y = "Platform", var.x = "Precision_novel", xlabel = "Precision_novel", jitter.color = "Lab", mycolor = Lab.colors)
B <- B + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
C <- fig.boxplots (genocode_eval_WTC11, sample = "WTC11", var.y = "Platform", var.x = "Sensitivity_known", xlabel = "Sensitivty_known", jitter.color = "Lab", mycolor = Lab.colors)
C <- C + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
D <- fig.boxplots (genocode_eval_WTC11, sample = "WTC11", var.y = "Platform", var.x = "Sensitivity_novel", xlabel = "Sensitivty_novel", jitter.color = "Lab", mycolor = Lab.colors)
D <- D + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

A2 <- fig.boxplots (genocode_eval_mouse, sample = "Mouse ES", var.y = "Platform", var.x = "Precision_known", xlabel = "Precision_known", jitter.color = "Lab", mycolor = Lab.colors)
B2 <- fig.boxplots (genocode_eval_mouse, sample = "Mouse ES", var.y = "Platform", var.x = "Precision_novel", xlabel = "Precision_novel", jitter.color = "Lab", mycolor = Lab.colors)
B2 <- B2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
C2 <- fig.boxplots (genocode_eval_mouse, sample = "Mouse ES", var.y = "Platform", var.x = "Sensitivity_known", xlabel = "Sensitivity_known", jitter.color = "Lab", mycolor = Lab.colors)
C2 <- C2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
D2 <- fig.boxplots (genocode_eval_mouse, sample = "Mouse ES", var.y = "Platform", var.x = "Sensitivity_novel", xlabel = "Sensitivity_novel", jitter.color = "Lab", mycolor = Lab.colors)
D2 <-D2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

suppl = "52"
figureS52 <- ggarrange(A,B,C,D,A2,B2,C2,D2, labels = c("a)", "b)", "c)", "d)", "e)", "f)", "g)"),
                       ncol = 4, nrow = 2, common.legend = TRUE)

mylegend <- paste0("     Extended Data Fig. ", suppl, ". Performance on GENCODE manually curated data by Platform.")
pdf(paste0(outdir, "/Extended_Fig._",suppl,"_gencode-manual-by-plat",".pdf"))
annotate_figure(figureS52,  bottom = text_grob(mylegend, hjust = 0,  x = 0,  size = 9))
dev.off()


# Local Variables:
# ess-indent-offset: 2
# End:
