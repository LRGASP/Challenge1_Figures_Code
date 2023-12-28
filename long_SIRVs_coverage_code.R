library(ggplot2)
library(scales)
library(ggpubr)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)

outdir = "output/sirvs"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)


pub_theme <- theme_pubclean(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=8),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=8) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=10), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15.5)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  theme(legend.position = "bottom")

### keep it simple

depths <- list()
for (i in c("PacBio_cDNA", "PacBio_CapTrap","ONT_cDNA", "ONT_CapTrap", "ONT_R2C2","ONT_dRNA")){
  depths[[i]] <- read.csv(paste0("Challenge1_Figures_Data/SIRV-coverage/",i,".depth.tsv"), sep="\t", header=F)
}


plot_depths <- function(x){
  x <- x %>% mutate(no_depth=ifelse(V3==0,
                                    "Not covered",
                                    "Covered"))
  x$V1 <- factor(x$V1 , levels=c("SIRV4001", "SIRV4002", "SIRV4003",
                                 "SIRV6001", "SIRV6002", "SIRV6003",
                                 "SIRV8001", "SIRV8002", "SIRV8003",
                                 "SIRV10001", "SIRV10002", "SIRV10003",
                                 "SIRV12001", "SIRV12002", "SIRV12003"))
  ggplot(x, aes(x=V2, y=V3, color = no_depth, group=1))+
    geom_line() +
    facet_wrap(~V1, ncol = 3, scales = "free") +
    scale_color_manual(values = c("#6BAED6", "#EE6A50")) +
    labs(y="Coverage", x="length (kb)") +
    scale_x_continuous(label = unit_format(unit = "", scale = 0.001))+
    scale_y_continuous(label = unit_format(unit = "x", scale = 1))+
    pub_theme
}



plots_depths <- map(depths, plot_depths)

for (j in c("PacBio_cDNA", "PacBio_CapTrap","ONT_cDNA", "ONT_CapTrap", "ONT_R2C2","ONT_dRNA")){
  pdf(paste0(outdir, "/", "sirv-cover-", j, ".pdf"))
  print(plots_depths[[j]]+ggtitle(j))
  dev.off()
}


#### Idea to put all in the same plot

get_sirv_length <- function(x){
  x$type <- apply(x,1,function(y){str_split(y["V1"], pattern="00[1-3]", simplify = T)[1]})
  x$length <- apply(x,1,function(y){str_split(y["type"], pattern="SIRV", simplify = T)[2]})
  x$sirv_length <- paste0(x$length, "kb SIRVs")
  x
}

new_depths <- map(depths, get_sirv_length)

### make it step wise across the length of the SIRV

calculate_mean_cov <- function(df_depth){
  df_depth$breaks=cut(df_depth$V2, breaks=20)
  df_depth <- df_depth %>% group_by(breaks) %>% mutate(mean_cov=mean(V3))
  df_depth
}

mean_cov_by_sirv <- function(df_all_depth){
  list_sirvs <- split(df_all_depth, f=~V1)
  new_list <- map(list_sirvs,calculate_mean_cov)
  bind_rows(new_list, .id = "V1")
}

new_depths <- map(new_depths, mean_cov_by_sirv)

all_depths <- bind_rows(new_depths, .id = "Data_type")
all_depths$type <- factor(all_depths$type, levels=c("SIRV4", "SIRV6", "SIRV8","SIRV10","SIRV12"))

library(MetBrewer)
library(RColorConesa)

p.all_sirvs.mean <- ggplot(all_depths, aes(x=V2, y=mean_cov, color = V1, group=V1))+
  geom_line() +
  facet_wrap(facets=Data_type~type, scales = "free", ncol = 5, ) +
  scale_color_conesa(palette = "complete", reverse = F, name="Spike-ins") +
  labs(y="Coverage", x="SIRV length (kb)") +
  scale_x_continuous(label = unit_format(unit = "", scale = 0.001,
                                         expand=expansion(mult=c(0,0))))+
  scale_y_continuous(label = unit_format(unit = "x", scale = 1,
                                         expand=expansion(mult=c(0,0))))+
  pub_theme +
  theme(strip.text = element_blank(), 
        legend.position = "bottom")+
  guides(color=guide_legend(nrow=2))

p.all_sirvs.cov <- ggplot(all_depths, aes(x=V2, y=V3, color = V1, group=V1))+
  geom_line() +
  facet_wrap(facets=Data_type~type, scales = "free", ncol = 5, ) +
  scale_color_conesa(palette = "complete", reverse = F, name="Spike-ins") +
  labs(y="Coverage", x="SIRV length (kb)") +
  scale_x_continuous(label = unit_format(unit = "", scale = 0.001,
                                         expand=expansion(mult=c(0,0))))+
  scale_y_continuous(label = unit_format(unit = "x", scale = 1,
                                         expand=expansion(mult=c(0,0))))+
  pub_theme +
  theme(strip.text = element_blank(), 
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=2))


ggsave(paste0(outdir, "/", "Coverage_all_SIRVs.steps.svg"), p.all_sirvs.mean, width = 12, height = 7)
ggsave(paste0(outdir, "/", "Coverage_all_SIRVs.cov.svg"), p.all_sirvs.cov, width = 12, height = 7)
