rm(list=ls())
library(data.table)
library(tidyverse)
options(bitmapType="cairo")

setwd("/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/")

gtex = data.table()

for (id in seq(30000,30300, 10)){
  
  for (i in c(1:22)){
    
    if(file.exists(paste0("results/Coloc_results_V2/GTEx/",id,"/chr",i,"_coloc_results.txt"))){
      
      gtex_chr = fread(paste0("results/Coloc_results_V2/GTEx/",id,"/chr",i,"_coloc_results.txt"))
      gtex = rbind(gtex, gtex_chr)
      
    } else{
      next()
    }
    
  }
  
}

gtex = gtex %>%
  filter(pval_nominal < 1e-03, pval < 1e-05)

interval = data.table()

for (id in seq(30000,30300, 10)){
  
  for (i in c(1:22)){
    
    if(file.exists(paste0("results/Coloc_results_V2/Interval/",id,"/chr",i,"_coloc_results.txt"))){
      
      interval_chr = fread(paste0("results/Coloc_results_V2/Interval/",id,"/chr",i,"_coloc_results.txt"))
      interval = rbind(interval, interval_chr, fill=T)
      
    } else{
      next()
    }
    
  }
  
}

interval = interval %>%
  filter(pval_nominal < 1e-03, pval < 1e-05)

# barplot comparing number of colocalizations per GWAS
gtex_barplot = gtex %>%
  group_by(gwas) %>%
  summarise(n = n()) %>%
  drop_na()  %>% mutate(study = "GTEx")

interval_barplot = interval %>%
  group_by(gwas) %>%
  summarise(n = n()) %>%
  drop_na() %>% mutate(study = "interval")


df = rbind(gtex_barplot, interval_barplot)

ggplot(df, aes(x = gwas, y = n, color = study)) +
  geom_bar(stat = "identity",  position = "dodge")  +
  coord_flip()

library(dplyr)
library(ggplot2)
library(forcats)
library(scales)
library(grid)   # for unit()

# 1) Order GWAS by total count across studies (largest at top)
df_plot <- df %>%
  group_by(gwas) %>%
  mutate(total_n = sum(n, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(gwas = fct_reorder(gwas, total_n, .desc = TRUE))

# 2) Colorblind-safe palette (Okabe–Ito; extend if you have >8 studies)
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

# 3) Positioning for nice dodge and label placement
pd <- position_dodge2(width = 0.72, padding = 0.05, preserve = "single")

png("Post_coloc_analysis/plots/GTEx_vs_Interval_barplot.png", width = 6, height = 4, units = "in", res = 300)

# 4) Plot
# rebuild the plot using fully-qualified scales calls
p <- ggplot(df_plot, aes(x = gwas, y = n, fill = study)) +
  geom_col(position = pd, width = 0.64) +
  geom_text(
    aes(label = scales::comma(n)),
    position = position_dodge2(width = 0.72, padding = 0.05),
    hjust = -0.1, size = 2.8
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = okabe_ito) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    labels = scales::label_comma()
  ) +
  guides(fill = guide_legend(title = "Study", nrow = 1, byrow = TRUE)) +
  labs(x = NULL, y = "Number of loci") +
  theme_classic(base_size = 9) +
  theme(
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.length = unit(2, "pt"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(margin = margin(r = 4)),
    axis.text.x = element_text(margin = margin(t = 2)),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.height = unit(6, "pt"),
    legend.key.width  = unit(12, "pt"),
    plot.margin = margin(5.5, 14, 5.5, 5.5)
  )

p
dev.off()

