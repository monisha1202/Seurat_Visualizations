#Load packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)


seurat_obj <- readRDS("seurat_obj.RDS")
gene_list <- c("gene1","gene2","gene3")

#DotPlot with average expression and percentage expressed bars at bottom

Idents(seurat_obj) <- "CellName"

plot <- DotPlot(seurat_obj, features = gene_list, cluster.idents = F)+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  theme(text = element_text(size=10,face="bold"))+
  guides(colour = guide_colorbar(order=1, ticks.colour = NA, barheight = 1, barwidth = 6,
                                 title.position="top", title.hjust = 0.5),
         size = guide_legend(order=2, reverse=F, title.position="top", title.hjust = 0.5,
                             label.position = "bottom")) +
  scale_color_gradientn(name = "Avg Exp", colours = rev(brewer.pal(n = 8, name = "RdBu"))) +
  scale_size_continuous(name = "% Exp") +
  theme(axis.title = element_blank(),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.spacing = unit(1, "points"),
        legend.box.spacing = unit(0.75, "points"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing.x = unit(0.1, 'points'), legend.justification = c("left")) + coord_flip()

#View Plot
View(plot)

#Save figure
ggsave(file="plot.svg", plot=plot, width=10, height=8) #custom size
