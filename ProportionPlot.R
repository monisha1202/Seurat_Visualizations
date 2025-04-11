#Proportion Plot to find how are cells distributed within groups

#Load packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)


seurat_obj <- readRDS("seurat_obj.RDS")


# Proportion / cell number composition per cluster in each Group
ggData = data.frame(prop.table(table(seurat_obj$Group, 
                                     seurat_obj$CellName), 
                               margin = 2))
colnames(ggData) = c("library", "cluster", "value")

#Viewing the table
View(ggData)

#Making plot
prop_plot <- ggplot(ggData, aes(cluster, value, fill = library)) +
  geom_col() +
  ylab("Proportion of Cells (%)") +
  theme(axis.title.y=element_blank()) +
  coord_flip() +
  ggtitle("Propotion in seurat_obj")+
  theme(axis.title.y=element_blank())

#Print plot
View(prop_plot)
