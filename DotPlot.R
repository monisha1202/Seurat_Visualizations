#Load packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)


seurat_obj <- readRDS("seurat_obj.RDS")
gene_list <- c("gene1","gene2","gene3")

#DotPlot
DotPlot(seurat_obj, features = gene_list) 

