#ClusterProfiler

#Load packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

#Load seurat object
seurat_obj <- readRDS("seurat_obj.RDS")

Idents(seurat_obj) <- "Cell_Day" #Example

# Differential markers 

neuron_Day0 <- FindMarkers(seurat_obj, ident.1 = "Neuron_D0", min.pct = 0.2, only.pos = T)
neuron_Day7 <- FindMarkers(seurat_obj, ident.1 = "Neuron_D7", min.pct = 0.2, only.pos = T)
neuron_Day21 <- FindMarkers(seurat_obj, ident.1 = "Neuron_D21", min.pct = 0.2, only.pos = T)

# Top differential markers

neuron_Day0 <- neuron_Day0 %>% top_n(200, avg_log2FC)
neuron_Day7 <- neuron_Day7 %>% top_n(200, avg_log2FC)
neuron_Day21 <- neuron_Day21 %>% top_n(200, avg_log2FC)

# Convert gene symbol to Enterez Id

neuron_Day0 <- neuron_Day0 %>% rownames_to_column(var = "gene")
x = bitr(neuron_Day0$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
neuron_Day0 <- merge(neuron_Day0, x, by.x = "gene", by.y = "SYMBOL", all = FALSE)


neuron_Day7 <- neuron_Day7 %>% rownames_to_column(var = "gene")
x= bitr(neuron_Day7$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
neuron_Day7 <- merge(neuron_Day7, x, by.x = "gene", by.y = "SYMBOL", all = FALSE)


neuron_Day21 <- neuron_Day21 %>% rownames_to_column(var = "gene")
x= bitr(neuron_Day21$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
neuron_Day21 <- merge(neuron_Day21, x, by.x = "gene", by.y = "SYMBOL", all = FALSE)

# Make a list of vectors that need to be comapred
comp <-
  list(Day0 = c(neuron_Day0$ENTREZID),
       Day7 = c(neuron_Day7$ENTREZID),
       Day21 = c(neuron_Day21$ENTREZID))


####----------------  Compare GO Pathways ----------------####

go_neuron_comp <-  compareCluster(geneCluster = comp, 
                                  fun = "enrichGO",
                                  keyType="ENTREZID", 
                                  OrgDb = org.Mm.eg.db) #organism is mouse

go_neuron_comp@compareClusterResult <- go_neuron_comp[go_neuron_comp@compareClusterResult$ID %in% 
                                                        names(go_neuron_comp@compareClusterResult$ID))[table(go_neuron_comp@compareClusterResult$ID) == 1], ]

# Plot GO Pathways
dotplot(go_neuron_comp,showCategory = 15) + ggtitle("GO")


go_neuron_comp@compareClusterResult$ID




####----------------  Compare KEGG Pathways ----------------####

kegg_neuron_comp <- compareCluster(geneCluster = comp,
                                   fun = "enrichKEGG", 
                                   organism="mmu")

kegg_neuron_comp@compareClusterResult

kegg_neuron_comp@compareClusterResult$Description <- sub(" - .*", "", kegg_neuron_comp@compareClusterResult$Description)

# Define colors for specific labels
label_colors <- c("MAPK signaling pathway" = "red", "cAMP signaling pathway" = "red", "Adrenergic signaling in cardiomyocytes" = "red")

dotplot(kegg_neuron_comp, showCategory = 20) + ggtitle("Differential pathways in Nerve (DRG)")  + theme(
  axis.text.y = element_text(color = ifelse(levels(factor(kegg_neuron_comp@compareClusterResult$Description)) %in% names(label_colors),
                                            label_colors[levels(factor(data$category))],
                                            "black")))


####----------------  Compare WIKI Pathways ----------------####

wiki_neuron_comp <- compareCluster(geneCluster =comp,
                                    fun = "enrichWP",
                                   organism="Mus musculus", 
                                   pvalueCutoff=0.05)

dotplot(wiki_neuron_comp) + ggtitle("Wiki Pathways")


####----------------  Compare ENRICH Pathways ----------------####

enrich_neuron_comp <- compareCluster(geneCluster =comp,
                                    fun = "enrichPathway", 
                                    pvalueCutoff=0.05)

dotplot(enrich_neuron_comp) + ggtitle("path")


