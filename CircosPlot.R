# Load packages
library(Seurat)
library(tidyverse)
library(circlize)
library(RColorBrewer)

# Load Seurat object
seurat_obj <- readRDS("seurat_obj.RDS")

# Prepare data for circos plot
ggData <- data.frame(
  Group = seurat_obj$Group,
  CellType = seurat_obj$CellName
) %>% 
  count(Group, CellType) %>% 
  group_by(Group) %>% 
  mutate(proportion = n/sum(n))

# Create matrix for circos plot (group-celltype relationships)
circos_matrix <- ggData %>% 
  select(-n) %>% 
  pivot_wider(names_from = Group, values_from = proportion, values_fill = 0) %>% 
  column_to_rownames("CellType") %>% 
  as.matrix()

# Set up color scheme
group_colors <- brewer.pal(n = length(unique(seurat_obj$Group)), name = "Set3")
celltype_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(nrow(circos_matrix))

# Initialize circos plot
pdf("circos_plot.pdf", width = 8, height = 8)
par(cex = 0.8, mar = c(0, 0, 0, 0))

# Create chord diagram
chordDiagram(
  x = circos_matrix,
  grid.col = c(setNames(group_colors, colnames(circos_matrix)),
               setNames(celltype_colors, rownames(circos_matrix))),
  transparency = 0.25,
  directional = 1,
  direction.type = "arrows",
  link.arr.type = "big.arrow",
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.2)
)

# Add labels
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(
      mean(xlim), 
      ylim[1], 
      sector.index, 
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.8
    )
  },
  bg.border = NA
)

# Add legend
legend(
  x = "bottomleft", 
  legend = colnames(circos_matrix),
  fill = group_colors,
  border = NA,
  bty = "n",
  title = "Groups"
)

dev.off()
