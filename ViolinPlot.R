library(Seurat)
library(scCustomize)

# Load Seurat object
seurat_obj <- readRDS("seurat_obj.RDS")

# Single feature plot
VlnPlot(seurat_obj, features = "CD79A")

# Multiple features with 2 columns
VlnPlot(seurat_obj, features = c("CD79A", "MS4A1"), ncol = 2)

# Customized plot with statistical comparisons
VlnPlot(seurat_obj, features = "CD8A",
        pt.size = 0.1,          # Point size (0 to hide)
        cols = c("red", "blue"), # Custom colors
        split.by = "Group",     # Split by metadata column
        log = TRUE,             # Log-transform y-axis
        raster = TRUE           # Rasterize points for large datasets
) + 
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Group1", "Group2")),
    label = "p.signif"
  )

################################################################################

# Horizontal stacking
VlnPlot(seurat_obj, 
        features = c("CD79A", "MS4A1", "CD8A"), 
        stack = TRUE,
        flip = TRUE  # Switch axes
)


################################################################################



# Install if needed: remotes::install_github("samuel-marsh/scCustomize")

# Stacked plot with metadata variables
Stacked_VlnPlot(
  seurat_object = seurat_obj,
  features = c("CD79A", "percent_mito", "nFeature_RNA"),
  x_lab_rotate = TRUE,
  plot_spacing = 0.3,    # Space between plots
  colors_use = c("red", "blue", "green")
)

################################################################################

# Combine with DimPlot
library(patchwork)

p1 <- VlnPlot(seurat_obj, "CD79A", pt.size = 0)
p2 <- DimPlot(seurat_obj, reduction = "umap")

p1 + p2 + plot_layout(widths = c(1, 2))
