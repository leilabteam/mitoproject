
########### Load R packages
getwd()
setwd("E:/Rpractice/mitocodrial/plotting")
getwd()
install.packages("BiocManager")
BiocManager::install("ggplot2", force = TRUE)
BiocManager::install("dplyr")
BiocManager::install("Seurat")
BiocManager::install("cowplot")
BiocManager::install("harmony")
BiocManager::install("tidyverse")

na
library("ggplot2")
library("dplyr")
library("Seurat")
library("cowplot")
library("harmony")
library("tidyverse")
dplyr()
x1 = ggplot()
install.packages("ggplot2")
### Custom R functions
```{r}
function_dir <- "E:/Rpractice/mitocodrial/plotting"
getwd()
setwd("E:/Rpractice/mitocodrial/plotting")
source(paste(function_dir, "/Enrichment_GO_KEGG_GSEA_etc.R", sep=""))
source(paste(function_dir, "\\fundmental_function.R", sep=""))
source(paste(function_dir, "\\Planarian_special_function.R", sep=""))
source(paste(function_dir, "\\t-test_deseq2_annova_etc.R", sep=""))
### Set common variables
tmp.pre.data <- format(Sys.Date(), "%Y%m%d-")
project.name <- "mito.high.low-harmony"
## Define fixed gene set
```{r}
# Genes to remove: rRNA (first four, uncertain), and potentially incorrect RNA (mapped to mitochondria)
rRNA.genes_and_wrong.RNA <- c(
  "SMED30032887", "SMED30012309", "SMED30032663", "SMED30027845",
  "SMED30001886", "SMED30002990", "SMED30004246", "SMED30004656",
  "SMED30005051", "SMED30033740", "SMED30014993", "SMED30025308", 
  "SMED30004214", "SMED30005334", "SMED30012906"
)

# 12 mitochondrial-encoded proteins:
mito.genes <- c(
  "SMED30000702", "SMED30001799", "SMED30003686", "SMED30008414", 
  "SMED30010253", "SMED30010375", "SMED30012212", "SMED30013571", 
  "SMED30019959", "SMED30026531", "SMED30031308", "SMED30034623"
)
getwd()
# Import data:
high.data <- Read10X(data.dir = "./high_filtered_feature_bc_matrix/")
dim(high.data)
low.data <- Read10X(data.dir = "./low_filtered_feature_bc_matrix/")
dim(low.data)

# Label cell names
colnames(high.data) <- paste("high", colnames(high.data), sep = "")
colnames(low.data) <- paste("low", colnames(low.data), sep = "")
tmp.all.cell <- c(paste("", colnames(high.data), sep = ""), paste("", colnames(low.data), sep = ""))
names(tmp.all.cell) <- c(rep("C-high", ncol(high.data)), rep("A-low", ncol(low.data)))

# Merge multiple samples
tmp.count <- cbind(high.data, low.data)
view("tmp.count")
print("tmp.count")
view(tmp.count)
dim(tmp.count)
#tmp.count <- high.data
# Remove rRNA
tmp.count <- tmp.count[-match(c("SMED30032663", "SMED30027845", "SMED30032887", "SMED30015328", "SMED30012309", "SMED30027385", "SMED30023781"), rownames(tmp.count)), ]
dim(tmp.count)

# Create Seurat object and perform quality control:
data.combined <- CreateSeuratObject(counts = tmp.count, project = "mito", min.cells = 3, min.features = 200)
dim(data.combined)
# Input mitochondrial proportion:
data.combined[["percent.mt"]] <- PercentageFeatureSet(data.combined, features = c("SMED30031308", "SMED30034623", "SMED30019959", "SMED30010375", "SMED30012212", "SMED30008414", "SMED30026531", "SMED30001799", "SMED30005014", "SMED30013300", "SMED30017284", "SMED30024106", "SMED30012974", "SMED30022128", "SMED30003672", "SMED30024888", "SMED30024665", "SMED30000702"))
data.combined <- subset(data.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
dim(data.combined)
data.combined <- subset(data.combined, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15)
data.combined <- subset(data.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

dim(data.combined)
view(data.combined)
# Input stim metadata
data.combined@meta.data$stim <- names(tmp.all.cell[na.omit(match(as.character(colnames(data.combined)), tmp.all.cell))])
# Normalization
data.combined <- NormalizeData(data.combined, normalization.method = "LogNormalize", scale.factor = 10000)
data.combined <- FindVariableFeatures(data.combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data.combined)
data.combined <- ScaleData(data.combined, features = all.genes, verbose = FALSE) 
gc()
data.combined <- RunPCA(data.combined, features = VariableFeatures(object = data.combined), verbose = FALSE)

# Harmony analysis
options(repr.plot.height = 2.5, repr.plot.width = 6)
data.combined <- data.combined %>% 
  RunHarmony("stim", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(data.combined, 'harmony')
harmony_embeddings[1:5, 1:5]

pca_num <- 10
set.seed(42)

view(data.combined)

# Downstream analysis
data.combined <- RunUMAP(data.combined, reduction = "harmony", dims = 1:10) 
data.combined <- FindNeighbors(data.combined, reduction = "harmony", dims = 1:pca_num) 
data.combined <- FindClusters(data.combined, resolution = 0.6) 
unique(Idents(data.combined))
getwd()

plot.width <- 1200
plot.height <- 1000

reduction <- "umap"
#data.combined <- RunUMAP(data.combined, dims = 1:pca_num)
data.combined <- RunUMAP(data.combined, reduction = "harmony", dims = 1:10)
p1 <- DimPlot(data.combined, reduction = reduction, label = F)  
p2 <- DimPlot(data.combined, reduction = reduction, label = TRUE)
p3combined <- DimPlot(data.combined, reduction = reduction, split.by = "orig.ident", label = T)
p5 <- DimPlot(data.combined, reduction = reduction, split.by = "orig.ident", label = F)
p4 <- DimPlot(data.combined, reduction = reduction, group.by = "orig.ident", label = F, order = F)
p6 <- DimPlot(data.combined, reduction = reduction, group.by = "orig.ident", label = T)
p7 = ggplot(data.combined, reduction = reduction, label = F) 
print(p3combined)
Idents(data.combined) = 'seurat_clusters'
p = DimPlot(data.combined, reduction = "umap", pt.size = 0.6, split.by = 'stim', label = TRUE)
pmerge = DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.4, label = TRUE)
print(pmerge)
pm6enuleadmerge = DimPlot(enucledm6amerge, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.4, label = TRUE)
print(pm6enuleadmerge)
pdf('umap.pdf', width = 14)
print(p)
print(p5)
print(p4)

####### FindMarkers Single-cell Annotation
cluster.markers <- FindMarkers(object = PRO, ident.1 = 3, ident.2 = 4, min.pct = 0.1, logfc.threshold = 0.25)
write.table(data.frame(gene_id = rownames(cluster.markers), cluster.markers), file = paste('out_cluster', '3VS4', '_diffgenes.xls', sep = ''), sep = '\t', quote = F, row.names = F)

# Plot output:

p1 <- DimPlot(data.combined, reduction = reduction, label = F)  
umap.name <- paste(reduction, "-", tmp.pre.data, project.name, plot.other, "no.label.png", sep = "")
png(umap.name, width = plot.width, height = plot.height, res = 72*3)
print(p1)
dev.off()
umap.name <- paste(reduction, "-", tmp.pre.data, project.name, plot.other, "label.png", sep = "")

png(umap.name, width = plot.width, height = plot.height, res = 72*3)
print(p2)
dev.off()
umap.name <- paste(reduction, "-", tmp.pre.data, project.name, plot.other, "split.stim.png", sep = "")
png(umap.name, width = plot.width*2.8, height = plot.height, res = 72*3)
print(p3)
print(p3combined)
dev.off()
umap.name <- paste(reduction, "-", tmp.pre.data, project.name, plot.other, "split.stim.no.label.png", sep = "")
png(umap.name, width = plot.width*2.8, height = plot.height, res = 72*3)
print(p5)
dev.off()
umap.name <- paste(reduction, "-", tmp.pre.data, project.name, plot.other, "group.stim.nolabel.png", sep = "")
png(umap.name, width = plot.width*1.2, height = plot.height, res = 72*3)
print(p4)
dev.off()
umap.name <- paste(reduction, "-", tmp.pre.data, project.name, plot.other, "group.stim.label.png", sep = "")
png(umap.name, width = plot.width*1.2, height = plot.height, res = 72*3)
print(p6)
dev.off()
umap.name <- paste(reduction, "-", tmp.pre.data, project.name, plot.other, "label-group.stim-combine.png", sep = "")
png(umap.name, width = plot.width*2, height = plot.height, res = 72*3)
print(plot_grid(p4, p2))
dev.off()

# Define the list of genes
genelist <- c('markergenes.csv')  # Replace with actual gene names

# Loop through each gene in the list
for (i in 1:length(genelist)){
  
  # Plot the feature for the current gene using UMAP reduction
  FeaturePlot(data.combined, features = genelist[i], reduction = 'umap', pt.size = 0.1, split.by = "orig.ident")
  
  # Save the plot to a PDF file named after the gene
  ggsave(paste(genelist[i], 'umap.pdf', sep = '_'), width = 10, height = 4.6)
}

######## cell cycle analysis
## Load required packages
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(multtest)

# Define cell cycle genes for G2M, G1, and S phases
G2Mgenecellcycle = c('SMED30009506','SMED30015678','SMED30006888','SMED30003246','SMED30014970','SMED30025245','SMED30004904','SMED30033670','SMED30015929','SMED30018221','SMED30004988','SMED30030198','SMED30033981','SMED30000949','SMED30025778','SMED30010251','SMED30022468','SMED30022580','SMED30011277','SMED30023082','SMED30019915','SMED30000277','SMED30033646','SMED30023953','SMED30016942','SMED30032426','SMED30030055','SMED30019616','SMED30027106','SMED30011934','SMED30033017','SMED30030505','SMED30006092','SMED30016301','SMED30011977','SMED30011091','SMED30024586','SMED30009683','SMED30010966','SMED30006356','SMED30019390','SMED30035314','SMED30012613','SMED30005609','SMED30024902','SMED30021702','SMED30033763','SMED30005185','SMED30004873','SMED30026475','SMED30010421','SMED30012777','SMED30019744','SMED30028681','SMED30023769','SMED30022293','SMED30026807','SMED30020865','SMED30006304','SMED30035063','SMED30029272','SMED30021639') %>% unique()

g1genecellcycle = c('SMED30032563','SMED30001802','SMED30013808','SMED30034500','SMED30002663','SMED30010350','SMED30035785','SMED30028898','SMED30011050','SMED30033019','SMED30035929','SMED30025384','SMED30027428','SMED30016879','SMED30010205','SMED30024929','SMED30007715','SMED30000968','SMED30001916','SMED30028883','SMED30021545','SMED30033482','SMED30022027','SMED30025635','SMED30024446','SMED30023388','SMED30025126','SMED30020940','SMED30000740','SMED30035118','SMED30030413','SMED30031519','SMED30026593','SMED30006028','SMED30021490','SMED30014999','SMED30005488','SMED30027370','SMED30022293','SMED30026807','SMED30020865','SMED30006304','SMED30009429','SMED30034254','SMED30027052','SMED30005066','SMED30004207','SMED30004214','SMED30015762','SMED30016919','SMED30012240','SMED30008690','SMED30020865','SMED30015642','SMED30023769','SMED30022293','SMED30026807','SMED30020865','SMED30006304','SMED30009429','SMED30034254','SMED30027052','SMED30005066','SMED30004207','SMED30004214','SMED30015762','SMED30016919','SMED30012240','SMED30008690','SMED30020865','SMED30015642','SMED30023769','SMED30033566','SMED30023769','SMED30008884','SMED30035485','SMED30033583','SMED30020397','SMED30035424','SMED30035357','SMED30011055','SMED30012984','SMED30008787','SMED30018508','SMED30008579','SMED30027813','SMED30026358','SMED30008911','SMED30016919','SMED30015120','SMED30003738','SMED30018762','SMED30019873','SMED30005982','SMED30022344','SMED30006268') %>% unique()

sgenecellcycle  = c('SMED30032563','SMED30022389','SMED30004470','SMED30000396','SMED30023821','SMED30029973','SMED30017135','SMED30020890','SMED30018837','SMED30010424','SMED30022437','SMED30021509','SMED30026247','SMED30010763','SMED30031773','SMED30034520','SMED30010517','SMED30033683','SMED30002348','SMED30025421','SMED30019555','SMED30017245','SMED30000541','SMED30006456','SMED30005215','SMED30003162') %>% unique()

# Set default assay and normalize data
DefaultAssay(data.combined) <- "RNA"
data.combined <- NormalizeData(data.combined)

# Cell cycle scoring
data.combined <- CellCycleScoring(data.combined, s.features = sgenecellcycle,
                                  g2m.features = G2Mgenecellcycle, set.ident = TRUE)

# UMAP plots by 'orig.ident'
DimPlot(data.combined, split.by = "orig.ident")

# Subset for G2M phase and plot UMAP
G2M <- subset(data.combined, subset = Phase == "G2M")
DimPlot(G2M, group.by = "seurat_clusters", split.by = "orig.ident")

# Subset for S phase and plot UMAP
S <- subset(data.combined, subset = Phase == "S")
DimPlot(S, group.by = "seurat_clusters", split.by = "orig.ident")

# Subset for G1 phase and plot UMAP
G1 <- subset(data.combined, subset = Phase == "G1")
DimPlot(G1, group.by = "seurat_clusters", split.by = "orig.ident")

##### Count and Calculate Proportions of Cell Types #####
# Count number of cells per cell type and sample
cell.number <- table(data.combined@active.ident, data.combined$orig.ident) %>% as.data.frame.array()
write.csv(cell.number,  file = "cell.number.csv")

# Convert to cell type proportions
cell.ratio <- prop.table(table(data.combined@active.ident, data.combined$orig.ident), margin = 2) %>%
  as.data.frame.array()
write.csv(cell.ratio,  file = "cell.ratio.csv")

##### Analysis for G2M Phase #####
# Subset for G2M phase
x = subset(data.combined, idents = "G2M")

# Count number of cells per cluster and sample
cell.number <- table(x$seurat_clusters, x$orig.ident) %>% as.data.frame.array()
write.csv(cell.number,  file = "G2M.cell.number.csv")

# Convert to cell type proportions
cell.ratio <- prop.table(table(x$seurat_clusters, x$orig.ident), margin = 2) %>%
  as.data.frame.array()
write.csv(cell.ratio,  file = "G2M.cell.ratio.csv")

# Save the modified Seurat object
save(data.combined, file = "data.combined.RData")


####### Subcluster analysis
#!/usr/bin/env Rscript
library(Seurat)

setwd('E:/Rpractice/mitocodrial/2024mito/subcluster/topmarker/muscle')

# Display the distribution of 'stim' values
table(data.mitohigh_low@meta.data$stim)

# Subset cells for C-high and A-low
data.mitohigh_low <- subset(datacombined, subset = stim %in% c("C-high", "A-low"))

# Subset cells for A-low
data.mitolow <- subset(datacombined, subset = stim == "A-low")

# Subset cells for clusters 2 and 3 (highlow23)
data.mitohigh_low23 <- subset(data.mitohigh_low, subset = seurat_clusters %in% c(2, 3))
head(data.mitohigh_low23@meta.data$seurat_clusters)
dim(data.mitohigh_low23)

# Subset cells for clusters 1, 5, 11 (highlowepider)
data.mitohigh_lowepider <- subset(data.mitohigh_low, subset = seurat_clusters %in% c(1, 5, 11))
head(data.mitohigh_lowepider@meta.data$seurat_clusters)
dim(data.mitohigh_lowepider)

# Subset cells for muscle clusters 0, 7, 8, 9, 10, 19
data.mitohigh_lowmuscle <- subset(data.mitohigh_low, subset = seurat_clusters %in% c(0, 7, 8, 9, 10, 19))
head(data.mitohigh_lowmuscle@meta.data$seurat_clusters)
dim(data.mitohigh_lowmuscle)

# Subset cells for low clusters 2 and 3
data.mitolow23 <- subset(data.mitolow, subset = seurat_clusters %in% c(2, 3))
head(data.mitolow23@meta.data$seurat_clusters)
dim(data.mitolow23)

# Subset cells for clusters 4, 15, 17, 20 (highlowcathepsin)
data.mitohigh_lowcathepsin <- subset(data.mitohigh_low, subset = seurat_clusters %in% c(4, 15, 17, 20))
head(data.mitohigh_lowcathepsin@meta.data$seurat_clusters)
dim(data.mitohigh_lowcathepsin)

# Subset cells for neuron clusters 18, 21 (highlowneuron)
data.mitohigh_lowneuron <- subset(data.mitohigh_low, subset = seurat_clusters %in% c(18, 21))
head(data.mitohigh_lowneuron@meta.data$seurat_clusters)
dim(data.mitohigh_lowneuron)

# Subset cells for intestine clusters 6, 13 (highlowintestine)
data.mitohigh_lowintestine <- subset(data.mitohigh_low, subset = seurat_clusters %in% c(6, 13))
head(data.mitohigh_lowintestine@meta.data$seurat_clusters)
dim(data.mitohigh_lowintestine)

# Subset cells for parenchymal clusters 12, 16, 22 (highlowparenchymal)
data.mitohigh_lowparenchymal <- subset(data.mitohigh_low, subset = seurat_clusters %in% c(12, 16, 22))
head(data.mitohigh_lowparenchymal@meta.data$seurat_clusters)
dim(data.mitohigh_lowparenchymal)

# Subset cells for protonephridia cluster 14 (highlowprotonephridia)
data.mitohigh_lowprotonephridia <- subset(data.mitohigh_low, subset = seurat_clusters == 14)
head(data.mitohigh_lowprotonephridia@meta.data$seurat_clusters)
dim(data.mitohigh_lowprotonephridia)

# Processing for clusters 2 and 3 (highlow23)
data.mitohigh_low23 <- NormalizeData(data.mitohigh_low23) %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA()

ElbowPlot(data.mitohigh_low23, ndims = 30)

pc.num = 1:20
data.mitohigh_low23 <- FindNeighbors(data.mitohigh_low23, dims = pc.num) %>%
        FindClusters(resolution = 0.8) %>%
        RunUMAP(dims = pc.num)

# Processing for clusters 1, 5, 11 (highlowepider)
data.mitohigh_lowepider <- NormalizeData(data.mitohigh_lowepider) %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA()

ElbowPlot(data.mitohigh_lowepider, ndims = 30)

pc.num = 1:20
data.mitohigh_lowepider <- FindNeighbors(data.mitohigh_lowepider, dims = pc.num) %>%
        FindClusters(resolution = 0.8) %>%
        RunUMAP(dims = pc.num)

# Processing for muscle clusters 0, 7, 8, 9, 10, 19 (muscle)
data.mitohigh_lowmuscle <- NormalizeData(data.mitohigh_lowmuscle) %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA()

ElbowPlot(data.mitohigh_lowmuscle, ndims = 30)

pc.num = 1:20
data.mitohigh_lowmuscle <- FindNeighbors(data.mitohigh_lowmuscle, dims = pc.num) %>%
        FindClusters(resolution = 0.8) %>%
        RunUMAP(dims = pc.num)

# Processing for low clusters 2 and 3 (low23)
data.mitolow23 <- NormalizeData(data.mitolow23) %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA()

ElbowPlot(data.mitolow23, ndims = 30)

pc.num = 1:20
data.mitolow23 <- FindNeighbors(data.mitolow23, dims = pc.num) %>%
        FindClusters(resolution = 0.8) %>%
        RunUMAP(dims = pc.num)

# Build cluster tree for muscle clusters
data.mitohigh_lowmuscle <- FindClusters(data.mitohigh_lowmuscle, resolution = seq(0.1, 1, 0.1))
data.mitohigh_lowmuscle <- BuildClusterTree(data.mitohigh_lowmuscle)

##### Analysis for Cluster 2 #####
# Subset for cluster 2
c2 = subset(data.combined, subset = seurat_clusters == "2")

# Count number of cells per phase and sample
cell.number <- table(c2$Phase, c2$orig.ident) %>% as.data.frame.array()
write.csv(cell.number,  file = "cluster2.cell.number.csv")

# Convert to cell type proportions
cell.ratio <- prop.table(table(c2$Phase, c2$orig.ident), margin = 2) %>%
  as.data.frame.array()
write.csv(cell.ratio,  file = "cluster2.cell.ratio.csv")

# Save the modified Seurat object
save(data.combined, file = "data.combined.RData")


######## average heatmap ######## Average Expression Heatmap

library(ggtext)
library(grid)
color_vector <- colorRampPalette(c('#0099FF', '#FFFFCC', '#FF3300'))(100)                                              
pdf('averageg2s.pdf')

source('function.R')

AverageHeatmap(object = data.combined, htCol = c('#FFFFCC', '#FF3300'), htRange = c(0, 1.5), border = F,
               markerGene = g1genecellcycle, fontsize = 6, width = 9, annoCol = TRUE,
               myanCol = c("#F8766D", "#00BFC4"))

######## Scale data ########

sub_34 = ScaleData(sub_34, features = g1genecellcycle)
sub_34 = ScaleData(sub_34, features = G2Mgenecellcycle)

max_cells_per_cluster <- 500

##### Loop through each cluster, selecting a maximum of 500 cells ######
clusters <- unique(Idents(sub_34))
subset_cells <- c()
for (cluster in clusters) {
  cluster_cells <- WhichCells(sub_34, idents = cluster)
  if (length(cluster_cells) > max_cells_per_cluster) {
    cluster_cells <- sample(cluster_cells, max_cells_per_cluster)
  }
  subset_cells <- c(subset_cells, cluster_cells)
}
sub_34 = subset(sub_34, cells = subset_cells)
levels(sub_34)
Idents(sub_34) = factor(Idents(sub_34), levels = c('other', 'fos34neg'))
Idents(sub_34) = 'type2'

#################### gene expression test #################### Gene Expression Detection #######
data.combined

rm(list = c("umap.name", "p1", "p2", "p3", "p4", "p5", "p6"))