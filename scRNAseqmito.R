
###########加载R包
getwd()
setwd("E:/Rpractice/mitocodrial/plotting")
getwd()
install.packages("BiocManager")
BiocManager::install("ggplot2",force = TRUE)
BiocManager::install("dplyr")
BiocManager::install("Seurat")
BiocManager::install("cowplot")
BiocManager::install("harmony")
BiocManager::install("tidyverse")
save.image("2022811")

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
###自定义R 函数
```{r}
function_dir <- "E:/Rpractice/mitocodrial/plotting"
getwd()
setwd("E:/Rpractice/mitocodrial/plotting")
source( paste(function_dir,"/Enrichment_GO_KEGG_GSEA_etc.R",sep="") )
source( paste(function_dir,"\\fundmental_function.R",sep="") )
source( paste(function_dir,"\\Planarian_special_function.R",sep="") )
source( paste(function_dir,"\\t-test_deseq2_annova_etc.R",sep="") )
### 常用变量设定
tmp.pre.data<- format(Sys.Date(), "%Y%m%d-")
project.name<- "mito.high.mid.low-harmony"
## 确定固定基因集
```{r}
#需要删掉的rRNA(前四个，不确定都是)，以及可能错误的RNA（比对到线粒体上）
rRNA.genes_and_wrong.RNA<-  c(
  "SMED30032887","SMED30012309","SMED30032663","SMED30027845",
  "SMED30001886", "SMED30002990", "SMED30004246", "SMED30004656",
  "SMED30005051", "SMED30033740", "SMED30014993", "SMED30025308", 
  "SMED30004214", "SMED30005334", "SMED30012906")

#12个线粒体编码蛋白：
mito.genes <-c(
  "SMED30000702", "SMED30001799", "SMED30003686", "SMED30008414", 
  "SMED30010253", "SMED30010375", "SMED30012212", "SMED30013571", 
  "SMED30019959", "SMED30026531", "SMED30031308", "SMED30034623")
getwd()
#导入数据：
high.data<- Read10X(data.dir = "./high_filtered_feature_bc_matrix/")
dim(high.data)

mid.data<- Read10X(data.dir = "./mid_filtered_feature_bc_matrix/")
dim(mid.data)

low.data<- Read10X(data.dir = "./low_filtered_feature_bc_matrix/")
dim(low.data)
save.image("20228011")
#对cell名打标签
colnames(high.data)<- paste("high",colnames(high.data),sep="")
colnames(mid.data)<- paste("mid",colnames(mid.data),sep="")
colnames(low.data)<- paste("low",colnames(low.data),sep="")
tmp.all.cell<- c(paste("",colnames(high.data),sep=""),paste("",colnames(mid.data),sep=""),paste("",colnames(low.data),sep=""))
names(tmp.all.cell)<- c(rep("C-high",ncol(high.data)),rep("B-mid",ncol(mid.data)),rep("A-low",ncol(low.data)))

#合并多个样本
tmp.count<- cbind(high.data, mid.data,low.data)
view ("tmp.count")
print("tmp.count")
view(tmp.count)
dim(tmp.count)
#tmp.count<- high.data
#去除rRNA
tmp.count<-tmp.count[-match(c("SMED30032663" ,"SMED30027845","SMED30032887","SMED30015328" ,"SMED30012309","SMED30027385","SMED30023781"),rownames(tmp.count)),]
dim(tmp.count)

#创建seurat对象并质控：
data.combined <- CreateSeuratObject(counts = tmp.count, project = "mito", min.cells = 3, min.features = 200 )
dim(data.combined)
#输入线粒体比例：
data.combined[["percent.mt"]] <- PercentageFeatureSet(data.combined, features = c("SMED30031308", "SMED30034623" ,"SMED30019959", "SMED30010375" ,"SMED30012212" ,"SMED30008414" ,"SMED30026531","SMED30001799" ,"SMED30005014" ,"SMED30013300","SMED30017284", "SMED30024106", "SMED30012974","SMED30022128","SMED30003672" ,"SMED30024888","SMED30024665","SMED30000702"))
data.combined<- subset(data.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 )
dim(data.combined)
data.combined<- subset(data.combined, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15 )
data.combined<- subset(data.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 )

dim(data.combined)
view(data.combined)
#输入stim 元变量
data.combined@meta.data$stim<- names(tmp.all.cell[na.omit(match(as.character(colnames(data.combined)),tmp.all.cell))])
#标椎化处理
data.combined <- NormalizeData(data.combined, normalization.method = "LogNormalize", scale.factor = 10000)
data.combined <- FindVariableFeatures(data.combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data.combined)
data.combined <- ScaleData(data.combined, features = all.genes,verbose = FALSE) 
gc()
data.combined <- RunPCA(data.combined, features = VariableFeatures(object = data.combined), verbose = FALSE)


#harmony分析
options(repr.plot.height = 2.5, repr.plot.width = 6)
data.combined <- data.combined %>% 
  RunHarmony("stim", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(data.combined, 'harmony')
harmony_embeddings[1:5, 1:5]


pca_num<-10
set.seed(42)

view(data.combined)

# downstrean analysis
data.combined <-RunUMAP(data.combined,reduction = "harmony", dims = 1:10) 
data.combined <-FindNeighbors(data.combined,reduction = "harmony", dims = 1:pca_num) 
data.combined <-FindClusters(data.combined,resolution = 0.6) 
unique(Idents(data.combined))
data.combined <-RunTSNE(data.combined,reduction = "harmony", dims = 1:10)
getwd()
# 导入seurat3.1.1 分析的结果：
```{r}
load("./fengcun-20200821.harmony.rdata")  # 此结果来自于seurat 3.1.1
data.combined@meta.data$orig.ident <- data.combined@meta.data$stim
# umap绘图
```{r }
graphics.off()
#Uamp和tsne绘图，如果都要绘制，则重复两次
#根据实际参数调整；
plot.other<-""
#绘图内容,二选一
reduction<- "tsne"

plot.width <- 1200
plot.height <- 1000

reduction<- "umap"
#data.combined <- RunUMAP(data.combined,  dims = 1:pca_num)
data.combined <-RunTSNE(data.combined,reduction = "harmony", dims = 1:10)
p1 <- DimPlot(data.combined, reduction = reduction, label = F)  
p2 <- DimPlot(data.combined, reduction = reduction, label = TRUE)
p3combined<-DimPlot(data.combined, reduction = reduction, split.by = "orig.ident",label = T)
p5<-DimPlot(data.combined, reduction = reduction, split.by = "orig.ident",label = F)
p4<-DimPlot(data.combined, reduction = reduction, group.by = "orig.ident",label = F,order=F)
p6<-DimPlot(data.combined, reduction = reduction, group.by = "orig.ident",label = T)
p7 = ggplot(data.combined, reduction = reduction, label = F) 
print(p3combined)
Idents(data.combined) = 'seurat_clusters'
p=DimPlot(data.combined, reduction = "umap",pt.size=0.6,split.by='stim',label = TRUE)
pmerge=DimPlot(data.combined, reduction = "umap",group.by = "seurat_clusters",pt.size=0.4,label = TRUE)
print(pmerge)
pm6enuleadmerge = DimPlot(enucledm6amerge, reduction = "umap",group.by = "seurat_clusters",pt.size=0.4,label = TRUE)
print(pm6enuleadmerge)
pdf('umap.pdf',width=14)
print(p)
print(p5)
print(p4)


#######单细胞细胞注释
cluster.markers <- FindMarkers(object=PRO,ident.1=3,ident.2=4,min.pct=0.1,logfc.threshold = 0.25)
write.table(data.frame(gene_id=rownames(cluster.markers),cluster.markers),file=paste('out_cluster','3VS4','_diffgenes.xls',sep=''),sep='\t',quote=F,row.names=F)




#绘图输出：
data.combined <-RunTSNE(data.combined,reduction = "harmony", dims = 1:10)
p1 <- DimPlot(data.combined, reduction = reduction, label = F)  
umap.name<- paste(reduction,"-",tmp.pre.data,project.name,plot.other,"no.label.png",sep="")
png(umap.name,width = plot.width,height = plot.height, res=72*3)
print(p1)
dev.off()
umap.name<- paste(reduction,"-",tmp.pre.data,project.name,plot.other,"label.png",sep="")
png(umap.name,width = plot.width,height = plot.height, res=72*3)
print(p2)
dev.off()
umap.name<- paste(reduction,"-",tmp.pre.data,project.name,plot.other,"split.stim.png",sep="")
png(umap.name,width = plot.width*2.8,height = plot.height, res=72*3)
print(p3)
print(p3combined)
dev.off()
umap.name<- paste(reduction,"-",tmp.pre.data,project.name,plot.other,"split.stim.no.label.png",sep="")
png(umap.name,width = plot.width*2.8,height = plot.height, res=72*3)
print(p5)
dev.off()
umap.name<- paste(reduction,"-",tmp.pre.data,project.name,plot.other,"group.stim.nolabel.png",sep="")
png(umap.name,width = plot.width*1.2,height = plot.height, res=72*3)
print(p4)
dev.off()
umap.name<- paste(reduction,"-",tmp.pre.data,project.name,plot.other,"group.stim.label.png",sep="")
png(umap.name,width = plot.width*1.2,height = plot.height, res=72*3)
print(p6)
dev.off()
umap.name<- paste(reduction,"-",tmp.pre.data,project.name,plot.other,"label-group.stim-combine.png",sep="")
png(umap.name,width = plot.width*2,height = plot.height, res=72*3)
print(plot_grid(p4, p2))
dev.off()


########循环出图#####
for (i in 1:length(data.combined)){
  
  FeaturePlot(L0L1,features = data.combined[i],reduction = 'umap',pt.size = 0.1,split.by = "orig.ident")
  ggsave(paste(data.combined[i],'umap.pdf',sep = '_'),width = 10,height = 4.6)
}

########平均表达值热图

library(ggtext)
library(grid)
color_vector <- colorRampPalette(c('#0099FF','#FFFFCC','#FF3300'))(100)                                              
pdf('averageg2s.pdf')


source('function.R')

AverageHeatmap(object = data.combined,htCol=c('#FFFFCC','#FF3300'),htRange=c(0,1.5),border=F,
               markerGene = g1genecellcycle,fontsize=6,width=9,annoCol = TRUE,
               myanCol = c("#F8766D","#00BFC4"))

########如果要scale data

sub_34 = ScaleData(sub_34,features=g1genecellcycle)
sub_34 = ScaleData(sub_34,features=G2Mgenecellcycle)


max_cells_per_cluster <- 500

#####循环每个cluster都是500个细胞
clusters <- unique(Idents(sub_34))
subset_cells <- c()
for (cluster in clusters) {
  cluster_cells <- WhichCells(sub_34, idents = cluster)
  if (length(cluster_cells) > max_cells_per_cluster) {
    cluster_cells <- sample(cluster_cells, max_cells_per_cluster)
  }
  subset_cells <- c(subset_cells, cluster_cells)
}
sub_34 = subset(sub_34,cells=subset_cells)
levels(sub_34)
Idents(sub_34) = factor(Idents(sub_34),levels=c('other','fos34neg'))
Idents(sub_34)= 'type2'









####################基因表达量检测
data.combined



rm(list=c("umap.name","p1","p2","p3","p4","p5","p6"))

getwd()
#检查nFeature,mtRNA分布是否合理
```{r}
#某些特征的分布绘图
reduction <- "tsne"
#reduction <- "umap"
if (reduction ==  "tsne"){
  plot.x <- "tSNE_1"
  plot.y <- "tSNE_2"
}else{
  plot.x <- "UMAP_1"
  plot.y <- "UMAP_2"
}
library(methods)
meta_matrix <- cbind(data.combined@reductions$tsne@cell.embeddings  ) 
meta_matrix<- cbind(meta_matrix,data.combined@meta.data)
library("ggplot2")
# 对cluster区分颜色进行绘图
#ggplot2中可以使用!!sym(plot.x) 来调用外部变量名，但是不适用于常规的R语言。
colnames(meta_matrix)

#mtRNA meta_matrix
g1 <- ggplot(meta_matrix, mapping = aes(x = !!sym(plot.x), y = !!sym(plot.y), colour = percent.mt ) )+
  geom_point(size = 0.4)+
  theme_light()+  scale_color_gradientn(colours = c("blue","white","red"))
png("g1.png")
pdf("g1.pdf")
print(g1.pdf)
print(g1)
ggsave(g1,filename = "meta_matrix")
dev.off()

