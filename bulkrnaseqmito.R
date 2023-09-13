## Introduction
This analysis of RNA-Seq was quality checked by FastQC 0.11.8/MultiQC 1.7 and 
trimmed by Trimmomatic 0.39 for HEADCROP 10.
After that, all samples were quantified by kallisto 0.46.0 with bias correction, 
100 times bootstrap and 42 seed.
```{r reading, cache=TRUE}
load("txi.rda")

## QC

How many genes were detected:
  
  ```{r qc,cache=TRUE,echo=FALSE}
apply(txi$counts,2,function(x) c(length(x[x > 0]),length(x[x >= 1]),
                                 length(x[x >= 5]),length(x[x >= 10]))) %>>% 
  rbind %>>% {rownames(.) <- c("> 0 gene numbers",">= 1 gene numbers",
                               ">= 5 gene numbers",">= 10 gene numbers");.}
```

Correlation heatmap

```{r correlation,cache=TRUE,echo=FALSE}
p <- cor(txi$abundance) %>>% 
  (pheatmap::pheatmap(.,color=colorRampPalette(c("grey88","orangered"))(100),
                      clustering_distance_rows="canberra",
                      clustering_distance_cols="canberra",
                      clustering_method="single",display_numbers=T))
print(p)
pdf("cor.pdf")
print(p)
dev.off()
```

Dimension reduction

```{r dr,cache=TRUE,echo=FALSE}
t(txi$abundance) %>>% dist %>>% cmdscale %>>%
  # edgeR::cpm(txi$counts) %>>% t %>>% dist %>>% cmdscale %>>% 
  {plot(.,xlab="PC1",ylab="PC2",pch=19,
        col=rep(c("orangered","cornflowerblue","darkgreen"),each=3));
    text(.[,1] + c(rep(1000,3),rep(0,3)),.[,2] + c(rep(0,3),rep(1000,3)),
         labels=colnames(txi$abundance))
  }
legend(12500,-5000,c("E","opa1"),pch=19,col=c("orangered","cornflowerblue"))
pdf("pca.pdf")
t(txi$abundance) %>>% dist %>>% cmdscale %>>%
  # edgeR::cpm(txi$counts) %>>% t %>>% dist %>>% cmdscale %>>% 
  {plot(.,xlab="PC1",ylab="PC2",pch=19,
        col=rep(c("orangered","cornflowerblue","darkgreen"),each=3));
    text(.[,1] + c(rep(1000,3),rep(0,3)),.[,2] + c(rep(0,3),rep(1000,3)),
         labels=colnames(txi$abundance))
  }
legend(12500,-5000,c("E","opa1"),pch=19,col=c("orangered","cornflowerblue"))
dev.off()
```

DEG (opa1 vs. E)

```{r deg1,cache=TRUE,echo=FALSE,message=FALSE}
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi,
                                data.frame(group=rep(c("E",
                                                       "opa1"),
                                                     each=3),
                                           row.names=colnames(txi$counts)),
                                design=~group)
des <- DESeq(dds)
res <- results(des,contrast=c("group","opa1","E"))
head(subset(res,padj < 0.05 & abs(log2FoldChange) >= 1))
nrow(subset(res,padj < 0.05 & abs(log2FoldChange) >= 1))
subset(res,padj < 0.05 & abs(log2FoldChange) >= 1) %>>% 
  write.csv(file="opa1.vs.E_deg.csv",quote=F)
```

Volcano plot (opa1 vs. E)

```{r volcano1,cache=TRUE,echo=FALSE}
d <- as.data.frame(res) %>>% {.[complete.cases(.),]}
d$Sig <- d$padj < 0.05
d$Sig[d$Sig == F] <- "NS"
d$Sig[d$log2FoldChange >= 1 & d$Sig == T] <- "Up regulated"
d$Sig[d$log2FoldChange <= -1 & d$Sig == T] <- "Down regulated"
d$Sig[d$Sig != "Up regulated" & d$Sig != "Down regulated"] <- "NS"
d$Sig <- relevel(as.factor(d$Sig),"NS")
d$Sig <- relevel(d$Sig,"Up regulated")
library(ggplot2)
library(ggrepel)
label_d <- rbind(head(d[order(d$log2FoldChange),],10),
                 d[which.max(d$log2FoldChange),],
                 d[which(d$log2FoldChange > 0),] %>>% 
                   (.[which.max(-log10(.$padj)),]))
p <- ggplot(d,aes(log2FoldChange,-log10(padj),color=Sig)) + 
  geom_point(shape=19,alpha=0.5) + 
  geom_text_repel(aes(log2FoldChange,-log10(padj),label=rownames(label_d)),
                  label_d,inherit.aes=F) +
  theme_light()
print(p)
pdf("vol.pdf")
print(p)
dev.off()
```

GO

```{r go,cache=TRUE,echo=FALSE,message=FALSE,warning=FALSE,fig.height=8}
deg <- read.csv("opa1.vs.E_deg.csv",row.names=1) %>>% rownames

library(clusterProfiler)
all_gene_go <-
  read.table("transcriptSmedwochong4DataLastGOAnotationIntergrate.txt",header=T)

ego <- enricher(deg,TERM2GENE=all_gene_go,minGSSize=3) %>>% as.data.frame

go_infor <- read.csv("GOInformation4Col.csv")

index <- match(ego[,1],go_infor[,1])
ego$ontology <- go_infor[index,2]
ego$Description <- as.character(go_infor[index,3])
ego$`function` <- go_infor[index,5]
ego <- ego[order(ego$p.adjust),]

write.csv(ego,"opa1.vs.E_deg_GO_enrich.csv",row.names=F,quote=F)

# library(dplyr)
library(data.table)
library(ggplot2)
# d <- ego[complete.cases(ego),] %>% arrange(p.adjust) %>%
#   group_by(ontology) %>% slice(1:10)
d <- as.data.table(ego)
d <- d[complete.cases(d),][order(p.adjust),head(.SD,10),keyby=ontology]
p <- ggplot(d,aes(forcats::fct_inorder(Description),Count,fill=ontology)) +
  geom_col(position="dodge",alpha=0.8) +
  theme_light() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        plot.margin=unit(c(0,0,0,2),"cm")) +
  labs(title="GO enrichment for opa1 vs. E DEGs",x="Description")
print(p)
pdf("go.pdf")
print(p)
dev.off()

g <- ggplot(d,aes(Count,forcats::fct_inorder(Description),size=p.adjust)) +
  geom_point(aes(color=ontology)) + 
  theme_light() +
  theme(axis.text.x=element_text(angle=90,hjust=0),
        axis.text.y=element_text(angle=-45,hjust=1,vjust=1,size=6),
        title=element_text(size=10,face="bold"),
        plot.margin=unit(c(1,0,0,2),"cm")) +
  labs(title="GO enrichment for opa1 vs. E DEGs",
       x="Count",y="Description")
print(g)
pdf("go_bubble.pdf",width=14)
print(g)
dev.off()
```

KEGG

```{r kegg,cache=TRUE,echo=FALSE,message=FALSE,warning=FALSE,fig.height=6}
all_gene_kegg <- read.csv("SMEDProteinKeggKO.csv") %>>%
  (.[c(2,4,1,3)])

ekg <- enricher(deg,TERM2GENE=all_gene_kegg,minGSSize=1,pvalueCutoff=0.5,
                qvalueCutoff=0.5) %>>% as.data.frame

kegg_infor <- read.csv("KEGGFunctionName.csv")

index <- match(ekg[,1],kegg_infor[,1])
ekg$Description <- as.character(kegg_infor[index,2])
ekg <- ekg[order(ekg$p.adjust),]

write.csv(ekg,"opa1.vs.E_deg_KEGG_enrich.csv",row.names=F,quote=F)

# d <- ekg[complete.cases(ekg),] %>% arrange(p.adjust) %>%
#   slice(1:20)
d <- as.data.table(ekg)
d <- d[complete.cases(d),][order(p.adjust),head(.SD,20)]
p <- ggplot(d,aes(forcats::fct_inorder(Description),Count)) +
  geom_col(position="dodge",fill="orangered",color="orangered") +
  theme_light() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        plot.margin=unit(c(0,0,0,1.7),"cm")) +
  labs(title="KEGG enrichment for opa1 vs. E DEGs",x="Description")
print(p)
pdf("kegg.pdf",width=14)
print(p)
dev.off()

g <- ggplot(d,aes(Count,forcats::fct_inorder(Description),size=p.adjust)) +
  geom_point() + 
  theme_light() +
  theme(axis.text.x=element_text(angle=90,hjust=0),
        axis.text.y=element_text(angle=-45,hjust=1,vjust=1,size=6),
        title=element_text(size=10,face="bold"),
        plot.margin=unit(c(1,0,0,2),"cm")) +
  labs(title="KEGG enrichment for opa1 vs. E DEGs",
       x="Count",y="Description")
print(g)
pdf("kegg_bubble.pdf",width=14)
print(g)
dev.off()
```


