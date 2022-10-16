library(DESeq2)
library(dplyr)
library(AnnotationHub)
library(ensembldb)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(EnhancedVolcano)

### tumores control vs alto MRP4 (H samples) T2

star_counts_T2 <- read.table("star counts.tab", header=TRUE)

count.matrixT2 <- select(star_counts, c("T1","T2","T3","T5","T6","T7"))
colnames(count.matrixT2) <- c("control_T1","control_T2","control_T3","MRP4_T1", "MRP4_T2","MRP4_T3")
count.matrixT2 <- as.matrix(count.matrixT2)

sample.tableT2 <- as.data.frame(matrix(c("control","control","control","MRP4","MRP4","MRP4"), ncol = 1))
rownames(sample.tableT2) <- c("control_T1","control_T2","control_T3","MRP4_T1", "MRP4_T2","MRP4_T3")
colnames(sample.tableT2) <- c("group")

sample.tableT2$group <- factor(sample.tableT2$group)

levels(sample.tableT2$group) 

#chequeo que las tablas estan bien
all(rownames(sample.tableT2) %in% colnames(count.matrixT2))
all(rownames(sample.tableT2) == colnames(count.matrixT2))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
ddsT2 <- DESeqDataSetFromMatrix(countData=count.matrixT2, 
                               colData=sample.tableT2, 
                               design=~group)

# al objeto dds (que es el dataset) le agregamos los analisis
ddsT2 <- estimateSizeFactors(ddsT2)
ddsT2 <- DESeq(ddsT2)

resT2 <- results(ddsT2)

table(resT2$padj != "NA")
table(resT2$padj < 0.05)

genesT2 <- star_counts_T2$GeneID.Ensembl
genesT2 <- substr(genesT2,1,15)

#le agrego los genes a la tabla de resultados
resT2$geneid <- genesT2

resT2$symbol_id <- mapIds(org.Hs.eg.db, keys = genesT2, keytype = "ENSEMBL", column="SYMBOL")

resT2$entrez_id <-  mapIds(org.Hs.eg.db, keys = genesT2, keytype = "ENSEMBL", column="ENTREZID")

head(resT2)
dim(resT2)

resT2[which(resT2$symbol_id=="ABCC4"),]

#A summary of the results can be generated:
summary(resT2)

EnhancedVolcano(resT2,
                lab = resT2$symbol_id,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Tumores mock vs MRP4',
                subtitle = "5366 DEG (FDR<0,05, log2FC>1)",
                selectLab = c('ABCC4'),
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 1.5,
                labSize = 5,
                col=c('black', 'gray', 'green', 'red3'),
                colAlpha = 1,
                xlim=c(-15,15))

### Visualizing results
#The MA-plot provides a global view of the differential genes, with the log2 
#fold change on the y-axis over the mean of normalized counts:
plotMA(resT2, ylim=c(-4,4)) #en rojo se ven los genes que tienen padj<0.1

#A sorted results table:
resSortT2 <- resT2[order(resT2$padj),]
head(resSortT2) 

signif_padjT2 <- resSortT2[which(resSortT2$padj<0.05 & resSortT2$padj!= "NA"),]
signif_padjT2$UP_DOWN <- ifelse(signif_padjT2$log2FoldChange>1,"UP",ifelse(signif_padjT2$log2FoldChange< -1, "DOWN","..."))
head(signif_padjT2)
dim(signif_padjT2)

signif_padjT2[which(signif_padjT2$symbol_id=="ABCC4"),]

write.csv(signif_padjT2, file = "total signif genes tumores C4vsH7.csv")

signif.genes_T2 <- signif_padjT2$symbol_id[which(signif_padjT2$UP_DOWN!="...")] %>% na.omit()

comunTT2 <- signif.genes_T2[which(signif.genes_T2%in% signif_genesT)]

#A PCA plot and a heatmap of the top genes: we need to use the rld object
rldT2 <- rlog(ddsT2)

plotPCA(rldT2, intgroup="group")

matT2 <- assay(rldT2)[which(resT2$symbol_id%in% signif.genes_T2),]
matT2 <- matT2 - rowMeans(matT2) # hago el z score al restar cada valor menos la media
head(matT2)
pheatmap(matT2, fontsize = 15)

#Functional analysis with clusterProfiler

# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genesT2 <- genesT2

# Run GO enrichment analysis ALL genes
egoT2_BP <- enrichGO(gene = signif_padjT2$geneid[which(signif_padjT2$UP_DOWN!="...")],
                   universe = all_genesT2,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

egoT2_MF <- enrichGO(gene = signif_padjT2$geneid[which(signif_padjT2$UP_DOWN!="...")],
                   universe = all_genesT2,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

egoT2_CC <- enrichGO(gene = signif_padjT2$geneid[which(signif_padjT2$UP_DOWN!="...")],
                   universe = all_genesT2,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

# Output results from GO analysis to a table
GO_BP.T2 <- data.frame(egoT2_BP) %>% mutate(GO=rep("BP", nrow(egoT2_BP)))
GO_MF.T2 <- data.frame(egoT2_MF)%>% mutate(GO=rep("MF", nrow(egoT2_MF)))
GO_CC.T2 <- data.frame(egoT2_CC) %>% mutate(GO=rep("CC", nrow(egoT2_CC)))

GO_BP.T2[1:50,2]
GO_MF.T2[1:50,2]
GO_CC.T2[1:50,2]

dotplot(egoT2_BP, showCategory=20)
dotplot(egoT2_MF, showCategory=20)
dotplot(egoT2_CC, showCategory=20)

#para achicar las categorias
drop_egoT2_BP <- gofilter(ego_BP.T2, level = 4)

GO_BP_dropT2 <- data.frame(drop_egoT2_BP)

dotplot(drop_egoT2_BP, showCategory=20)

#armo tabla de todos los resultados GO

GOresultsT2 <- rbind(GO_BP.T2, GO_MF.T2,GO_CC.T2)

write.csv(GOresultsT2, file="GO results tumores C4vsH7.csv")

# Run GO enrichment analysis UP genes
egoT2_BP.UP <- enrichGO(gene = signif_padjT2$geneid[which(signif_padjT2$UP_DOWN=="UP")],
                      universe = all_genesT2,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

egoT2_MF.UP <- enrichGO(gene = signif_padjT2$geneid[which(signif_padjT2$UP_DOWN=="UP")],
                      universe = all_genesT2,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db,
                      ont = "MF",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

egoT2_CC.UP <- enrichGO(gene = signif_padjT2$geneid[which(signif_padjT2$UP_DOWN=="UP")],
                      universe = all_genesT2,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

# Output results from GO analysis to a table
GO_BP.T2UP <- data.frame(egoT2_BP.UP) %>% mutate(GO=rep("BP", nrow(egoT2_BP.UP)))
GO_MF.T2UP <- data.frame(egoT2_MF.UP)%>% mutate(GO=rep("MF", nrow(egoT2_MF.UP)))
GO_CC.T2UP <- data.frame(egoT2_CC.UP) %>% mutate(GO=rep("CC", nrow(egoT2_CC.UP)))

GO_BP.T2UP[1:50,2]
GO_MF.T2UP[1:50,2]
GO_CC.T2UP[1:50,2]

dotplot(egoT2_BP.UP, showCategory=20)
dotplot(egoT2_MF.UP, showCategory=20)
dotplot(egoT2_CC.UP, showCategory=20)

#para achicar las categorias
drop_egoT2_BP.UP <- gofilter(egoT2_BP.UP, level = 4)

GO_BP_drop.T2UP <- data.frame(drop_egoT2_BP.UP)

dotplot(drop_egoT2_BP.UP, showCategory=20)

#armo tabla de todos los resultados GO

GOresultsT2.UP <- rbind(GO_BP.T2UP, GO_MF.T2UP,GO_CC.T2UP)

write.csv(GOresultsT2.UP, file="GO results tumores UP C4vsH7.csv")

# Run GO enrichment analysis DOWN genes
egoT2_BP.DOWN <- enrichGO(gene = signif_padjT2$geneid[which(signif_padjT2$UP_DOWN=="DOWN")],
                        universe = all_genesT2,
                        keyType = "ENSEMBL",
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        readable = TRUE)

egoT2_MF.DOWN <- enrichGO(gene = signif_padjT2$geneid[which(signif_padjT2$UP_DOWN=="DOWN")],
                        universe = all_genesT2,
                        keyType = "ENSEMBL",
                        OrgDb = org.Hs.eg.db,
                        ont = "MF",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        readable = TRUE)

egoT2_CC.DOWN <- enrichGO(gene = signif_padjT2$geneid[which(signif_padjT2$UP_DOWN=="DOWN")],
                        universe = all_genesT2,
                        keyType = "ENSEMBL",
                        OrgDb = org.Hs.eg.db,
                        ont = "CC",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        readable = TRUE)

# Output results from GO analysis to a table
GO_BP.T2DOWN <- data.frame(egoT2_BP.DOWN) %>% mutate(GO=rep("BP", nrow(egoT2_BP.DOWN)))
GO_MF.T2DOWN <- data.frame(egoT2_MF.DOWN)%>% mutate(GO=rep("MF", nrow(egoT2_MF.DOWN)))
GO_CC.T2DOWN <- data.frame(egoT2_CC.DOWN) %>% mutate(GO=rep("CC", nrow(egoT2_CC.DOWN)))

GO_BP.T2DOWN[1:50,]
GO_MF.T2DOWN[1:50,]
GO_CC.T2DOWN[1:50,]

dotplot(egoT2_BP.DOWN, showCategory=20)
dotplot(egoT2_MF.DOWN, showCategory=20)
dotplot(egoT2_CC.DOWN, showCategory=20)

#para achicar las categorias
drop_egoT2_BP.DOWN <- gofilter(egoT2_BP.DOWN, level = 4)

GO_BPdrop.T2DOWN <- data.frame(drop_egoT2_BP.DOWN)

dotplot(drop_egoT2_BP.DOWN, showCategory=20)

#armo tabla de todos los resultados GO

GOresultsT2.DOWN <- rbind(GO_BP.T2DOWN, GO_MF.T2DOWN,GO_CC.T2DOWN)

write.csv(GOresultsT2.DOWN, file="GO results cultivos DOWN C4vsH7.csv")
