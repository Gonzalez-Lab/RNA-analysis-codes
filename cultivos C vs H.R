library(DESeq2)
library(dplyr)
library(AnnotationHub)
library(ensembldb)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(pheatmap)
library(EnhancedVolcano)

### cultivos control vs alto MRP4 (H samples)

star_counts <- read.table("star counts.tab", header=TRUE)

count.matrix <- select(star_counts, c("C2","C4","H1","H2"))
colnames(count.matrix) <- c('control_1',"control_2","MRP4_1","MRP4_2")
count.matrix <- as.matrix(count.matrix)

sample.table <- as.data.frame(matrix(c("control","control","MRP4","MRP4"), ncol = 1))
rownames(sample.table) <- c("control_1","control_2","MRP4_1","MRP4_2")
colnames(sample.table) <- c("group")

sample.table$group <- factor(sample.table$group)

levels(sample.table$group) 

#chequeo que las tablas estan bien
all(rownames(sample.table) %in% colnames(count.matrix))
all(rownames(sample.table) == colnames(count.matrix))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
dds <- DESeqDataSetFromMatrix(countData=count.matrix, 
                              colData=sample.table, 
                              design=~group)

# al objeto dds (que es el dataset) le agregamos los analisis
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

res <- results(dds)

table(res$padj != "NA")
table(res$padj < 0.05)
table(res$padj < 0.05 & res$log2FoldChange< -0.5) # downregulados
table(res$padj < 0.05 & res$log2FoldChange> 0.5) # upregulados

#le agrego los genes a la tabla de resultados
res$geneid <- star_counts$GeneID.Ensembl

# remuevo los ultimos tres caracteres de ensembl porque no lo lee mapIds sino
genes <- star_counts$GeneID.Ensembl
genes <- substr(genes,1,15)

#le agrego los genes a la tabla de resultados
res$geneid <- genes

res$symbol_id <- mapIds(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", column="SYMBOL")

res$entrez_id <-  mapIds(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", column="ENTREZID")

head(res)
dim(res)

res[which(res$symbol_id=="ABCC4"),]

#A summary of the results can be generated:
summary(res)

EnhancedVolcano(res,
                lab = res$symbol_id,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Cultivos mock vs MRP4',
                subtitle = "1593 DEGs (FDR<0,05, log2FC>1)",
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
plotMA(res, ylim=c(-4,4)) #en rojo se ven los genes que tienen padj<0.1

#A sorted results table:
resSort <- res[order(res$padj),]
head(resSort) 

signif_padj <- resSort[which(resSort$padj<0.05 & resSort$padj!= "NA"),]
signif_padj$UP_DOWN <- ifelse(signif_padj$log2FoldChange>1,"UP",ifelse(signif_padj$log2FoldChange< -1,"DOWN","..."))
head(signif_padj)
dim(signif_padj)

signif.genes_C <- signif_padj$symbol_id[which(signif_padj$UP_DOWN!="...")] %>% na.omit()
signif.genes_CUP <- signif_padj$symbol_id[which(signif_padj$UP_DOWN=="UP")] %>% na.omit()
signif.genes_CDOWN <- signif_padj$symbol_id[which(signif_padj$UP_DOWN=="DOWN")] %>% na.omit()

comunCT <- signif.genes_C[which(signif.genes_C %in% signif_genesT)]

write.csv(signif_padj, file = "total signif genes cultivos C4 vs H7.csv")

#A PCA plot and a heatmap of the top genes: we need to use the rld object
rld <- rlog(dds)

plotPCA(rld, intgroup="group")

mat <- assay(rld)[which(res$symbol_id %in% signif.genes_C),]
mat <- mat - rowMeans(mat) # hago el z score al restar cada valor menos la media
head(mat)
pheatmap(mat, fontsize = 15)

#Functional analysis with clusterProfiler

# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- genes

# Run GO enrichment analysis ALL genes
ego_BP <- enrichGO(gene = signif_padj$geneid[which(signif_padj$UP_DOWN!="...")],
                    universe = all_genes,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

ego_MF <- enrichGO(gene = signif_padj$geneid[which(signif_padj$UP_DOWN!="...")],
                    universe = all_genes,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "MF",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

ego_CC <- enrichGO(gene = signif_padj$geneid[which(signif_padj$UP_DOWN!="...")],
                    universe = all_genes,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "CC",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

# Output results from GO analysis to a table
GO_BP <- data.frame(ego_BP) %>% mutate(GO=rep("BP", nrow(ego_BP)))
GO_MF <- data.frame(ego_MF)%>% mutate(GO=rep("MF", nrow(ego_MF)))
GO_CC <- data.frame(ego_CC) %>% mutate(GO=rep("CC", nrow(ego_CC)))

GO_BP[1:50,]
GO_MF[1:50,]
GO_CC[1:50,]

dotplot(ego_BP, showCategory=20)
dotplot(ego_MF, showCategory=20)
dotplot(ego_CC, showCategory=20)

#para achicar las categorias
drop_ego_BP <- gofilter(ego_BP, level = 4)

GO_BP_drop <- data.frame(drop_ego_BP)

dotplot(drop_ego_BP, showCategory=20)

#armo tabla de todos los resultados GO

GOresults <- rbind(GO_BP, GO_MF,GO_CC)

write.csv(GOresults, file="GO results cultivos C4vsH7.csv")

# Run GO enrichment analysis UP genes
ego_BP.UP <- enrichGO(gene = signif_padj$geneid[which(signif_padj$UP_DOWN=="UP")],
                       universe = all_genes,
                       keyType = "ENSEMBL",
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)

ego_MF.UP <- enrichGO(gene = signif_padj$geneid[which(signif_padj$UP_DOWN=="UP")],
                       universe = all_genes,
                       keyType = "ENSEMBL",
                       OrgDb = org.Hs.eg.db,
                       ont = "MF",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)

ego_CC.UP <- enrichGO(gene = signif_padj$geneid[which(signif_padj$UP_DOWN=="UP")],
                       universe = all_genes,
                       keyType = "ENSEMBL",
                       OrgDb = org.Hs.eg.db,
                       ont = "CC",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)

# Output results from GO analysis to a table
GO_BP.UP <- data.frame(ego_BP.UP) %>% mutate(GO=rep("BP", nrow(ego_BP.UP)))
GO_MF.UP <- data.frame(ego_MF.UP)%>% mutate(GO=rep("MF", nrow(ego_MF.UP)))
GO_CC.UP <- data.frame(ego_CC.UP) %>% mutate(GO=rep("CC", nrow(ego_CC.UP)))

GO_BP.UP[1:50,]
GO_MF.UP[1:50,]
GO_CC.UP[1:50,]

dotplot(ego_BP.UP, showCategory=20)
dotplot(ego_MF.UP, showCategory=20)
dotplot(ego_CC.UP, showCategory=20)

#para achicar las categorias
drop_ego_BP.UP <- gofilter(ego_BP.UP, level = 4)

GO_BP_drop.UP <- data.frame(drop_ego_BP.UP)

dotplot(drop_ego_BP.UP, showCategory=20)

#armo tabla de todos los resultados GO

GOresults.UP <- rbind(GO_BP.UP, GO_MF.UP,GO_CC.UP)

write.csv(GOresults.UP, file="GO results cultivos UP C4vsH7.csv")

# Run GO enrichment analysis DOWN genes
ego_BP.DOWN <- enrichGO(gene = signif_padj$geneid[which(signif_padj$UP_DOWN=="DOWN")],
                         universe = all_genes,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

ego_MF.DOWN <- enrichGO(gene = signif_padj$geneid[which(signif_padj$UP_DOWN=="DOWN")],
                         universe = all_genes,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db,
                         ont = "MF",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

ego_CC.DOWN <- enrichGO(gene = signif_padj$geneid[which(signif_padj$UP_DOWN=="DOWN")],
                         universe = all_genes,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db,
                         ont = "CC",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

# Output results from GO analysis to a table
GO_BP.DOWN <- data.frame(ego_BP.DOWN) %>% mutate(GO=rep("BP", nrow(ego_BP.DOWN)))
GO_MF.DOWN <- data.frame(ego_MF.DOWN)%>% mutate(GO=rep("MF", nrow(ego_MF.DOWN)))
GO_CC.DOWN <- data.frame(ego_CC.DOWN) %>% mutate(GO=rep("CC", nrow(ego_CC.DOWN)))

GO_BP.DOWN[1:50,]
GO_MF.DOWN[1:50,]
GO_CC.DOWN[1:50,]

dotplot(ego_BP.DOWN, showCategory=20)
dotplot(ego_MF.DOWN, showCategory=20)
dotplot(ego_CC.DOWN, showCategory=20)

#para achicar las categorias
drop_ego_BP.DOWN <- gofilter(ego_BP.DOWN, level = 4)

GO_BPdrop.DOWN <- data.frame(drop_ego_BP.DOWN)

dotplot(drop_ego_BP.DOWN, showCategory=20)

#armo tabla de todos los resultados GO

GOresults.DOWN <- rbind(GO_BP.DOWN, GO_MF.DOWN,GO_CC.DOWN)

write.csv(GOresults.DOWN, file="GO results cultivos DOWN C4vsH7.csv")
