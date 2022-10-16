library(DESeq2)
library(tidyverse)
library(AnnotationHub)
library(ensembldb)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(EnhancedVolcano)
library(DOSE)

### tumores control vs alto MRP4 (H samples) pasados por xenofilter T

star_counts_T <- read.delim("Counts de Homo Xenofiltereadas.txt", header=TRUE)

count.matrixT <- select(star_counts_T, c("T1","T2","T3","T5","T6","T7"))
colnames(count.matrixT) <- c("control_T1","control_T2","control_T3","MRP4_T1", "MRP4_T2","MRP4_T3")
count.matrixT <- as.matrix(count.matrixT)

sample.tableT <- as.data.frame(matrix(c("control","control","control","MRP4","MRP4","MRP4"), ncol = 1))
rownames(sample.tableT) <- c("control_T1","control_T2","control_T3","MRP4_T1", "MRP4_T2","MRP4_T3")
colnames(sample.tableT) <- c("group")

sample.tableT$group <- factor(sample.tableT$group)

levels(sample.tableT$group) 

#chequeo que las tablas estan bien
all(rownames(sample.tableT) %in% colnames(count.matrixT))
all(rownames(sample.tableT) == colnames(count.matrixT))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
ddsT <- DESeqDataSetFromMatrix(countData=count.matrixT, 
                               colData=sample.tableT, 
                               design=~group)

# al objeto dds (que es el dataset) le agregamos los analisis
ddsT <- estimateSizeFactors(ddsT)
ddsT <- DESeq(ddsT)

resT <- results(ddsT)

table(resT$padj != "NA")
table(resT$padj < 0.05)

#le agrego los genes a la tabla de resultados
resT$geneid <- star_counts_T$Ensembl.Gene.ID

genes <- star_counts_T$Ensembl.Gene.ID

resT$symbol_id <- mapIds(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", column="SYMBOL")

resT$entrez_id <-  mapIds(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", column="ENTREZID")

head(resT)
dim(resT)

resT[which(resT$symbol_id=="GLI3"),]

#A summary of the results can be generated:
summary(resT)

EnhancedVolcano(resT,
                lab = resT$symbol_id,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Tumores mock vs MRP4',
                subtitle = "xxxx DEG (FDR<0,05, log2FC>0,5)",
                selectLab = c('ABCC4'),
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 1.5,
                labSize = 5,
                col=c('black', 'gray', 'green', 'red3'),
                colAlpha = 1,
                xlim=c(-15,15))

### Visualizing results
#The MA-plot provides a global view of the differential genes, with the log2 
#fold change on the y-axis over the mean of normalized counts:
plotMA(resT, ylim=c(-4,4)) #en rojo se ven los genes que tienen padj<0.1

#A sorted results table:
resSortT <- resT[order(resT$padj),]
head(resSortT) 

signif_padjT <- resSortT[which(resSortT$padj<0.05 & resSortT$padj!= "NA"),]
signif_padjT$UP_DOWN <- ifelse(signif_padjT$log2FoldChange>1,"UP",ifelse(signif_padjT$log2FoldChange< -1, "DOWN","..."))
head(signif_padjT)
dim(signif_padjT)

signif_padjT[which(signif_padjT$symbol_id=="KDM1A"),]

signif_genesT <- signif_padjT$symbol_id[which(signif_padjT$UP_DOWN!="...")]
signif_genesT <- unique(signif_genesT) %>% unlist() %>% na.omit()

signif_genesTUP <- signif_padjT$symbol_id[which(signif_padjT$UP_DOWN=="UP")] 
signif_genesTUP <- unique(signif_genesTUP) %>% unlist() %>% na.omit()
signif_genesT.UP <- cbind(signif_genesTUP, UP.DOWN=rep("UP",length(signif_genesTUP)))

signif_genesTDOWN <- signif_padjT$symbol_id[which(signif_padjT$UP_DOWN=="DOWN")] 
signif_genesTDOWN <- unique(signif_genesTDOWN) %>% unlist() %>% na.omit()
signif_genesT.DOWN <- cbind(signif_genesTDOWN, UP.DOWN=rep("DOWN",length(signif_genesTDOWN)))

write.csv(rbind(signif_genesT.UP,signif_genesT.DOWN), file = "total signif genes tumores C4 vs H7 xenofilter.csv")

which(signif_genesT=="MDFIC")

#A PCA plot and a heatmap of the top genes: we need to use the rld object
rldT <- rlog(ddsT)

plotPCA(rldT, intgroup="group")

matT <- assay(rldT)[which(resT$symbol_id %in% signif_genesT),]
matT <- matT - rowMeans(matT) # hago el z score al restar cada valor menos la media
head(matT)
pheatmap(matT, fontsize = 15)

#Functional analysis with clusterProfiler

# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genesT <- genes

# Run GO enrichment analysis ALL genes
egoT_BP <- enrichGO(gene = signif_padjT$geneid[which(signif_padjT$UP_DOWN!="...")],
                 universe = all_genesT,
                 keyType = "ENSEMBL",
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 readable = TRUE)

egoT_MF <- enrichGO(gene = signif_padjT$geneid[which(signif_padjT$UP_DOWN!="...")],
                    universe = all_genesT,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "MF",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

egoT_CC <- enrichGO(gene = signif_padjT$geneid[which(signif_padjT$UP_DOWN!="...")],
                    universe = all_genesT,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "CC",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

# Output results from GO analysis to a table
GO_BP_T <- data.frame(egoT_BP) %>% mutate(GO=rep("BP", nrow(egoT_BP)))
GO_MF_T <- data.frame(egoT_MF)%>% mutate(GO=rep("MF", nrow(egoT_MF)))
GO_CC_T <- data.frame(egoT_CC) %>% mutate(GO=rep("CC", nrow(egoT_CC)))

GO_BP_T[1:50,]
GO_MF_T[1:50,]
GO_CC_T[1:50,]

dotplot(egoT_BP, showCategory=20)
dotplot(egoT_MF, showCategory=20)
dotplot(egoT_CC, showCategory=20)

#para achicar las categorias
drop_egoT_BP <- gofilter(egoT_BP, level = 4)

GO_BP_Tdrop <- data.frame(drop_egoT_BP)

dotplot(drop_egoT_BP, showCategory=20)

#armo tabla de todos los resultados GO

GOresultsT <- rbind(GO_BP_T, GO_MF_T,GO_CC_T)

write.csv(GOresultsT, file="GO results tumores HUMAN xenof.csv")

# Run GO enrichment analysis UP genes
egoT_BP.UP <- enrichGO(gene = signif_padjT$geneid[which(signif_padjT$UP_DOWN=="UP")],
                    universe = all_genesT,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

egoT_MF.UP <- enrichGO(gene =signif_padjT$geneid[which(signif_padjT$UP_DOWN=="UP")],
                    universe = all_genesT,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "MF",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

egoT_CC.UP <- enrichGO(gene = signif_padjT$geneid[which(signif_padjT$UP_DOWN=="UP")],
                    universe = all_genesT,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "CC",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

# Output results from GO analysis to a table
GO_BP_T.UP <- data.frame(egoT_BP.UP) %>% mutate(GO=rep("BP", nrow(egoT_BP.UP)))
GO_MF_T.UP <- data.frame(egoT_MF.UP)%>% mutate(GO=rep("MF", nrow(egoT_MF.UP)))
GO_CC_T.UP <- data.frame(egoT_CC.UP) %>% mutate(GO=rep("CC", nrow(egoT_CC.UP)))

GO_BP_T.UP[1:50,]
GO_MF_T.UP[1:50,]
GO_CC_T.UP[1:50,]

dotplot(egoT_BP.UP, showCategory=20)
dotplot(egoT_MF.UP, showCategory=20)
dotplot(egoT_CC.UP, showCategory=20)

#para achicar las categorias
drop_egoT_BP.UP <- gofilter(egoT_BP.UP, level = 4)
GO_BP_Tdrop.UP <- data.frame(drop_egoT_BP.UP)

drop_egoT_MF.UP <- gofilter(egoT_MF.UP, level = 4)
GO_MF_Tdrop.UP <- data.frame(drop_egoT_MF.UP)

drop_egoT_CC.UP <- gofilter(egoT_CC.UP, level = 4)
GO_CC_Tdrop.UP <- data.frame(drop_egoT_CC.UP)

dotplot(drop_egoT_BP.UP, showCategory=20)
dotplot(drop_egoT_MF.UP, showCategory=20)
dotplot(drop_egoT_CC.UP, showCategory=20)

#armo tabla de todos los resultados GO

GOresultsT.UP <- rbind(GO_BP_T.UP, GO_MF_T.UP,GO_CC_T.UP)
GOresultsT.UP <- rbind(GO_BP_Tdrop.UP, GO_MF_Tdrop.UP,GO_CC_Tdrop.UP)

write.csv(GOresultsT.UP, file="GO results tumores HUMAN xenof UP.csv")

# Run GO enrichment analysis DOWN genes
egoT_BP.DOWN <- enrichGO(gene = signif_padjT$geneid[which(signif_padjT$UP_DOWN=="DOWN")],
                       universe = all_genesT,
                       keyType = "ENSEMBL",
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)

egoT_MF.DOWN <- enrichGO(gene =signif_padjT$geneid[which(signif_padjT$UP_DOWN=="DOWN")],
                       universe = all_genesT,
                       keyType = "ENSEMBL",
                       OrgDb = org.Hs.eg.db,
                       ont = "MF",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)

egoT_CC.DOWN <- enrichGO(gene = signif_padjT$geneid[which(signif_padjT$UP_DOWN=="DOWN")],
                       universe = all_genesT,
                       keyType = "ENSEMBL",
                       OrgDb = org.Hs.eg.db,
                       ont = "CC",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)

# Output results from GO analysis to a table
GO_BP_T.DOWN <- data.frame(egoT_BP.DOWN) %>% mutate(GO=rep("BP", nrow(egoT_BP.DOWN)))
GO_MF_T.DOWN <- data.frame(egoT_MF.DOWN)%>% mutate(GO=rep("MF", nrow(egoT_MF.DOWN)))
GO_CC_T.DOWN <- data.frame(egoT_CC.DOWN) %>% mutate(GO=rep("CC", nrow(egoT_CC.DOWN)))

GO_BP_T.DOWN[1:50,]
GO_MF_T.DOWN[1:50,]
GO_CC_T.DOWN[1:50,]

dotplot(egoT_BP.DOWN, showCategory=20)
dotplot(egoT_MF.DOWN, showCategory=20)
dotplot(egoT_CC.DOWN, showCategory=20)

#para achicar las categorias
drop_egoT_BP.DOWN <- gofilter(egoT_BP.DOWN, level = 4)
GO_BP_Tdrop.DOWN <- data.frame(drop_egoT_BP.DOWN)

dotplot(drop_egoT_BP.DOWN, showCategory=20)

#armo tabla de todos los resultados GO

GOresultsT.DOWN <- rbind(GO_BP_T.DOWN, GO_MF_T.DOWN,GO_CC_T.DOWN)

write.csv(GOresultsT.DOWN, file="GO results tumores HUMAN xenof DOWN.csv")
