library(DESeq2)
library(dplyr)
library(AnnotationHub)
library(AnnotationDbi)
library(ensembldb)
library(org.Mm.eg.db)
library(clusterProfiler)
library(pathview)
library(pheatmap)
library(EnhancedVolcano)


### cultivos control vs alto MRP4 (H samples) counts de raton xenofilter TM

star_counts_Tmouse <- read.delim("MouseCountsXenofiltered.txt", header=TRUE)

count.matrixTmouse <- select(star_counts_Tmouse, c("T1","T2","T3","T5","T6","T7"))
colnames(count.matrixTmouse) <- c("control_T1","control_T2","control_T3","MRP4_T1", "MRP4_T2","MRP4_T3")
count.matrixTmouse <- as.matrix(count.matrixTmouse)

sample.tableTmouse <- as.data.frame(matrix(c("control","control","control","MRP4","MRP4","MRP4"), ncol = 1))
rownames(sample.tableTmouse) <- c("control_T1","control_T2","control_T3","MRP4_T1", "MRP4_T2","MRP4_T3")
colnames(sample.tableTmouse) <- c("group")

sample.tableTmouse$group <- factor(sample.tableTmouse$group)

levels(sample.tableTmouse$group) 

#chequeo que las tablas estan bien
all(rownames(sample.tableTmouse) %in% colnames(count.matrixTmouse))
all(rownames(sample.tableTmouse) == colnames(count.matrixTmouse))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
ddsTM <- DESeqDataSetFromMatrix(countData=count.matrixTmouse, 
                               colData=sample.tableTmouse, 
                               design=~group)

# al objeto dds (que es el dataset) le agregamos los analisis
ddsTM <- estimateSizeFactors(ddsTM)
ddsTM <- DESeq(ddsTM)

resTM <- results(ddsTM)

table(resTM$padj != "NA")
table(resTM$padj < 0.1)

#le agrego los genes a la tabla de resultados
genesTM <- star_counts_Tmouse$Ensembl.Mouse.Transcripts.ID

resTM$geneid <- genesTM

resTM$symbol_id <- mapIds(org.Mm.eg.db, keys = genesTM, keytype = "ENSEMBL", column="SYMBOL")

resTM$entrez_id <- mapIds(org.Mm.eg.db, keys = genesTM, keytype = "ENSEMBL", column="ENTREZID")

head(resTM)
dim(resTM)

resTM[which(resTM$symbol_id=="Xkr4"),]

#A summary of the results can be generated:
summary(resTM)

EnhancedVolcano(resTM,
                lab = resT$symbol_id,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Tumores mock vs MRP4: estroma mouse',
                subtitle = "xxxx DEG (FDR<0,05, log2FC>0,5)",
                selectLab = c('ABCC4'),
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 1.5,
                labSize = 5,
                col=c('black', 'gray', 'green', 'red3'),
                colAlpha = 1,
                xlim=c(-10,10))

### Visualizing results
#The MA-plot provides a global view of the differential genes, with the log2 
#fold change on the y-axis over the mean of normalized counts:
plotMA(resTM, ylim=c(-4,4)) #en rojo se ven los genes que tienen padj<0.1

#A sorted results table:
resSortTM <- resTM[order(resTM$padj),]
head(resSortTM) 

signif_padjTM <- resSortTM[which(resSortTM$padj<0.1 & resSortTM$padj!= "NA"),]
signif_padjTM$UP_DOWN <- ifelse(signif_padjTM$log2FoldChange>0.5,"UP",ifelse(signif_padjTM$log2FoldChange< -0.5, "DOWN","..."))
head(signif_padjTM)
dim(signif_padjT)

signif_padjTM[which(signif_padjTM$symbol_id=="Nprl3"),]

signif_genesTM <- signif_padjTM$symbol_id
signif_genesTM <- unique(signif_genesTM) %>% unlist()
signif_genesTM <- na.omit(signif_genesTM)

write.csv(signif_genesTM, file = "total signif genes tumores C4 vs H7 RATON.csv")

#A PCA plot and a heatmap of the top genes: we need to use the rld object
rldTM <- rlog(ddsTM)

plotPCA(rldTM, intgroup="group")

matTM <- assay(rldTM)[which(resTM$padj<0.1 & resTM$padj!="NA"),]
matTM <- matTM - rowMeans(matTM) # hago el z score al restar cada valor menos la media
head(matTM)

pheatmap(matTM, fontsize = 15, show_rownames = FALSE)

#Functional analysis with clusterProfiler

# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genesTM <- genesTM

# Run GO enrichment analysis ALL genes
egoTM_BP <- enrichGO(gene = signif_padjTM$geneid,
                    universe = all_genesTM,
                    keyType = "ENSEMBL",
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

egoTM_MF <- enrichGO(gene = signif_padjTM$geneid,
                    universe = all_genesTM,
                    keyType = "ENSEMBL",
                    OrgDb = org.Mm.eg.db,
                    ont = "MF",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

egoTM_CC <- enrichGO(gene = signif_padjTM$geneid,
                    universe = all_genesTM,
                    keyType = "ENSEMBL",
                    OrgDb = org.Mm.eg.db,
                    ont = "CC",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

# Output results from GO analysis to a table
GO_BP_TM <- data.frame(egoTM_BP) %>% mutate(GO=rep("BP", nrow(egoTM_BP)))
GO_MF_TM <- data.frame(egoTM_MF)%>% mutate(GO=rep("MF", nrow(egoTM_MF)))
GO_CC_TM <- data.frame(egoTM_CC) %>% mutate(GO=rep("CC", nrow(egoTM_CC)))

GO_BP_TM[1:50,]
GO_MF_TM[1:50,]
GO_CC_TM[1:50,]

dotplot(egoTM_BP, showCategory=20)
dotplot(egoTM_MF, showCategory=20)
dotplot(egoTM_CC, showCategory=20)

#armo tabla de todos los resultados GO

GOresultsTM <- rbind(GO_BP_TM, GO_MF_TM,GO_CC_TM)

write.csv(GOresultsTM, file="GO results tumores MOUSE xenof.csv")

# Extract significant results
signif_resTMUP <- signif_padjTM[signif_padjTM$log2FoldChange> 0.5,]
signif_genesTM_UP <- signif_resTMUP$geneid %>% unique() %>% na.omit()
length(signif_genesTM_UP)

signif_resTMDOWN <- signif_padjTM[signif_padjTM$log2FoldChange< -0.5,]
signif_genesTM_DOWN <- signif_resTMDOWN$geneid %>% unique() %>% na.omit()
length(signif_genesTM_DOWN)

TM_UP <- cbind(unique(signif_resTMUP$symbol_id), rep("UP", length(unique(signif_resTMUP$symbol_id))))
TM_DOWN <- cbind(unique(signif_resTMDOWN$symbol_id), rep("DOWN", length(unique(signif_resTMDOWN$symbol_id))))

signif_genesTM.UP.DOWN <- rbind(TM_UP,TM_DOWN)

write.csv(signif_genesTM.UP.DOWN, "signif genes T UP DOWN xenof mouse.csv")

# Run GO enrichment analysis UP genes
egoTM_BP.UP <- enrichGO(gene = signif_genesTM_UP,
                       universe = all_genesTM,
                       keyType = "ENSEMBL",
                       OrgDb = org.Mm.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.1,
                       readable = TRUE)

egoTM_MF.UP <- enrichGO(gene =signif_genesTM_UP,
                       universe = all_genesTM,
                       keyType = "ENSEMBL",
                       OrgDb = org.Mm.eg.db,
                       ont = "MF",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.1,
                       readable = TRUE)

egoTM_CC.UP <- enrichGO(gene = signif_genesTM_UP,
                       universe = all_genesTM,
                       keyType = "ENSEMBL",
                       OrgDb = org.Mm.eg.db,
                       ont = "CC",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.1,
                       readable = TRUE)

# Output results from GO analysis to a table
GO_BP_TM.UP <- data.frame(egoTM_BP.UP) %>% mutate(GO=rep("BP", nrow(egoTM_BP.UP)))
GO_MF_TM.UP <- data.frame(egoTM_MF.UP)%>% mutate(GO=rep("MF", nrow(egoTM_MF.UP)))
GO_CC_TM.UP <- data.frame(egoTM_CC.UP) %>% mutate(GO=rep("CC", nrow(egoTM_CC.UP)))

GO_BP_TM.UP[1:50,]
GO_MF_TM.UP[1:50,]
GO_CC_TM.UP[1:50,]

dotplot(egoTM_BP.UP, showCategory=20)
dotplot(egoTM_MF.UP, showCategory=20)
dotplot(egoTM_CC.UP, showCategory=20)

#armo tabla de todos los resultados GO

GOresultsTM.UP <- rbind(GO_BP_TM.UP, GO_MF_TM.UP,GO_CC_TM.UP)

write.csv(GOresultsTM.UP, file="GO results tumores MOUSE xenof UP.csv")

# Run GO enrichment analysis DOWN genes
egoTM_BP.DOWN <- enrichGO(gene = signif_genesTM_DOWN,
                         universe = all_genesTM,
                         keyType = "ENSEMBL",
                         OrgDb = org.Mm.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

egoTM_MF.DOWN <- enrichGO(gene =signif_genesTM_DOWN,
                         universe = all_genesTM,
                         keyType = "ENSEMBL",
                         OrgDb = org.Mm.eg.db,
                         ont = "MF",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

egoTM_CC.DOWN <- enrichGO(gene = signif_genesTM_DOWN,
                         universe = all_genesTM,
                         keyType = "ENSEMBL",
                         OrgDb = org.Mm.eg.db,
                         ont = "CC",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

# Output results from GO analysis to a table
GO_BP_TM.DOWN <- data.frame(egoTM_BP.DOWN) %>% mutate(GO=rep("BP", nrow(egoTM_BP.DOWN)))
GO_MF_TM.DOWN <- data.frame(egoTM_MF.DOWN)%>% mutate(GO=rep("MF", nrow(egoTM_MF.DOWN)))
GO_CC_TM.DOWN <- data.frame(egoTM_CC.DOWN) %>% mutate(GO=rep("CC", nrow(egoTM_CC.DOWN)))

GO_BP_TM.DOWN[1:50,]
GO_MF_TM.DOWN[1:50,]
GO_CC_TM.DOWN[1:50,]

dotplot(egoTM_BP.DOWN, showCategory=20)
dotplot(egoTM_MF.DOWN, showCategory=20)
dotplot(egoTM_CC.DOWN, showCategory=20)

#armo tabla de todos los resultados GO

GOresultsTM.DOWN <- rbind(GO_BP_TM.DOWN, GO_MF_TM.DOWN,GO_CC_TM.DOWN)

write.csv(GOresultsTM.DOWN, file="GO results tumores MOUSE xenof DOWN.csv")








#################################################################################
# gene set enrichment analysis (GSEA) using clusterProfiler and Pathview

genelistTM <- signif_padjTM$entrez_id
genelistTM <- genelistTM[which(genelistTM!="NA")] %>% unlist()
genelist <- as.numeric(genelist)


KEGG <- enrichKEGG(
  genelistTM,
  organism = "mmu",
  keyType = "kegg",
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH",
  universe,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)
  
# Extract the GSEA results
KEGG_results <- KEGG@result

#################################################

DAVID <- enrichDAVID(
  genelist,
  idType = "ENTREZ_GENE_ID",
  universe,
  minGSSize = 10,
  maxGSSize = 500,
  annotation = "GOTERM_BP_FAT",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  species = NA,
  david.user
)
