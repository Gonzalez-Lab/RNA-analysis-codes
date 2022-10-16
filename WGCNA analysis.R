library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

count_matrix_T <- count.matrixT
rownames(count_matrix_T) <- res$symbol_id
count_matrix_T <- count_matrix_T[which(rownames(count_matrix_T)!="NA"),]

datExpr <- t(count_matrix_T) # samples in row, genes in column

# Calculate sample distance and cluster the samples
sampleTree = hclust(dist(datExpr), method = "average")

# plot sample tree
par(cex = 0.5)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1,cex.axis = 1, cex.main = 1)

# Plot a line to show the cut
abline(h = 1.25e+06, col = "red")

# Choose a set of soft threshold parameters
powers = c(c(1:20), seq(from = 22, to=30, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 

# Scale-free topology fit index as a function of the soft-thresholding power
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#===============================================================================
#
#  Turn data expression into topological overlap matrix
#
#===============================================================================

# Turn data expression into topological overlap matrix
power=sft$powerEstimate

adjacency = adjacency(datExpr, power = power)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency) # TOM topological overlap matrix
dissTOM = 1-TOM # dissimilarity for clustering

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

par(mfrow = c(1,1))
# Plot the resulting clustering tree (dendrogram)
plot(geneTree, 
     xlab="", sub="", 
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, 
                            distM = dissTOM,
                            deepSplit = 2, 
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, 
                    dynamicColors, 
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
plot(METree, 
     main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.5
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(datExpr, 
                          dynamicColors, 
                          cutHeight = MEDissThres, 
                          verbose = 3)

# The merged module colors
mergedColors = merge$colors

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

plotDendroAndColors(geneTree, 
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

tabla_mod <- as.data.frame(cbind(colnames(datExpr), dynamicColors, mergedColors))
colnames(tabla_mod) = c("gene","clust","mergedclust")

write.csv(tabla_mod,"modulos genes PAAD alto vs bajo MRP4.csv")

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs