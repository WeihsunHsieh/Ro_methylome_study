### Weighted Gene Co-expression Network Analysis (WGCNA) ###
 # WGCNA was performed in R using normalized expression values derived from DESeq2.

# 1. Data Import: Load RNA-seq data into DESeq2 with DEGs filtered (FDR < 0.05, fold change > 2), and import the Anthocyanin data accordingly.

RoData=read.csv("DEseq2NormCounts.txt",sep="\t",header=TRUE);
RoPheno=read.csv("Anthocyanins.txt",sep="\t",header=TRUE);
RoPhenoSpec7=read.csv("metabolites.txt",sep="\t",header=TRUE);

# 2. Transpose matrices (samples × genes)

datExpr0=as.data.frame(t(RoData[,-c(1:1)]));
names(datExpr0)=RoData$ID;
rownames(datExpr0)=names(RoData)[-c(1:1)];

## a. Removing Outlier Genes
gsg=goodSamplesGenes(datExpr0, verbose=3);
gsg$allOK

if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse=",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing sample:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=",")));
  datExpr0=datExpr0[gsg$goodSamples, gsg$goodGenes]
}

## b. Removing Outlier Samples (Outliers were excluded using hierarchical clustering)
sampleTree=hclust(dist(datExpr0), method='average');
sizeGrWindow(12,9)
#pdf(file="stage_hclust.pdf", width = 12, height = 9) 
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
#dev.off()

## c. Retaining the Rest of Samples
clust=cutreeStatic(sampleTree, cutHeight=0, minSize=10)
keepSamples=(clust==0)
datExpr=datExpr0[keepSamples,]
nGenes=ncol(datExpr)
nSamples=nrow(datExpr)

## d. Plotting the samples and trait heatmap
diagData=diag(1, nrow=14);
datTraits=as.data.frame(diagData); 
RoSamples = rownames(datExpr); 
names(datTraits)=rownames(datTraits)=RoSamples;
collectGarbage();

sampleTree2=hclust(dist(datExpr), method="average")
sizeGrWindow(12,9)
#pdf(file="stage_ClusteringAndHeatmap.pdf", width = 12, height = 9)
traitColors=numbers2colors(datTraits, signed=FALSE);
plotDendroAndColors(sampleTree2, traitColors, 
                    groupLabels=names(datTraits), 
                    main="Sample dendrogram and trait heatmap", sub="", xlab="")
#dev.off()

# 3. Network Construction (Soft-threshold power selection)

## a. Set network types and use Pearson correlation as the selection criterion
powers = c(c(1:10), seq(from = 12, to = 30, by=2))

### AdjMatToFindPower
similarity.pearson = cor((datExpr), use = 'p', method = 'pearson' );
sft.sm.ps= pickSoftThreshold.fromSimilarity(similarity = similarity.pearson, RsquaredCut = 0.7, powerVector = powers)
save(datExpr0, sampleTree, datExpr, datTraits, sampleTree2, sft.sm.ps, 
     file = "inputdataPower_AdjMat.RData")

## b. Plot the scale-free topology fit index against soft-thresholding power to identify the optimal β for downstream analysis.
sft <- sft.sm.ps
#pdf(file="ChoosingSoftThresholdPower_AdjMat_R0pt7.pdf", width = 12, height = 9)
par(mfrow = c(1,2));
cex1 = 0.70;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2", type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1, col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90, col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
#dev.off()

# 4. Use the optimal β to construct the adjacency matrix and Topological Overlap Matrix (TOM) to assess gene–gene similarity.

## a. adjacency: Correlation values are raised to the optimal soft-thresholding power (β) and stored in an adjacency matrix.
softPower = 20;
adjacency=adjacency(datExpr, power=softPower, type = "signed")

## b. TOM: adjacency matrix is transformed in a Topological Overlapping Matrix (TOM)
TOM = TOMsimilarity(adjacency, TOMType = "signed", TOMDenom = "min");

## c. Save all adjacency and TOM in R.Data
save(sft, softPower, adjacency, TOM, dissTOM, 
     file = "softpowerCutOff0pt7_softPower20_way1.RData")

# 5. Module detection using dynamic tree cutting

## a. Hierarchical clustering is applied to the Topological Overlap Matrix (TOM) dissimilarity to identify network modules.
lnames=load(file="softpowerCutOff0pt7_softPower20_way1.RData");
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");

## b. Module construction
minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8,6)

#pdf(file="GeneDendrogramModuleColors.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#dev.off()

## c. Module eigengenes are calculated to characterize and summarize the expression of each network module.
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");

sizeGrWindow(20, 10)

#pdf(file="ClusteringModuleEigengenes.pdf", width = 20, height = 10)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
#dev.off()

## d. Merging of modules whose expression profiles are very similar (Using module eigengene dissimilarity, build a hierarchical cluster tree. 
Add a cut-off line at the desired height to merge modules with similar expression profiles.)

# Modules-Height Cutoff <- 0
MEDissThres <- 0.0
sizeGrWindow(20, 10)
pdf(file="HeightCutoff0pt0_ClusteringModuleEigengenes.pdf", width = 20, height = 10)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
length(mergedMEs)

#pdf(file="stage_DynamicTreeCut_MergedDynamic_0pt0.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# 6. Module–trait correlation analysis

## a. Associate modules with stages using Pearson correlation
nGenes=ncol(datExpr);
nSamples=nrow(datExpr);
MEs0=moduleEigengenes(datExpr, dynamicColors)$eigengenes
MEs=orderMEs(MEs0)
names(MEs)

moduleTraitCor=cor(MEs, diagData, use="p");
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(12,9)
#pdf(file="HeightCutoff0pt0_38module_Heatmap_ModuleStageAssociations.pdf", width = 12, height = 9)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.lab.x = 0.75,
               cex.lab.y = 0.75,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Stage Associations"))
#dev.off()

## b. output moduleGeneList
modNames=substring(names(MEs),3)

for (i in c(1:length(MEs))){
  moduleGeneList <- names(datExpr)[dynamicColors==modNames[i]]
  write.table(moduleGeneList, quote =FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names =FALSE,
              col.names = FALSE, file = (paste("HeightCutoff0pt0_38module_",
                                               "ME", modNames[i], ".id", sep="")))
}

## c. Associate modules with Anthocyanin phenotypes using Pearson correlation
dim(RoPheno)
names(RoPheno)
RoSamples = rownames(datExpr);
traitRows = match(RoSamples, RoPheno$sample);
RoTraits=RoPheno[traitRows, -1];
rownames(RoTraits)=RoPheno[traitRows, 1];
dim(RoTraits)
collectGarbage();

#=======================================================================
###cor(MEs, RoTraits)
#=======================================================================

nGenes=ncol(datExpr);
nSamples=nrow(datExpr);
MEs0=moduleEigengenes(datExpr, dynamicColors)$eigengenes
MEs=orderMEs(MEs0)
names(MEs)

moduleTraitCor.RoTraits=cor(MEs, RoTraits, use="p");
moduleTraitPvalue.RoTraits=corPvalueStudent(moduleTraitCor.RoTraits, nSamples);
sizeGrWindow(30,9)
pdf(file="HeightCutoff0pt0_38module89phenotype_Heatmap_ModuleAnthocyaninsAssociations.pdf", width = 30, height = 9)
textMatrix =  paste(signif(moduleTraitCor.RoTraits, 2), "\n(",
                    signif(moduleTraitPvalue.RoTraits, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor.RoTraits)
par(mar = c(6, 8.5, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor.RoTraits,
               xLabels = names(RoTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.lab.x = 0.75,
               cex.lab.y = 0.75,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Anthocyanins Associations"))
dev.off()

## d. Identify hub genes and calculate intramodular connectivity for each module.
numberME38 <- table(moduleColors)
hubgeneME38Top <- chooseTopHubInEachModule(datExpr, colorh=moduleColors, omitColors = "grey", power = 10, type = "signed")
hubgeneME38Top
write.table(hubgeneME38Top, quote =FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names =T, 
            col.names = F, file = (paste("HeightCutoff0pt0_ME38_hubgene.txt")))

#=======================================================================
###Connectivity
#=======================================================================

CNT.adj <- intramodularConnectivity(adjacency, colors=moduleColors, scaleByMax = FALSE)
write.table(CNT.adj , quote =FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names =T, 
            col.names = T, file = (paste("HeightCutoff0pt0_ME38_connectivity_adj.txt")))

## e. Associate modules with 7 specific Anthocyanin phenotypes using Pearson correlation 

dim(RoPhenoSpec7)
names(RoPhenoSpec7)
RoSamples = rownames(datExpr);
traitRows = match(RoSamples, RoPhenoSpec7$sample);
RoTraits=RoPhenoSpec7[traitRows, -1];
rownames(RoTraits)=RoPhenoSpec7[traitRows, 1];
dim(RoTraits)
collectGarbage();

#=======================================================================
###cor(MEs, RoTraits)
#=======================================================================

nGenes=ncol(datExpr);
nSamples=nrow(datExpr);
MEs0=moduleEigengenes(datExpr, dynamicColors)$eigengenes
MEs=orderMEs(MEs0)
names(MEs)

moduleTraitCor.RoTraits=cor(MEs, RoTraits, use="p");
moduleTraitPvalue.RoTraits=corPvalueStudent(moduleTraitCor.RoTraits, nSamples);
sizeGrWindow(5,15)
pdf(file="HeightCutoff0pt0_Heatmap_ModuleAnthocyaninsAssociations.pdf", width = 30, height = 9)
textMatrix =  paste(signif(moduleTraitCor.RoTraits, 2), "\n(",
                    signif(moduleTraitPvalue.RoTraits, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor.RoTraits)
par(mar = c(16, 8, 2, 2));

labeledHeatmap(Matrix = moduleTraitCor.RoTraits,
               xLabels = names(RoTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.lab.x = 0.75,
               cex.lab.y = 0.75,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Anthocyanins Associations"))
dev.off()

#7. Module robustness by permutation test

## a. Load data (import RNA-seq data into DESeq2 with module labels and load the TOM data.)
observedTOM <- as.matrix(read.table("TOM_with_gene_names.txt", sep= "\t", header=T))
moduleExp <- read.table("Ro_exprGenes_modLabels.txt", sep= "\t", header=T)

## b. Calculate average Topological overlap (TO) for each observed module
mods <- unique(moduleExp$moduleLabels_final)
modNum <- length(mods)

# separate modules
observedModules <- list()
for (i in 1:modNum){  
  observedModules[[i]] <- subset(moduleExp, moduleExp$moduleLabels_final == i)
}

# get gene numbers for each module
modGeneNums <- c()
for (j in 1: modNum){
  modGeneNums[j] <- nrow(observedModules[[j]])
}

# Calculate average TO
meanTOvec <- c()
for (k in 1: modNum){
  geneids <- c()
  geneids <- rownames(subset(moduleExp, moduleExp$moduleLabels_final == k))
  TO_sub <- as.numeric(observedTOM[geneids, geneids])
  meanTOvec[k] <- mean(TO_sub)
}

save(mods, modNum, observedModules, modGeneNums, geneids, TO_sub, meanTOvec, file = "WGCNA_module_validation_Step1.RData")

## c. Calculate TO for random modules using 100,000 permutations
randomMeanTO <- function(permutationNum, modNum){
  geneids <- randomResult[[permutationNum]][[modNum]]
  TO_sub <- as.numeric(observedTOM[ geneids, geneids ])
  meanTO <- mean(TO_sub)
  return (meanTO)
}

permutationNum <- 100000

genelist_complete <- row.names(moduleExp)

randomSample <- list()
randomResult <- list()

for (x in 1:permutationNum){
  genelist <- genelist_complete
  for (y in 1: modNum){
    cat ("Generating random Module#", y, ", iteration#", x,  "...\n")
    randomSample[[y]] <- sample(genelist, size = modGeneNums[y])
    genelist <- setdiff(genelist, randomSample[[y]])
    randomResult[[x]]<- randomSample
  }
}

mat <- matrix(0, nrow=permutationNum, ncol=modNum)
NumGreaterTO <- rep(0, ncol(mat))
for (ii in 1: permutationNum){
  for (jj in 1: modNum ){
    cat ("Calculating mean TOs for Module#", jj, ", iteration#", ii,  "...\n")
    mat[ii,jj] <- randomMeanTO(permutationNum = ii, modNum = jj)
    if (mat[ii,jj] >= meanTOvec[jj]){
      NumGreaterTO[jj] <- NumGreaterTO[jj]+1
    }
  }
}

colnames(mat) <- paste("Module", 1:modNum, sep= "")
rownames(mat) <- paste("Iteration", 1:permutationNum, sep= "")
meanVecRandom <- apply(mat, 2, mean)

maxMeanTOs <- c()
for (kk in 1:modNum){
  maxMeanTOs[kk] <- max(mat[ ,kk])
}

save(randomMeanTO, permutationNum, genelist_complete, randomResult, NumGreaterTO, mat, maxMeanTOs, file = "WGCNA_module_validation_Step2.RData")

## d. Obtain results from the above steps to calculate permutation-based empirical p-values and save R.Data.
#mean TO for observed modules
print(meanTOvec)

#mean TO random for each module
print(meanVecRandom)

#max TO random for each module
print(maxMeanTOs)

#empirical p-values
print(NumGreaterTO/permutationNum)

save(meanTOvec, meanVecRandom, maxMeanTOs, NumGreaterTO, permutationNum, file = “WGCNA_module_validation_Ro_Result.RData")
