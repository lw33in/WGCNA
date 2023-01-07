# Applying WGCNA to Human Alopecia Areata Skin Biopsy Samples
# Public GEO data: GSE68801 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68801)
# WGCNA: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559

# AA Sample Summary:
# AAP (patchy type disease, n=22 patients; lesional samples n=22, nonlesional samples n=20)
# AAP.T (transient patchy type disease, disease duration less than 1 year, n=6; lesional samples n=6, nonlesional samples n=6)
# AU (alopecia universalis, n=23 patients)
# AT (alopecia totalis, n=9 patients)
# Normal (healthy controls, n=36 patients)

#=====================================================================================
# Directories 
#=====================================================================================
geodir = "/opt/projects/../AA_transcriptomics/GSE68801/expression_data/"
datadir = "/opt/projects/../AA_transcriptomics/GSE68801/WGCNA/Data_WGCNA/"
resultdir = "/opt/projects/../AA_transcriptomics/GSE68801/WGCNA/Result_WGCNA/"

#=====================================================================================
# Prep
#=====================================================================================
install.packages("BiocManager")
BiocManager::install("WGCNA")
library(WGCNA)
library(dplyr)

options(stringsAsFactors = FALSE)

#=====================================================================================
# Load input data files 
#=====================================================================================
setwd(datadir)

expA = read.table('GSE68801_WGCNA_expA.tsv',sep="\t",header=T) 
expN = read.table('GSE68801_WGCNA_expN.tsv',sep="\t",header=T)
clinicsA = read.table('GSE68801_WGCNA_clinicsA.tsv',sep="\t",header=T)
geneanno = read.table('GSE68801_WGCNA_geneanno.tsv',sep="\t",header=T)

# form multi-set expression data
setLabels = c("AA_GSE68801", "Normal_GSE68801")
shortLabels = c("AA", "Normal")
nSets = 2
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(t(expA[-c(1:2)])))
names(multiExpr[[1]]$data) = expA$ProbeName
rownames(multiExpr[[1]]$data) = names(expA)[-c(1:2)]
multiExpr[[2]] = list(data = as.data.frame(t(expN[-c(1:2)])))
names(multiExpr[[2]]$data) = expN$ProbeName
rownames(multiExpr[[2]]$data) = names(expN)[-c(1:2)]
exprSize = checkSets(multiExpr)
clinics = rbind(clinicsA,clinicsN)
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, clinics$name);
  Traits[[set]] = list(data = clinics[traitRows, -1]);
  rownames(Traits[[set]]$data) = clinics[traitRows, 1];
}
collectGarbage()
nGenes = exprSize$nGenes
nSamples = exprSize$nSamples

# save interim data output
setwd(resultdir)
save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize,
     file = "GSE68801_WGCNA_Consensus-dataInput.RData")

#=====================================================================================
# Load data for consensus network analysis 
#=====================================================================================
enableWGCNAThreads()
lnames = load(file = "GSE68801_WGCNA_Consensus-dataInput.RData")

#=====================================================================================
# Choosing the soft-thresholding power: analysis of network topology 
#=====================================================================================
# soft-thresholding powers β
powers = c(seq(4,10,by=1), seq(12,20, by=2))
# initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets)
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]])
collectGarbage()
colors = c("black", "red")
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity")
# minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
    
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
sizeGrWindow(8, 6)
pdf(file = "./GSE68801_WGCNA_scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()

#=====================================================================================
# Network construction and consensus module detection - AA specific
#=====================================================================================
setwd(datadir)

expA = read.table('GSE68801_WGCNA_expA.tsv',sep="\t",header=T) 
clinicsA = read.table('GSE68801_WGCNA_clinicsA.tsv',sep="\t",header=T)
geneanno = read.table('GSE68801_WGCNA_geneanno.tsv',sep="\t",header=T)

datExpr = as.data.frame(t(expA[, -c(1:2)]))
names(datExpr) = expA$ProbeName
rownames(datExpr) = names(expA)[-c(1:2)]

AASamples = rownames(datExpr)
traitRows = match(AASamples, clinicsA$name)

datTraits = clinicsA[traitRows, ]
rownames(datTraits) = clinicsA[traitRows, 1]
collectGarbage()
setwd(resultdir)
save(datExpr, datTraits, file = "GSE68801_WGCNA_AAspecific-Consensus-dataInput.RData")

lnames = load(file = "GSE68801_WGCNA_AAspecific-Consensus-dataInput.RData")
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
pdf(file = "./GSE68801_WGCNA_AAspecific-Scale-free-topology-fit.pdf", wi = 12, he = 6)
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
cor = WGCNA::cor
bwnet = blockwiseModules(datExpr, maxBlockSize = 2000,
                         power = 6, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "AATOM-blockwise",
                         verbose = 3)
# bwnet$colors contains the module assignment
# bwnet$MEs contains the module eigengenes of the modules

# Relabel blockwise modules
bwLabels = matchLabels(bwnet$colors, moduleLabels)
table(bwLabels)
# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)

# Here we show a more flexible way of plotting several trees and colors on one page
sizeGrWindow(12,6)
pdf(file = "./GSE68801_WGCNA_AAspecific-BlockwiseGeneDendrosAndColors.pdf", wi = 12, he = 6);
# Use the layout function for more involved screen sectioning
layout(matrix(c(1:4), 2, 2), heights = c(0.8, 0.2), widths = c(1,1))
#layout.show(4);
nbwBlocks = length(bwnet$dendrograms)
# Plot the dendrogram and the module colors underneath for each block
for (block in 1:nbwBlocks)
  plotDendroAndColors(bwnet$dendrograms[[block]], bwModuleColors[bwnet$blockGenes[[block]]],
                      "Module colors",
                      main = paste("Gene dendrogram and module colors in block", block),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      setLayout = FALSE)
dev.off()

moduleLabels = bwnet$colors
moduleColors = labels2colors(bwnet$colors)
MEs = bwnet$MEs
geneTree = bwnet$dendrograms[[1]]

# save interim data output
setwd(resultdir)
save(MEs, moduleLabels, moduleColors, geneTree,file = "GSE68801_WGCNA_AAspecific-networkConstruction-auto.RData")

#=====================================================================================
#  Network construction and consensus module detection - Automatic, one-step network construction
#=====================================================================================
# call the function blockwiseConsensusModules
net = blockwiseConsensusModules(
  multiExpr, power = 6, minModuleSize = 30, deepSplit = 2,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)
names(net)
# module labels contained in the component colors
# the module eigengenes for each dataset contained in multiMEs
# gene dendrogram (clustering tree) in dendrograms[[1]]
consMEs = net$multiMEs
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]]
sizeGrWindow(8,6)
pdf(file = "./GSE68801_WGCNA_ConsensusDendrogram-auto.pdf", wi = 10, he = 8)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(consMEs, moduleLabels, moduleColors, consTree, file = "GSE68801_WGCNA_Consensus-NetworkConstruction-auto.RData")

#=====================================================================================
#  Network construction and consensus module detection - block-wise network construction
#=====================================================================================
bnet = blockwiseConsensusModules(
  multiExpr, maxBlockSize = 2000, power = 6, minModuleSize = 30,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)

bwLabels = matchLabels(bnet$colors, moduleLabels, pThreshold = 1e-7)
bwColors = labels2colors(bwLabels)

sizeGrWindow(12,6)
pdf(file = "./GSE68801_WGCNA_BlockwiseGeneDendrosAndColors-6.pdf", wi = 12, he = 6)
layout(matrix(c(1:4), 2, 2), heights = c(0.8, 0.2), widths = c(1,1))
nBlocks = length(bnet$dendrograms)
for (block in 1:nBlocks)
  plotDendroAndColors(bnet$dendrograms[[block]], moduleColors[bnet$blockGenes[[block]]],
                      "Module colors",
                      main = paste("Gene dendrogram and module colors in block", block),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      setLayout = FALSE)
dev.off()

# save interim data output
setwd(resultdir)
save(consMEs, moduleLabels, moduleColors, consTree, file = "GSE68801_WGCNA_Consensus-NetworkConstruction-auto.RData")

#=====================================================================================
# setting up the R session for Relating consensus modules to AA set-specific modules
#=====================================================================================
lnames = load(file = "GSE68801_WGCNA_Consensus-NetworkConstruction-auto.RData")
# "consMEs"      "moduleLabels" "moduleColors" "consTree"  
lnames

# We have loaded the variables multiExpr and Traits containing the expression and trait data, respectively. 
# Further, expression data dimensions are stored in nGenes and nSamples

#=====================================================================================
# Relating consensus modules to AA set-specific modules
#=====================================================================================
lnames = load("GSE68801_WGCNA_AAspecific-networkConstruction-auto.RData")
#lnames
# Rename variables to avoid conflicts
AALabels = moduleLabels;
AAColors = moduleColors;
AATree = geneTree;
AAMEs = orderMEs(MEs, greyName = "ME0");
# Next we load the results of the consensus module identification:
lnames = load("GSE68801_WGCNA_Consensus-NetworkConstruction-auto.RData")
#lnames 
# "consMEs"      "moduleLabels" "moduleColors" "consTree"  

# Isolate the module labels in the order they appear in ordered module eigengenes
AAModuleLabels = substring(names(AAMEs), 3)
consModuleLabels = substring(names(consMEs[[1]]$data), 3)
# Convert the numeric module labels to color labels
AAModules = labels2colors(as.numeric(AAModuleLabels))
consModules = labels2colors(as.numeric(consModuleLabels))
# Numbers of female and consensus modules
nAAMods = length(AAModules)
nConsMods = length(consModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nAAMods, ncol = nConsMods);
CountTbl = matrix(0, nrow = nAAMods, ncol = nConsMods);
# Execute all pairwaise comparisons
for (AAmod in 1:nAAMods)
  for (cmod in 1:nConsMods)
  {
    AAMembers = (AAColors == AAModules[AAmod]);
    consMembers = (moduleColors == consModules[cmod]);
    pTable[AAmod, cmod] = -log10(fisher.test(AAMembers, consMembers, alternative = "greater")$p.value);
    CountTbl[AAmod, cmod] = sum(AAColors == AAModules[AAmod] & moduleColors ==
                                  consModules[cmod])
  }
# To display the p-value and count tables in an informative way
# create a color-coded table of the intersection counts. The colors will indicate the p-value significance
# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
AAModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)
# Actual plotting
sizeGrWindow(10,7 );
pdf(file = "./GSE68801_WGCNA_ConsensusVsAAModules.pdf", wi = 30, he = 15)
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(8, 10.4, 2.7, 1)+0.3)
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", AAModules),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("AA ", AAModules, ": ", AAModTotals, sep=""),
               textMatrix = CountTbl,
               colors = blueWhiteRed(100)[50:100],
               main = "Correspondence of AA set-specific and AA-NORMAL consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)
dev.off()
