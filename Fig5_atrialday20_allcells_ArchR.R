library(ArchR)
library(parallel)
##setting threads = 1
addArchRThreads(threads = 1)
##setting genome to hg38
addArchRGenome("hg38")

##CREATING ARROW FILES -----------------------------------------------------------------------------------------
inputDir <- "/path/to/arrowfiles"
smplDir <- list.dirs(inputDir, recursive = FALSE)
inputFiles <- file.path(smplDir, "/fragments.tsv.gz") 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = basename(smplDir),
  minTSS = 4, #May want to increase later
  minFrags = 1000, #using default
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

##DOUBLET IDENTIFICATION -----------------------------------------------------------------------------------------
##doublet identification - wt1 = 0.98773; het1 = 0.98698; hom1 = 0.9794, wt3 = 0.96946, het3 = 0.96796, hom3 = 0.96778; 
#ideally these values should be > 0.9, we do doublet filtering later
doubScores <- addDoubletScores(input = ArrowFiles,
                               k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                               knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
                               LSIMethod = 1)

##CREATE AND SAVE ARCHR PROJECT -----------------------------------------------------------------------------------------
atrialatac <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "/path/to/outputdir",
  copyArrows = TRUE) #This is recommended so that if you modify the Arrow files you have an original copy for later usage
##save project
saveArchRProject(ArchRProj = atrialatac, outputDirectory = "/path/to/outputdir", load = TRUE)

##SAMPLE QC STATISTICS -----------------------------------------------------------------------------------------
##Make a ridge plot for each sample for the TSS enrichment scores
p1 <- plotGroups(ArchRProj = atrialatac, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment",plotAs = "ridges")
p1
##Make a violin plot for each sample for the TSS enrichment scores.
p2 <- plotGroups(ArchRProj = atrialatac, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
p2
##Make a ridge plot for each sample for the log10(unique nuclear fragments).
p3 <- plotGroups(ArchRProj = atrialatac, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)",plotAs = "ridges")
p3
##plot fragment sizes
p4 <- plotFragmentSizes(ArchRProj = atrialatac)
p4
##TSS enrichment profiles should show a clear peak in the center and a smaller shoulder peak right-of-center 
##which is caused by the well-positioned +1 nucleosome.
addArchRThreads(threads = 2)
p5 <- plotTSSEnrichment(ArchRProj = atrialatac)
p5
#save editable pdf
plotPDF(p1,p2,p3,p4,p5, name = "QC-Sample-Statistics.pdf", ArchRProj = atrialatac, addDOC = FALSE, width = 5, height = 5)
#set back to threads = 1
addArchRThreads(threads = 1)


##DOUBLET FILTERING-----------------------------------------------------------------------------------------
##first look at how many cells before filtering = 94,485
atrialatac
# ##remove doublets
atrialatac <- filterDoublets(atrialatac)
# Hom-1 : 1482 of 22889 (6.5%)
# Het-3 : 1431 of 11966 (12%)
# Hom-3 : 2269 of 15064 (15.1%)
# Het-1 : 1919 of 13853 (13.9%)
# WT-1 : 2078 of 14418 (14.4%)
# WT-3 : 2655 of 16295 (16.3%)
##look at how many cells remain
atrialatac
##82,651 cells remaining

##ADDITIONAL QC - FILTERED MORE CELLS BY UNIQUE FRAGS AND TSS ENRICHMENT -----------------------------------------------------------------------------------------
##first look at how many cells before filtering = 82,651
#Get nFrags and TSSEnrichment scores
df <- getCellColData(atrialatac, select = c("log10(nFrags)", "TSSEnrichment"))
df
#Build plot
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 1.0)+0.5),
  ylim = c(0, quantile(df[,2], probs = 1.0)+5)
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
#View plot
p
#Save plot
ggsave("atrialatac_TSSbyLog10nFrags_initial.pdf", plot = p, device = "pdf", path = "/path/to/outputdir", width = 5, height = 5, dpi = 300)
#Before filtering, visualize subset of cells with new cutoffs
#TSS Enrichiment >= 7, TSS Enrichiment <= 30 nFrags >= log10 (2000) = 3.3
idxPass <- which(atrialatac$TSSEnrichment >= 7 & atrialatac$nFrags >=2000)
cellsPass <- atrialatac$cellNames[idxPass]
df3 <- getCellColData(atrialatac[cellsPass, ], select = c("log10(nFrags)", "TSSEnrichment"))
p3 <- ggPoint(
  x = df3[,1], 
  y = df3[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df3[,1], probs = 1.0)+0.5),
  ylim = c(0, quantile(df3[,2], probs = 1.0)+5)
) + geom_hline(yintercept = 7, lty = "dashed") + geom_hline(yintercept = 30, lty = "dashed") + geom_vline(xintercept = 3.3, lty = "dashed")
#View plot
ggAlignPlots(p,p3, type="h")
#Save plot
ggsave("atrialatac_TSSbyLog10nFrags_filtered.pdf", plot = p3, device = "pdf", path = "/path/to/outputdir", width = 5, height = 5, dpi = 300)
#Looks good, now filter these cells from project 
atrialatac <- atrialatac[cellsPass, ]
atrialatac
##remaining cells = 38,420

##DIMENSIONALITY REDUCTION ---------------------------------------------------------------------------------------------------------------------------------------
##dimensionality reduction; can run multiple resolutions if length of # of iterations = clusterParams + 1
atrialatac <- addIterativeLSI(ArchRProj = atrialatac, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 11,force = TRUE,seed = 1,
                              clusterParams = list(resolution = c(0.1,0.2,0.4,0.5,0.6,0.9,1.0,1.2,1.4,1.5), sampleCells = 10000, n.start = 10), varFeatures = 25000, dimsToUse = 1:30)

##Add columns for batch and genotype ------------------------------------------------------------------------------------------------------------------------------
bioNames <- gsub(".*-","",atrialatac$Sample)
head(bioNames)
tail(bioNames)
atrialatac$batch <- bioNames

bioNames <- gsub("-.*","",atrialatac$Sample)
head(bioNames)
tail(bioNames)
atrialatac$Genotype <- bioNames

##BATCH Correction -----------------------------------------------------------------------------------------------------------------------------------------------
atrialatac <- addHarmony(ArchRProj = atrialatac,reducedDims = "IterativeLSI",name = "Harmony",groupBy = "batch", force = T)

##CLUSTERING AND UMAPS -------------------------------------------------------------------------------------------------------------------------------------------
atrialatac <- addClusters(input = atrialatac, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res0.1",resolution = 0.1, seed = 1, force = T)
atrialatac <- addClusters(input = atrialatac, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res0.2",resolution = 0.2, seed = 1, force = T)
atrialatac <- addClusters(input = atrialatac, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res0.4",resolution = 0.4, seed = 1, force = T)
atrialatac <- addClusters(input = atrialatac, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res0.5",resolution = 0.5, seed = 1, force = T)
atrialatac <- addClusters(input = atrialatac, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res0.6",resolution = 0.6, seed = 1, force = T)
atrialatac <- addClusters(input = atrialatac, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res0.9",resolution = 0.9, seed = 1, force = T)
atrialatac <- addClusters(input = atrialatac, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res1.0",resolution = 1.0, seed = 1, force = T)
atrialatac <- addUMAP(ArchRProj = atrialatac, reducedDims = "Harmony", name = "UMAP", nNeighbors = 20, minDist = 0.2, metric = "cosine", seed = 1, force=T)
p6 <- plotEmbedding(ArchRProj = atrialatac, colorBy = "cellColData", name = "Clusters_res0.4", embedding = "UMAP")
pb <- plotEmbedding(ArchRProj = atrialatac, colorBy = "cellColData", name = "batch", embedding = "UMAP")
test <- c("Het" = "red","Hom" = "steelblue3","WT" = "black")
p7 <- plotEmbedding(ArchRProj = atrialatac, colorBy = "cellColData", name = "Genotype", embedding = "UMAP", pal = test)
test2 <- c("Het-1" = "red","Het-3" = "magenta","Hom-1" = "steelblue3","Hom-3" = "blue4","WT-1" = "black","WT-3"="green4")
p8 <- plotEmbedding(ArchRProj = atrialatac, colorBy = "cellColData", name = "Sample", embedding = "UMAP", pal = test2)
ggAlignPlots(p6, p7, type = "h")
plotPDF(p6,p7,p8,pb, name = "Plot-UMAP-Genotype-Batch-Clusters-res0.4.pdf", ArchRProj = atrialatac, addDOC = FALSE, width = 5, height = 5)


##MAGIC IMPUTATION ----------------------------------------------------------------------------------------------------------------------------------
##impute gene scores with MAGIC to smooth signal over nearby cells and improve sparse data
atrialatac <- addImputeWeights(atrialatac)

##FEATUREPLOTs -------------------------------------------------------------------------------------------------------
##single gene feature plot
p3 <- plotEmbedding(ArchRProj = atrialatac, colorBy = "GeneScoreMatrix", name = "COL3A1", embedding = "UMAP", rastr = FALSE,plotAs = "points",size = 0.5,
                    quantCut = c(0.01, 0.95),imputeWeights = getImputeWeights(atrialatac))
p3
##multiple feature plots
celltypegenes <- c("TNNT2","COL1A1","EPCAM","NEFM","SERPINA1")
p9 <- plotEmbedding(ArchRProj = atrialatac, colorBy = "GeneScoreMatrix", name = celltypegenes, embedding = "UMAP",rastr = FALSE,plotAs = "points",
                    quantCut = c(0.01, 0.95),imputeWeights = getImputeWeights(atrialatac))
##to plot a single gene from above list
p9$TNNT2
##to plot all
p10 <- lapply(p9, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())})
do.call(cowplot::plot_grid, c(list(ncol = 3),p10))
plotPDF(plotList = p10, name = "Featureplots_celltypes.pdf", ArchRProj = atrialatac, addDOC = FALSE, width = 5, height = 5)

##SAVE ARCHR PROJECT -----------------------------------------------------------------------------------------
saveArchRProject(ArchRProj = atrialatac, outputDirectory = "/path/to/outputdir", load = TRUE)





