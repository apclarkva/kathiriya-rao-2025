library(ArchR)
library(parallel)
##setting threads = 1
addArchRThreads(threads = 1)
##setting genome to hg38
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

##subset cm clusters from res 0.4---------------------------------------------------------------------------------------------------------------------------------------
#load previously made ArchR object containing all cells
atrialatac <- loadArchRProject(path="/path/to/ArchRproject/")
c2 <- BiocGenerics::which(atrialatac$Clusters_res0.4 %in% "C2")
c3 <- BiocGenerics::which(atrialatac$Clusters_res0.4 %in% "C3")
c4 <- BiocGenerics::which(atrialatac$Clusters_res0.4 %in% "C4")
c6 <- BiocGenerics::which(atrialatac$Clusters_res0.4 %in% "C6")
cellsPass2 <- atrialatac$cellNames[c2]
cellsPass3 <- atrialatac$cellNames[c3]
cellsPass4 <- atrialatac$cellNames[c4]
cellsPass6 <- atrialatac$cellNames[c6]
cellsPass <- c(cellsPass2, cellsPass3, cellsPass4, cellsPass6)
atrialatac[cellsPass, ]
atrialatacms <- subsetArchRProject(ArchRProj = atrialatac, cells = cellsPass, outputDirectory = "ArchRSubset",
                                   dropCells = TRUE, logFile = NULL, threads = getArchRThreads(),force = TRUE)
#load new subsetted object
atrialatacms <- loadArchRProject(path="/path/to/ArchRSubset")

##DIMENSIONALITY REDUCTION ---------------------------------------------------------------------------------------------------------------------------------------
##dimensionality reduction; can run multiple resolutions if length of # of iterations = clusterParams + 1
atrialatacms<- addIterativeLSI(ArchRProj = atrialatacms, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 8,force = TRUE,seed = 1,
                               clusterParams = list(resolution = c(0.1,0.2,0.4,0.5,0.6,0.9,1.0,1.1,1.2,1.4,1.5), sampleCells = 10000, n.start = 10), varFeatures = 25000, dimsToUse = 1:30)


##Add columns for batch and genotype ------------------------------------------------------------------------------------------------------------------------------
bioNames <- gsub(".*-","",atrialatacms$Sample)
head(bioNames)
tail(bioNames)
atrialatacms$batch <- bioNames

bioNames <- gsub("-.*","",atrialatacms$Sample)
head(bioNames)
tail(bioNames)
atrialatacms$Genotype <- bioNames

##BATCH Correction -----------------------------------------------------------------------------------------------------------------------------------------------
atrialatacms<- addHarmony(ArchRProj = atrialatacms,reducedDims = "IterativeLSI",name = "Harmony",groupBy = "batch", force = T)

##CLUSTERING AND UMAPS -------------------------------------------------------------------------------------------------------------------------------------------
atrialatacms<- addClusters(input = atrialatacms, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res0.1",resolution = 0.1, seed = 1, force = T)
atrialatacms<- addClusters(input = atrialatacms, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res0.2",resolution = 0.2, seed = 1, force = T)
atrialatacms<- addClusters(input = atrialatacms, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res0.4",resolution = 0.4, seed = 1, force = T)
atrialatacms<- addClusters(input = atrialatacms, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res0.5",resolution = 0.5, seed = 1, force = T)
atrialatacms<- addClusters(input = atrialatacms, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res0.6",resolution = 0.6, seed = 1, force = T)
atrialatacms<- addClusters(input = atrialatacms, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res0.9",resolution = 0.9, seed = 1, force = T)
atrialatacms<- addClusters(input = atrialatacms, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res1.0",resolution = 1.0, seed = 1, force = T)
atrialatacms<- addClusters(input = atrialatacms, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res1.1",resolution = 1.1, seed = 1, force = T)
atrialatacms<- addClusters(input = atrialatacms, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res1.2",resolution = 1.2, seed = 1, force = T)
atrialatacms<- addClusters(input = atrialatacms, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res1.4",resolution = 1.4, seed = 1, force = T)
atrialatacms<- addClusters(input = atrialatacms, reducedDims = "Harmony", method = "Seurat",name = "Clusters_res1.5",resolution = 1.5, seed = 1, force = T)
#add UMAP
atrialatacms <- addUMAP(ArchRProj = atrialatacms, reducedDims = "Harmony", name = "UMAP", nNeighbors = 20, minDist = 0.2, metric = "cosine", seed = 1, force=T)
p6 <- plotEmbedding(ArchRProj = atrialatacms, colorBy = "cellColData", name = "Clusters_res0.4", embedding = "UMAP")
pb <- plotEmbedding(ArchRProj = atrialatacms, colorBy = "cellColData", name = "batch", embedding = "UMAP")
test <- c("Het" = "red","Hom" = "steelblue3","WT" = "black")
p7 <- plotEmbedding(ArchRProj = atrialatacms, colorBy = "cellColData", name = "Genotype", embedding = "UMAP", pal = test)
test2 <- c("Het-1" = "red","Het-3" = "magenta","Hom-1" = "steelblue3","Hom-3" = "blue4","WT-1" = "black","WT-3"="green4")
p8 <- plotEmbedding(ArchRProj = atrialatacms, colorBy = "cellColData", name = "Sample", embedding = "UMAP", pal = test2)
ggAlignPlots(p6, p7, type = "h")
plotPDF(p6,p7,p8,pb, name = "Plot-UMAP-Genotype-Batch-Clusters-res1.4.pdf", ArchRProj = atrialatacms, addDOC = FALSE, width = 5, height = 5)


##MAGIC IMPUTATION----------------------------------------------------------------------------------------------------------------------------------
##impute gene scores with MAGIC to smooth signal over nearby cells and improve sparse data
atrialatacms<- addImputeWeights(atrialatacms)

##FEATUREPLOTS -------------------------------------------------------------------------------------------------------
##single gene feature plot
p3 <- plotEmbedding(ArchRProj = atrialatacms, colorBy = "GeneScoreMatrix", name = "DES", embedding = "UMAP", rastr = FALSE,plotAs = "points",size = 0.5,
                    quantCut = c(0.01, 0.95),imputeWeights = getImputeWeights(atrialatacms))
p3
##multiple feature plots
markerGenes  <- c("MYH6","PITX2","HAMP","NPPA","COL3A1","HBD")
p9 <- plotEmbedding(ArchRProj = atrialatacms, colorBy = "GeneScoreMatrix", name = markerGenes, embedding = "UMAP",rastr = FALSE,plotAs = "points",
                    quantCut = c(0.01, 0.95),imputeWeights = getImputeWeights(atrialatacms))
##to plot a single gene from above list
p9$PITX2
##to plot all
p10 <- lapply(p9, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())})
do.call(cowplot::plot_grid, c(list(ncol = 3),p10))
plotPDF(plotList = p10, name = "Featureplots_cmsubset.pdf", ArchRProj = atrialatacms, addDOC = FALSE, width = 5, height = 5)


##SAVE PROJECT-------------------------------------------------------------------------------------------------------
saveArchRProject(ArchRProj = atrialatacms, outputDirectory = "~/Desktop/atrial-deep/nonorm_sct/ArchRSubset", load = TRUE)


##CELL NUMBER TABLES-------------------------------------------------------------------------------------------------------
#cellnumbers per genotype in clusters
atrialatacms_cellnumbers_res1.2 <- table(atrialatacms$Sample, atrialatacms$Clusters_res1.2)
write.csv(atrialatacms_cellnumbers_res1.2, file = "atrialatacms_cellnumbers_res1.2.csv")
#BY GENOTYPE
table(atrialatacms$Genotype, atrialatacms$Clusters_res0.4)
#       C1   C2   C3   C4   C5   C6
# Het  247 1701 2209   47  239 2858
# Hom    2   11   18 3828 2312  135
# WT   576 3797 4812    5  139  750
#BY SAMPLE
table(atrialatacms$Sample, atrialatacms$Clusters_res0.4)
#         C1   C2   C3   C4   C5   C6
# Het-1  117 1431 1464   44   55 1480
# Het-3  130  270  745    3  184 1378
# Hom-1    0    0    0 3797   83    9
# Hom-3    2   11   18   31 2229  126
# WT-1   347 2939 2717    3   47  280
# WT-3   229  858 2095    2   92  470
#BY GENOTYPE
table(atrialatacms$Genotype, atrialatacms$Clusters_res1.4)
#       C1  C10  C11  C12  C13  C14   C2   C3   C4   C5   C6   C7   C8   C9
# Het    4 1455 1533  307  984  740   21   19  871 1111  187   46    2   21
# Hom 1043   21   22    2    3    0 1476 1314  110   25 1675  612    0    3
# WT     0 1136  943  625 3895 2944    0    6  204  134  111   26    7   48


##PSEUDOBULK REPLICATES ----------------------------------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg38)
##ArchR makes pseudobulk replicates in a sample aware fashion
##more cells in the group than minCells * minRep
##at least 1 sample has minCells total
##at least minRep samples have minCells total
##at least maxRep samples have minCells total
atrialatacms<- addGroupCoverages(ArchRProj = atrialatacms, groupBy = "ClusterByGenotype",useLabels = TRUE,
                                 minCells = 500, maxCells = 1000,
                                 maxFragments = 25 * 10^6, minReplicates = 2, maxReplicates = 2,
                                 sampleRatio = 0.8,kmerLength = 6,threads = getArchRThreads(), returnGroups = FALSE,
                                 parallelParam = NULL, force = TRUE,verbose = TRUE,logFile = createLogFile("addGroupCoverages"))


##CALL PEAKS ---------------------------------------------------------------------------------------------------------
##Note: with reproducibility = 2, 2 OUT OF 2 pseodobulk reps for each group will be required to contain a peak call at the locus
pathToMacs2 <- "/Library/Frameworks/Python.framework/Versions/3.9/bin/MACS2"
atrialatacms<- addReproduciblePeakSet(ArchRProj = atrialatacms, groupBy = "ClusterByGenotype",
                                      reproducibility = "2",peaksPerCell = 500,maxPeaks = 150000, minCells = 20,
                                      excludeChr = c("chrM","chrY","chrX"),pathToMacs2 = pathToMacs2,shift = -75,extsize = 150,cutOff = 0.1,extendSummits = 250,
                                      promoterRegion = c(2000, 100),plot = TRUE,force = TRUE,verbose = TRUE,
                                      logFile = createLogFile("addReproduciblePeakSet"))
# Group    nCells nCellsUsed nReplicates nMin nMax maxPeaks
#C1_x_Het    247        242           2  206  209   121000
#C1_x_Hom      2          2           2    2    2     1000
#C1_x_WT    576        490            2  331  339   150000
#C2_x_Het   1701       1500           2  500 1000   150000
#C2_x_Hom     11         11           2   11   11     5500
#C2_x_WT   3797       1858           2   858 1000   150000
#C3_x_Het   2209       1745           2  745 1000   150000
#C3_x_Hom     18         18           2   18   18     9000
#C3_x_WT   4812       2000           2  1000 1000   150000
#C4_x_Het     47         47           2   47   47    23500
#C4_x_Hom   3828       1500           2  500 1000   150000
#C4_x_WT      5          5           2    5    5     2500
#C5_x_Het    239        238           2  201  207   119000
#C5_x_Hom   2312       1500           2  500 1000   150000
#C5_x_WT    139        139           2   131  133    69500
#C6_x_Het   2858       2000           2  1000 1000   150000
#C6_x_Hom    135        135           2  128  128    67500
#C6_x_WT    750        680           2   500  500   150000
##addPeakMatrix
atrialatacms<- addPeakMatrix(atrialatacms)
getAvailableMatrices(atrialatacms)

##PEAK CALL SUMMARY ---------------------------------------------------------------------------------------------------
peakDF <- metadata(atrialatacms@peakSet)$PeakCallSummary
pal <- c("Distal" = "#60BA64", "Exonic" = "#73C6FF", "Intronic" = "#620FA3", "Promoter" = "#FFC554") 
pal <- pal[!duplicated(names(pal))]

lengthMax <- split(peakDF$Freq, peakDF$Group) %>% lapply(sum) %>% unlist %>% max

colnames(peakDF)[colnames(peakDF)=="Var1"] <- "PeakAnno"

ggplot(peakDF, aes(x=Group, y=Freq, fill=PeakAnno)) + 
  geom_bar(stat = "identity") + 
  theme_ArchR(xText90 = TRUE) +
  ylab("Number of Peaks (x10^3)") +
  xlab("") + 
  theme(legend.position = "bottom", 
        legend.key = element_rect(size = 2), 
        legend.box.background = element_rect(color = NA)
  ) +
  scale_fill_manual(values=pal) +
  scale_y_continuous(
    breaks = seq(0, lengthMax * 2,50), 
    limits = c(0, lengthMax * 1.1), 
    expand = c(0,0)
  )


##DIFFERENTIAL PEAK ANALYSIS - het vs wt-------------------------------------------------------------------------------
markerTest_C6HetvsC2C3WT <- getMarkerFeatures(
  ArchRProj = atrialatacms,
  useMatrix = "PeakMatrix",
  groupBy = "ClusterByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C6_x_Het",
  bgdGroups = c("C2_x_WT","C3_x_WT"))

##sort the list by logfc and save as csv
markerTest_C6HetvsC2C3WT_upinHet <- getMarkers(markerTest_C6HetvsC2C3WT, cutOff = "FDR <= 0.1 & Log2FC >= 0.25", returnGR = TRUE)
#write.csv(markerTest_C6HetvsC2C3WT_upinHet, file = "markerTest_C6HetvsC2C3WT_upinHet.csv")
markerTest_C6HetvsC2C3WT_downinHet <- getMarkers(markerTest_C6HetvsC2C3WT, cutOff = "FDR <= 0.1 & Log2FC <= -0.25", returnGR = TRUE)
#write.csv(markerTest_C6HetvsC2C3WT_downinHet, file = "markerTest_C6HetvsC2C3WT_downinHet.csv")
# #Let’s convert the up list to a GRanges object, sorted by Log2FC
markerTest_C6HetvsC2C3WT_upinHet_gr <- markerTest_C6HetvsC2C3WT_upinHet$`C6_x_Het`
order_upinHet <- order(markerTest_C6HetvsC2C3WT_upinHet_gr$Log2FC, start(markerTest_C6HetvsC2C3WT_upinHet_gr))
sorted_markerTest_C6HetvsC2C3WT_upinHet_gr <- markerTest_C6HetvsC2C3WT_upinHet_gr[order_upinHet]
# #Let’s convert the down list to a GRanges object, sorted by Log2FC
markerTest_C6HetvsC2C3WT_downinHet_gr <- markerTest_C6HetvsC2C3WT_downinHet$`C6_x_Het`
order_downinHet <- order(markerTest_C6HetvsC2C3WT_downinHet_gr$Log2FC, start(markerTest_C6HetvsC2C3WT_downinHet_gr))
sorted_markerTest_C6HetvsC2C3WT_downinHet_gr <- markerTest_C6HetvsC2C3WT_downinHet_gr[order_downinHet]


##DIFFERENTIAL PEAK ANALYSIS - hom vs wt ---------------------------------------------------------------------------
markerTest_C5HomvsC2C3WT <- getMarkerFeatures(
  ArchRProj = atrialatacms,
  useMatrix = "PeakMatrix",
  groupBy = "ClusterByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C5_x_Hom",
  bgdGroups = c("C3_x_WT","C2_x_WT"))

##sort the list by logfc and save as csv
markerTest_C5HomvsC2C3WT_upinHom <- getMarkers(markerTest_C5HomvsC2C3WT, cutOff = "FDR <= 0.1 & Log2FC >= 0.25", returnGR = TRUE)
#write.csv(markerTest_C5HomvsC2C3WT_upinHom, file = "markerTest_C5HomvsC2C3WT_upinHom_fdr0.1.csv")
markerTest_C5HomvsC2C3WT_downinHom <- getMarkers(markerTest_C5HomvsC2C3WT, cutOff = "FDR <= 0.1 & Log2FC <= -0.25", returnGR = TRUE)
#write.csv(markerTest_C5HomvsC2C3WT_downinHom, file = "markerTest_C5HomvsC2C3WT_downinHom_fdr0.1.csv")
# #Let’s convert the up list to a GRanges object, sorted by Log2FC
markerTest_C5HomvsC2C3WT_upinHom_gr <- markerTest_C5HomvsC2C3WT_upinHom$`C5_x_Hom`
order_upinHom <- order(markerTest_C5HomvsC2C3WT_upinHom_gr$Log2FC, start(markerTest_C5HomvsC2C3WT_upinHom_gr))
sorted_markerTest_C5HomvsC2C3WT_upinHom_gr <- markerTest_C5HomvsC2C3WT_upinHom_gr[order_upinHom]
# #Let’s convert the up list to a GRanges object, sorted by Log2FC
markerTest_C5HomvsC2C3WT_downinHom_gr <- markerTest_C5HomvsC2C3WT_downinHom$`C5_x_Hom`
order_downinHom <- order(markerTest_C5HomvsC2C3WT_downinHom_gr$Log2FC, start(markerTest_C5HomvsC2C3WT_downinHom_gr))
sorted_markerTest_C5HomvsC2C3WT_downinHom_gr <- markerTest_C5HomvsC2C3WT_downinHom_gr[order_downinHom]



### check overlaps between het and hom peaks -------------------------------------------------------------------------------
# Find overlaps
overlaps_down_hethom <- findOverlaps(markerTest_C6HetvsC2C3WT_downinHet_gr, markerTest_C5HomvsC2C3WT_downinHom_gr, type = "any", ignore.strand = F)
overlaps_up_hethom <- findOverlaps(markerTest_C6HetvsC2C3WT_upinHet_gr, markerTest_C5HomvsC2C3WT_upinHom_gr, type = "any", ignore.strand = F)
# Extract matching data
results_down_hethom <- data.frame(chr_query = seqnames(markerTest_C6HetvsC2C3WT_downinHet_gr)[queryHits(overlaps_down_hethom)],
                                  start_query = start(markerTest_C6HetvsC2C3WT_downinHet_gr)[queryHits(overlaps_down_hethom)],
                                  end_query = end(markerTest_C6HetvsC2C3WT_downinHet_gr)[queryHits(overlaps_down_hethom)],
                                  chr_subject = seqnames(markerTest_C5HomvsC2C3WT_downinHom_gr)[subjectHits(overlaps_down_hethom)],
                                  start_subject = start(markerTest_C5HomvsC2C3WT_downinHom_gr)[subjectHits(overlaps_down_hethom)],
                                  end_subject = end(markerTest_C5HomvsC2C3WT_downinHom_gr)[subjectHits(overlaps_down_hethom)])
print(results_down_hethom)
results_up_hethom <- data.frame(chr_query = seqnames(markerTest_C6HetvsC2C3WT_upinHet_gr)[queryHits(overlaps_up_hethom)],
                                start_query = start(markerTest_C6HetvsC2C3WT_upinHet_gr)[queryHits(overlaps_up_hethom)],
                                end_query = end(markerTest_C6HetvsC2C3WT_upinHet_gr)[queryHits(overlaps_up_hethom)],
                                chr_subject = seqnames(markerTest_C5HomvsC2C3WT_upinHom_gr)[subjectHits(overlaps_up_hethom)],
                                start_subject = start(markerTest_C5HomvsC2C3WT_upinHom_gr)[subjectHits(overlaps_up_hethom)],
                                end_subject = end(markerTest_C5HomvsC2C3WT_upinHom_gr)[subjectHits(overlaps_up_hethom)])
print(results_up_hethom)
#save
write.csv(results_down_hethom, file = "commonpeaks_down_C6hetC5hom.csv")
write.csv(results_up_hethom, file = "commonpeaks_up_C6hetC5hom.csv")



##VOLCANO PLOT ---------------------------------------------------------------------------------------------------------------------------------------------------
volcano_markerTest_C6HetvsC2C3WT <- plotMarkers(seMarker = markerTest_C6HetvsC2C3WT, name = "C6_x_Het",
                                                cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25", plotAs = "Volcano") + scale_color_manual(values = c("red","grey","black"))
volcano_markerTest_C6HetvsC2C3WT
plotPDF(volcano_markerTest_C6HetvsC2C3WT, name = "newcolorsvolcano_markerTest_C6HetvsC2C3WT_fdr>0.1_logfc>0.25", width = 5, height = 5, ArchRProj = atrialatacms, addDOC = FALSE)

volcano_markerTest_C5HomvsC2C3WT <- plotMarkers(seMarker = markerTest_C5HomvsC2C3WT, name = "C5_x_Hom",
                                                cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25", plotAs = "Volcano") + scale_color_manual(values = c("steelblue","grey","black"))
volcano_markerTest_C5HomvsC2C3WT
plotPDF(volcano_markerTest_C5HomvsC2C3WT, name = "newcolorsvolcano_markerTest_C5HomvsC2C3WT_fdr>0.1_logfc>0.25", width = 5, height = 5, ArchRProj = atrialatacms, addDOC = FALSE)

##plot volcano using ggplot2 to change colors


##BROWSER TRACK ------------------------------------------------------------------------------------------------------------
##differential peaks - plot with ChIP track
library(readxl)
library(dplyr)
markerTest_C6HetvsC2C3WT_sig <- getMarkers(markerTest_C6HetvsC2C3WT, cutOff = "FDR <= 0.1 & abs(Log2FC) >= -0.25", returnGR = TRUE)
markerTest_C6HetvsC2C3WT_gr <- markerTest_C6HetvsC2C3WT_sig$`C6_x_Het`
markerTest_C5HomvsC2C3WT_sig <- getMarkers(markerTest_C5HomvsC2C3WT, cutOff = "FDR <= 0.1 & abs(Log2FC) >= -0.25", returnGR = TRUE)
markerTest_C5HomvsC2C3WT_gr <- markerTest_C5HomvsC2C3WT_sig$`C5_x_Hom`
chip_peaks <- read_excel("/path/to/TBX5-ChIP_d11_d20_d45_peaks.xlsx")
# Creating GRanges objects from data frames
gr_chip <- GRanges(seqnames = chip_peaks$chr,ranges = IRanges(start = chip_peaks$start, end = chip_peaks$end))
p <- plotBrowserTrack(ArchRProj = atrialatacms, groupBy = "ClusterByGenotype", 
                      geneSymbol = c("NR2F1","NPPA"),
                      features =  GRangesList(TrackA = markerTest_C5HomvsC2C3WT_gr,TrackB = gr_chip),
                      useGroups = c("C2_x_WT","C3_x_WT","C6_x_Het","C5_x_Hom"))
grid::grid.newpage()
grid::grid.draw(p$NR2F1)
grid::grid.draw(p$NPPA)


##Export Bigwig file for each group --------------------------------------------------------------------------------------------------------------
getGroupBW(ArchRProj = atrialatacms, groupBy = "ClusterByGenotype", normMethod = "ReadsInTSS",verbose = TRUE)


##MOTIF ENRICHMENT --------------------------------------------------------------------------------------------------------------
atrialatacms<- addMotifAnnotations(ArchRProj = atrialatacms, motifSet = "cisbp", name = "Motif", force=T)
motifsUp <- peakAnnoEnrichment(seMarker = markerTest_C5HomvsC2C3WT, ArchRProj = atrialatacms, peakAnnotation = "Motif",
                               cutOff = "FDR <= 0.05 & Log2FC >= 0.25")
df1 <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df1 <- df1[order(df1$mlog10Padj, decreasing = TRUE),]
df1$rank <- seq_len(nrow(df1))
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + geom_point(size = 1) + ggrepel::geom_label_repel(data = df[rev(seq_len(30)), ],
                                                                                                                 aes(x = rank, y = mlog10Padj, label = TF), size = 4, nudge_x = 6, color = "black", max.overlaps = 50) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") + xlab("Rank Sorted TFs Enriched") + scale_color_gradientn(colors = paletteContinuous(set = "comet"))
ggUp
motifsDo <- peakAnnoEnrichment(seMarker = markerTest_C5HomvsC2C3WT, ArchRProj = atrialatacms, peakAnnotation = "Motif",
                               cutOff = "FDR <= 0.05 & Log2FC <= -0.25")
df2 <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df2 <- df2[order(df2$mlog10Padj, decreasing = TRUE),]
df2$rank <- seq_len(nrow(df2))
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + geom_point(size = 1) + ggrepel::geom_label_repel(data = df[rev(seq_len(30)), ],
                                                                                                                 aes(x = rank, y = mlog10Padj, label = TF), size = 4, nudge_x = 6, color = "black", max.overlaps = 50) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") + xlab("Rank Sorted TFs Enriched") + scale_color_gradientn(colors = paletteContinuous(set = "comet"))
ggDo
plotPDF(ggUp, ggDo, name = "atrialatacms-Motifs-C6HetvsC2C3WT-FDR>0.1_logfc0.25", width = 5, height = 10, ArchRProj = atrialatacms, addDOC = FALSE)
#to plot without boxes around labels
library(ggplot2)
library(ggrepel)
#To label Top 10
labels <- df[1:30,]
myplot <- ggplot(data = df, aes(x = df$rank,y = df$mlog10Padj)) + geom_point(aes(color = mlog10Padj)) + geom_text_repel(data = labels, aes(x = rank, y = mlog10Padj, label = TF), nudge_x = 0.5, nudge_y = 0.5, max.overlaps = Inf) +
  labs(x = "Rank", y = "Log10 Adjusted P-Value", title = "Motif Enrichment") + theme_classic() + scale_color_gradientn(colors = paletteContinuous(set = "comet")) 
myplot
##combined motif plots
# Add a new column to each data frame to indicate the source
# Adjust the rank for df2 (Down dataset) to be negative
df2$rank <- -df2$rank  # Make the rank negative for "Down"
# Combine the two data frames with the updated rank
df_combined <- rbind(df1, df2)
# Calculate the range of the ranks for setting limits
# Adjust any mlog10Padj values that are 0 by adding 0.01
df_combined$mlog10Padj <- ifelse(df_combined$mlog10Padj == 0, 0.001, df_combined$mlog10Padj)
# Create a volcano-like ggplot
gg_volcano <- ggplot(df_combined, aes(rank, mlog10Padj)) +
  geom_point(size = 1, aes(color = mlog10Padj)) +
  #ggrepel::geom_label_repel(aes(x = rank, y = mlog10Padj), #I'm labeling points manually in Illustrator
  #size = 2, nudge_x = 10, max.overlaps = 100) +
  theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched (Up vs Down)") +
  geom_vline(xintercept = 0, linetype = "dashed")  # Add a vertical line at x=0 for symmetry

# Display the plot
print(gg_volcano)
plotPDF(gg_volcano, name = "atrialatacms-combinedMotifs-C6HetvsC2C3WT-FDR>0.1_logfc0.25", ArchRProj = atrialatacms, addDOC = FALSE)


##MOTIF-PEAK-----------------------------------------------------------------------------------------------------------------------
#how to find differential peaks that contain TBX/MEIS motifs
#Vector of motif names
motif_names <- c("TBX5_781", "TBX4_784", "TBX15_782", "MEIS1_486", "MEIS2_471", "MEIS3_411")
# Initialize empty lists to store all dataframes
all_motifs_down <- list()
all_motifs_up <- list()
# Loop through each motif
for (motif in motif_names) {
  tryCatch({
    # Get matches for current motif
    mp <- getMatches(ArchRProj = atrialatacms, name = "Motif", annoName = motif)
    # Subset matching peaks
    sub_mp <- mp[which(mp@assays@data@listData$matches), ]
    # Find overlaps for down and up regulated peaks
    match_motif_down <- findOverlaps(sorted_markerTest_C5HomvsC2C3WT_downinHom_gr, sub_mp)
    match_motif_up <- findOverlaps(sorted_markerTest_C5HomvsC2C3WT_upinHom_gr, sub_mp)
    # Get matched ranges
    matched_ranges_down <- sub_mp[match_motif_down@to]
    matched_ranges_up <- sub_mp[match_motif_up@to]
    # Create data frames for down-regulated peaks
    motif_down <- data.frame(
      chr = seqnames(matched_ranges_down),
      start = start(matched_ranges_down),
      end = end(matched_ranges_down),
      strand = strand(matched_ranges_down),
      nearestGene = matched_ranges_down@rowRanges$nearestGene,
      peakType = matched_ranges_down@rowRanges$peakType,
      Log2FC = sorted_markerTest_C5HomvsC2C3WT_downinHom_gr[match_motif_down@from]$Log2FC,
      motif = motif)  # Add motif name as a column
    # Create data frames for up-regulated peaks
    motif_up <- data.frame(
      chr = seqnames(matched_ranges_up),
      start = start(matched_ranges_up),
      end = end(matched_ranges_up),
      strand = strand(matched_ranges_up),
      nearestGene = matched_ranges_up@rowRanges$nearestGene,
      peakType = matched_ranges_up@rowRanges$peakType,
      Log2FC = sorted_markerTest_C5HomvsC2C3WT_upinHom_gr[match_motif_up@from]$Log2FC,
      motif = motif)  # Add motif name as a column
    # Store dataframes in lists
    all_motifs_down[[motif]] <- motif_down
    all_motifs_up[[motif]] <- motif_up
    # Optional: Print progress message
    message(sprintf("Processed motif: %s", motif))
  }, error = function(e) {
    message(sprintf("Error processing motif %s: %s", motif, e$message))
  })
}
# Combine all dataframes
combined_down <- do.call(rbind, all_motifs_down)
combined_up <- do.call(rbind, all_motifs_up)
# Write combined results to files
fwrite(combined_down, row.names = TRUE, file = "TBXMEIS_Motifs_NearestGene_Peaks-DowninC5HomvsC2C3WT.csv")
fwrite(combined_up, row.names = TRUE, file = "TBXMEIS_Motifs_NearestGene_Peaks-UpinC5HomvsC2C3WT.csv")
# Optional: Print summary statistics
message("Summary of results:")
message(sprintf("Total down-regulated peaks across all motifs: %d", nrow(combined_down)))
message(sprintf("Total up-regulated peaks across all motifs: %d", nrow(combined_up)))


##how to also extract peak coordinates
motif_TBX5_test <- data.frame(
  chr = seqnames(matched_ranges),  # Extract chromosome numbers
  start = start(matched_ranges),          # Extract start positions
  end = end(matched_ranges),              # Extract end positions
  strand = strand(matched_ranges),        # Extract strand information
  peakType = matched_ranges@rowRanges$peakType,  # Assuming peakType is stored here
  Log2FC = sorted_markerTest_C6HetvsC3WT_downinHet_gr[match_motif@from]$Log2FC,
  nearestGene = matched_ranges@rowRanges$nearestGene)  # Nearest gene information
# Display the dataframe to check contents
print(motif_TBX5_test)


##check overlap with TBX5 ChIP peaks-----------------------------------------------------------------------------------------------------------------------
library(readxl)
library(dplyr)
chip_peaks <- read_excel("/path/to/TBX5-ChIP_d11_d20_d45_peaks.xlsx")
# Creating GRanges objects from data frames
gr_chip <- GRanges(seqnames = chip_peaks$chr,ranges = IRanges(start = chip_peaks$start, end = chip_peaks$end))
# Find overlaps for het
chipoverlaps_hetdown <- findOverlaps(gr_chip, markerTest_C6HetvsC2C3WT_downinHet_gr, type = "any", ignore.strand = F)
chipoverlaps_hetup <- findOverlaps(gr_chip, markerTest_C6HetvsC2C3WT_upinHet_gr, type = "any", ignore.strand = F)
chipoverlaps_homdown <- findOverlaps(gr_chip, markerTest_C5HomvsC2C3WT_downinHom_gr, type = "any", ignore.strand = F)
chipoverlaps_homup <- findOverlaps(gr_chip, markerTest_C5HomvsC2C3WT_upinHom_gr, type = "any", ignore.strand = F)
# Extract matching data
results_down_het <- data.frame(chr_chip = seqnames(gr_chip)[queryHits(chipoverlaps_hetdown)],
                               start_chip = start(gr_chip)[queryHits(chipoverlaps_hetdown)],
                               end_chip = end(gr_chip)[queryHits(chipoverlaps_hetdown)],
                               chr_subject = seqnames(markerTest_C6HetvsC2C3WT_downinHet_gr)[subjectHits(chipoverlaps_hetdown)],
                               start_subject = start(markerTest_C6HetvsC2C3WT_downinHet_gr)[subjectHits(chipoverlaps_hetdown)],
                               end_subject = end(markerTest_C6HetvsC2C3WT_downinHet_gr)[subjectHits(chipoverlaps_hetdown)])
print(results_down_het)
results_up_het <- data.frame(chr_chip = seqnames(gr_chip)[queryHits(chipoverlaps_hetup)],
                             start_chip = start(gr_chip)[queryHits(chipoverlaps_hetup)],
                             end_chip = end(gr_chip)[queryHits(chipoverlaps_hetup)],
                             chr_subject = seqnames(markerTest_C6HetvsC2C3WT_upinHet_gr)[subjectHits(chipoverlaps_hetup)],
                             start_subject = start(markerTest_C6HetvsC2C3WT_upinHet_gr)[subjectHits(chipoverlaps_hetup)],
                             end_subject = end(markerTest_C6HetvsC2C3WT_upinHet_gr)[subjectHits(chipoverlaps_hetup)])
print(results_up_het)
# Add the 'direction' column
results_up_het$direction <- "up"
results_down_het$direction <- "down"
# Merge the two dataframes
merged_results_het <- rbind(results_up_het, results_down_het)
#save
write.csv(merged_results_het, file = "het_overlap_chip.csv")
# Find overlaps for hom
chipoverlaps_homdown <- findOverlaps(gr_chip, markerTest_C5HomvsC2C3WT_downinhom_gr, type = "any", ignore.strand = F)
chipoverlaps_homup <- findOverlaps(gr_chip, markerTest_C5HomvsC2C3WT_upinhom_gr, type = "any", ignore.strand = F)
chipoverlaps_homdown <- findOverlaps(gr_chip, markerTest_C5HomvsC2C3WT_downinHom_gr, type = "any", ignore.strand = F)
chipoverlaps_homup <- findOverlaps(gr_chip, markerTest_C5HomvsC2C3WT_upinHom_gr, type = "any", ignore.strand = F)
# Extract matching data
results_down_hom <- data.frame(chr_chip = seqnames(gr_chip)[queryHits(chipoverlaps_homdown)],
                               start_chip = start(gr_chip)[queryHits(chipoverlaps_homdown)],
                               end_chip = end(gr_chip)[queryHits(chipoverlaps_homdown)],
                               chr_subject = seqnames(markerTest_C5HomvsC2C3WT_downinHom_gr)[subjectHits(chipoverlaps_homdown)],
                               start_subject = start(markerTest_C5HomvsC2C3WT_downinHom_gr)[subjectHits(chipoverlaps_homdown)],
                               end_subject = end(markerTest_C5HomvsC2C3WT_downinHom_gr)[subjectHits(chipoverlaps_homdown)])
print(results_down_hom)
results_up_hom <- data.frame(chr_chip = seqnames(gr_chip)[queryHits(chipoverlaps_homup)],
                             start_chip = start(gr_chip)[queryHits(chipoverlaps_homup)],
                             end_chip = end(gr_chip)[queryHits(chipoverlaps_homup)],
                             chr_subject = seqnames(markerTest_C5HomvsC2C3WT_upinHom_gr)[subjectHits(chipoverlaps_homup)],
                             start_subject = start(markerTest_C5HomvsC2C3WT_upinHom_gr)[subjectHits(chipoverlaps_homup)],
                             end_subject = end(markerTest_C5HomvsC2C3WT_upinHom_gr)[subjectHits(chipoverlaps_homup)])
print(results_up_hom)
# Add the 'direction' column
results_up_hom$direction <- "up"
results_down_hom$direction <- "down"
# Merge the two dataframes
merged_results_hom <- rbind(results_up_hom, results_down_hom)
#save
write.csv(merged_results_hom, file = "hom_overlap_chip.csv")