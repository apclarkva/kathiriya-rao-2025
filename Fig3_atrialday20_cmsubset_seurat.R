library(Seurat)
library(ggplot2)

##subset all TNNT2+ clusters from full object
load("atrial.Robj")
Idents(atrial) <- "SCT_snn_res.0.1"
atrial20cms <- subset(atrial, idents = 0)
atrial20cms <- CreateSeuratObject(counts = atrial20cms@assays$RNA@counts, project = "atrial-d20")

#save
save(atrial20cms, file = "atrial20cms_unprocessed.Robj")

##re-normalizing and scaling using sc transform 
atrial20cms <- SCTransform(atrial20cms, method = "glmGamPoi",vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#dimensionality reduction
atrial20cms <- RunPCA(atrial20cms)
ElbowPlot(atrial20cms, ndims = 50)

##batch correction
atrial20cms <- RunHarmony(atrial20cms, group.by.vars = "batch")

#cluster cells
atrial20cms <- FindNeighbors(atrial20cms, dims = 1:35,reduction = "harmony")
atrial20cms <- FindClusters(atrial20cms, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))


#umap
atrial20cms <- RunUMAP(atrial20cms, reduction = "harmony",dims = 1:35, n.neighbors = 5, min.dist = 0.1)
Idents(atrial20cms) <- "SCT_snn_res.0.7"
DimPlot(atrial20cms, reduction = "umap", label = T, label.size = 6, pt.size = 0.4, raster = F) + NoLegend()
DimPlot(object = atrial20cms, reduction = "umap", group.by = "genotype",cols = c("red","steelblue2","black"), pt.size = 0.4, raster=F, shuffle = T) + NoLegend()
DimPlot(object = atrial20cms, reduction = "umap", group.by = "batch",pt.size = 0.4, raster=F, shuffle = T)

#featureplot
FeaturePlot(atrial20cms, features = c("MYH6","TOP2A","PITX2","HAMP","NPPA","MYH7","IRX4","RSPO3","COL3A1","HBD"),pt.size = 0.1,cols = c("grey","purple","darkblue"), raster=F)

##CELLTYPE LABELS #
Idents(atrial20cms) <- "SCT_snn_res.0.7"
#first print the levels - "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13"
#if they are out of order for some reason, label accordingly!!!
levels(atrial20cms)
new.cluster.ids <- c("CMEs-1", "CMEs-2","RA","Atrial Precursors","LA","NPPA+ CMs","Undetermined","Vent-like CMs", "HBD+", 
                     "AVCMs", "Vent-like CMs", "Dividing CMs","CMEs-3","Undetermined")
names(new.cluster.ids) <- levels(atrial20cms)
atrial20cms <- RenameIdents(atrial20cms, new.cluster.ids)
atrial20cms@meta.data$Celltype <- Idents(atrial20cms)
DimPlot(atrial20cms, reduction = "umap", label = TRUE, group.by = 'Celltype',pt.size = 0.6, label.size = 6, repel = T, 
        cols = c("magenta","salmon","darkgreen","green","olivedrab","dodgerblue4","grey","gold3","chocolate","purple","pink2","maroon3","grey")) + NoLegend()

##save
save(atrial20cms, file = "atriald20cms_sctharmony.Robj")

##cellnumbers
cellnumbersbygenotype_res0.7 <- table(atrial20cms@meta.data$genotype, atrial20cms@meta.data$SCT_snn_res.0.7)
write.csv(cellnumbersbygenotype_res0.7, file = "atriald20cms_cellnumbersbygenotype_res0.7.csv")


##diff gene test between MYH6+ clusters vs others
Idents(atrial20cms) <- "SCT_snn_res.0.2"
myh6cms <- subset(atrial20cms, idents = c(0,3))
othercms <- subset(atrial20cms, idents = c(1,2,4,5,6))
Idents(myh6cms) <- "genotype"
Idents(othercms) <- "genotype"
Idents(atrial20cms, cells = WhichCells(myh6cms, idents = "WT")) <- "wtmyh6cms"
Idents(atrial20cms, cells = WhichCells(othercms, idents = "Het")) <- "hetothercms"
Idents(atrial20cms, cells = WhichCells(othercms, idents = "Hom")) <- "homothercms"
wtmyh6cms_vs_hetothercms <- FindMarkers(atrial20cms, ident.1 = "hetothercms", ident.2 = "wtmyh6cms")
wtmyh6cms_vs_homothercms <- FindMarkers(atrial20cms, ident.1 = "homothercms", ident.2 = "wtmyh6cms")
write.csv(wtmyh6cms_vs_hetothercms, file = "hetothercms_vs_wtmyh6cms.csv")
write.csv(wtmyh6cms_vs_homothercms, file = "homothercms_vs_wtmyh6cms.csv")


#ENHANCED VOLCANO PLOT ------------------------------------------------------------------------------------------------------------
library(EnhancedVolcano)
atrial20cms_hetvswt <- FindMarkers(atrial20cms, ident.1 = "hetothercms", ident.2 = "wtmyh6cms", logfc.threshold = 0.1)
atrial20cms_homvswt <- FindMarkers(atrial20cms, ident.1 = "homothercms", ident.2= "wtmyh6cms", logfc.threshold = 0.1)
atrial20cms_hetvswt$X <- row.names(atrial20cms_hetvswt)
atrial20cms_homvswt$X <- row.names(atrial20cms_homvswt)
genes_to_label_het <- scan(text = "DES NPPA HAMP MYH6 NR2F2 NR2F1 HEY1 PDGFRA ANGPT1 PITX2 MEF2C BMP2 RSPO3 MYH7 ACTA1 COL3A1 FN1
                           FBN2 GOPC HCN4 TBX3 TBX2-AS1 LBH", what = "")
genes_to_label_hom <- scan(text = "NPPA	CKM	HAMP	TRH	MB	TECRL	EMC10	TFPI2	RRAD	NR2F2	GNG11	ANGPT1
                           	PEG10	TUBA1A	STMN1	TMSB4X	GPC3	FN1	ANXA2	RSPO3	S100A10	COL3A1", what = "")
EnhancedVolcano(
  atrial20cms_homvswt,
  lab = atrial20cms_homvswt$X,
  selectLab = genes_to_label_hom,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = "hom vs WT CMs",
  pCutoff = .05,
  FCcutoff = 0.25,
  pointSize = 2.0,
  labSize = 5.0,
  xlim = c(-3.5,3.5),
  ylim = c(0,115),
  axisLabSize = 10,
  titleLabSize = 10,
  labFace = "italic",
  col = c('grey','grey','grey','purple3'),
  max.overlaps = 50,
  boxedLabels = FALSE,
  drawConnectors = TRUE,
  widthConnectors = 0.25,
  directionConnectors = 'y',
  arrowheads = FALSE,
  legendPosition = 'none')
ggsave(paste(wd,st, proj,"_c9_Volcano_","log25",".pdf",sep = ""), width = 12,height = 12)

##dose dependent genes
##combine gene lists 
names(wtmyh6cms_vs_hetothercms)[names(wtmyh6cms_vs_hetothercms) == 'avg_log2FC'] <- 'avg_log2FC.het'
names(wtmyh6cms_vs_homothercms)[names(wtmyh6cms_vs_homothercms) == 'avg_log2FC'] <- 'avg_log2FC.hom'
common_hethom <- merge(wtmyh6cms_vs_hetothercms, wtmyh6cms_vs_homothercms, by.x = 0, by.y = 0)
write.csv(common_hethom, file = "atrial20cms_common_hethom.csv")

##removed non-signifcant genes from above common list
##now to plot with only common CHD and AF genes that are dose dependent
chd <- scan(text = "ACTB	CDKN1C	CFC1	COL2A1	GPC3	HAND1	JAG1	MEIS2	NR2F2	RABGAP1L	RPS24	SMAD6	VCAN	COL1A1	
                       COL3A1	MYH6	MYOCD	PDGFRA	PRKAG2	TBX3", what ="")
af <-  scan(text = "ACTG1	ACTN2	ADM	ANGPT1	ATP1A1	ATP1A2	ATP2A2	CACNA1D	CKM	CMYA5	CNTN5	CPNE5	CSRP3	
                  CST3	CTSL	EIF4EBP1	ERBB4	FBXO32	FN1	GAS5	GPX3	GSTP1	HAND2-AS1	HCN4	HSPB2 LGALS3	LRRC10	
                  MLLT3	MYL3	MYO18B	NPPA	PLAT	PLXNA4	SLC25A4	SLC8A1	SPEG	TAGLN	TIMP1	TNNI3	TNNT2	TPM1
                       TRIM55	VIM", what = "")
#run differential again to arrange by logfc of het
tempchd <- FindMarkers(atrial20cms, features = chd,ident.1 = "wtmyh6cms", ident.2 = "hetothercms")
##sort by logfc
# Sort the data frame by 'avg_log2FC' from largest to smallest
sorted_chd <- tempchd[order(tempchd$avg_log2FC, decreasing = TRUE), ]
#run differential again to arrange by logfc of het
tempaf <- FindMarkers(atrial20cms, features = af,ident.1 = "wtmyh6cms", ident.2 = "hetothercms")
##sort by logfc
# Sort the data frame by 'avg_log2FC' from largest to smallest
sorted_af <- tempaf[order(tempaf$avg_log2FC, decreasing = TRUE), ]
#merge <- c(row.names(sorted_af), row.names(sorted_chd))
subset <- subset(atrial20cms, idents = c("wtmyh6cms", "hetothercms", "homothercms"))
#plot one dotplot for chd
DotPlot(subset, features = row.names(sorted_chd), dot.min = 0, dot.scale = 6) + scale_y_discrete(expand = c(4,0), limits=c("wtmyh6cms","hetothercms","homothercms")) + scale_x_discrete(expand = c(0.025,0)) + theme(axis.text.x = element_text(angle = 90,size =8,face = 'italic')) + NoLegend()
#plot one dotplot for af
DotPlot(subset, features = row.names(sorted_af)[1:21], dot.min = 0, dot.scale = 6) + scale_y_discrete(expand = c(4,0), limits=c("wtmyh6cms","hetothercms","homothercms")) + scale_x_discrete(expand = c(0.025,0)) + theme(axis.text.x = element_text(angle = 90,size =8,face = 'italic')) + NoLegend()
DotPlot(subset, features = row.names(sorted_af)[22:43], dot.min = 0, dot.scale = 6) + scale_y_discrete(expand = c(4,0), limits=c("wtmyh6cms","hetothercms","homothercms")) + scale_x_discrete(expand = c(0.025,0)) + theme(axis.text.x = element_text(angle = 90,size =8,face = 'italic')) + NoLegend()
