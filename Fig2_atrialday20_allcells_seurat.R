library(Seurat)
library(ggplot2)

##read in 10X matrix of day 20 samples
all.data <- Read10X(data.dir = "~/path/to/filtered_feature_bc_matrix")
atrial <- CreateSeuratObject(counts = all.data, project = "atrial")

##add metadata column for timepoint
atrial@meta.data$timepoint <- "day20"
head(atrial@meta.data)
tail(atrial@meta.data)

#add gemgroup and names as metadata columns
classification.vec <- as.numeric(gsub(".*-","", (colnames(atrial@assays$RNA@data))))
names(classification.vec) <- colnames(atrial@assays$RNA@data)
atrial <- AddMetaData(atrial, classification.vec, "gem.group")
tail(atrial@meta.data)

#create a column with replicate names 
classification.vec1 <- as.numeric(gsub(".*-","", (colnames(atrial@assays$RNA@data))))
names(classification.vec1) <- colnames(atrial@assays$RNA@data)
# rename gem group assignments as something different 
classification.vec1[classification.vec1==1] <- "WT-1"
classification.vec1[classification.vec1==2] <- "WT-2"
classification.vec1[classification.vec1==3] <- "Het-1"
classification.vec1[classification.vec1==4] <- "Het-2"
classification.vec1[classification.vec1==5] <- "Hom-1"
classification.vec1[classification.vec1==6] <- "Hom-2"
atrial <- AddMetaData(atrial, classification.vec1, "replicate")
head(atrial@meta.data)
tail(atrial@meta.data)

#create a column with genotype names 
classification.vec1 <- as.numeric(gsub(".*-","", (colnames(atrial@assays$RNA@data))))
names(classification.vec1) <- colnames(atrial@assays$RNA@data)
# rename gem group assignments as something different 
classification.vec1[classification.vec1==1] <- "WT"
classification.vec1[classification.vec1==2] <- "WT"
classification.vec1[classification.vec1==3] <- "Het"
classification.vec1[classification.vec1==4] <- "Het"
classification.vec1[classification.vec1==5] <- "Hom"
classification.vec1[classification.vec1==6] <- "Hom"
atrial <- AddMetaData(atrial, classification.vec1, "genotype")
head(atrial@meta.data)
tail(atrial@meta.data)

#create a column with batch names 
classification.vec1 <- as.numeric(gsub(".*-","", (colnames(atrial@assays$RNA@data))))
names(classification.vec1) <- colnames(atrial@assays$RNA@data)
# rename gem group assignments as something different 
classification.vec1[classification.vec1==1] <- "exp1"
classification.vec1[classification.vec1==2] <- "exp2"
classification.vec1[classification.vec1==3] <- "exp1"
classification.vec1[classification.vec1==4] <- "exp2"
classification.vec1[classification.vec1==5] <- "exp1"
classification.vec1[classification.vec1==6] <- "exp2"
atrial <- AddMetaData(atrial, classification.vec1, "batch")
head(atrial@meta.data)
tail(atrial@meta.data)

#QC
# Calculate percent mt- MAKE SURE TO HAVE MT for human data, mt for mouse data!!!!
atrial[["percent.mt"]] <- PercentageFeatureSet(atrial, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(atrial, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships
plot1 <- FeatureScatter(atrial, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "genotype", raster = F) + scale_x_continuous(limits=c(0, 150000))
plot2 <- FeatureScatter(atrial, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "genotype", raster = F)
plot1+plot2
##filter cells with low nCount_RNA and high percent.mt
atrial <- subset(atrial, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 4000 & nCount_RNA < 100000)
##filtered from 67356 to 48427 cells

#save unscaled object
save(atrial, file = "atrialdeep_rep1and2_unprocessed.Robj")

##normalizing and scaling using sc transform
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("glmGamPoi")
atrial <- SCTransform(atrial, method = "glmGamPoi",vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

##save object
save(atrial, file = "atrialdeep_rep1and2_nonormsct.Robj")

#dimensionality reduction
atrial <- RunPCA(atrial)
ElbowPlot(atrial, ndims = 50)

##batch correction
atrial <- RunHarmony(atrial, group.by.vars = "batch")

#cluster cells
atrial <- FindNeighbors(atrial, dims = 1:35,reduction = "harmony")
atrial <- FindClusters(atrial, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))

#umap
atrial <- RunUMAP(atrial, reduction = "harmony",dims = 1:35, n.neighbors = 5, min.dist = 0.25)
Idents(atrial) <- "SCT_snn_res.0.5"
DimPlot(atrial, reduction = "umap", label = T, label.size = 6, pt.size = 0.4, raster = F) + NoLegend()
DimPlot(object = atrial, reduction = "umap", group.by = "genotype",cols = c("red","steelblue2","black"), pt.size = 0.4, raster=F, shuffle =T, split.by = "genotype") + NoLegend()
DimPlot(object = atrial, reduction = "umap", group.by = "batch", pt.size = 0.4, raster=F, shuffle = T)
DimPlot(object = atrial, reduction = "umap", pt.size = 0.4, raster=F, split.by = "genotype", ncol = 2) + NoLegend()

#feature plot
FeaturePlot(atrial, features = c("TNNT2","TOP2A","RGS5","DLK1","EPCAM","CLDN4","NEFM","MSX1","TAGLN","DCN","AFP","SERPINA1","TTR","CCNO"),
            pt.size = 0.1,cols = c("grey","purple","darkblue"), raster=F)


##cellnumbers
cellnumbersbygenotype_res0.1 <- table(atrial@meta.data$genotype, atrial@meta.data$Celltype)
write.csv(cellnumbersbygenotype_res0.1, file = "atriald20_cellnumbersbygenotype_res0.1.csv")

cellnumbersbyreplicate <- table(atrial@meta.data$replicate, atrial@meta.data$Celltype)
write.csv(cellnumbersbyreplicate, file = "atriald20_cellnumbersbyreplicate.csv")


##CELLTYPE LABELS
Idents(atrial) <- "SCT_snn_res.0.1"
new.cluster.ids <- c("CMs", "Fibroblasts", "Epithelial-like", "Dividing cells", "Neural-like", "Undetermined","Fibroblasts-2",
                     "Endoderm", "Endoderm-2")
names(new.cluster.ids) <- levels(atrial)
atrial <- RenameIdents(atrial, new.cluster.ids)
atrial@meta.data$Celltype <- Idents(atrial)
DimPlot(atrial, reduction = "umap", label = TRUE, group.by = "Celltype",pt.size = 0.4, label.size = 5, repel = T) + NoLegend()

Idents(atrial) <- "SCT_snn_res.0.5"
new.cluster.ids <- c("CMs-1", "CMEs","Fibroblasts", "Epithelial-like", "Fibroblasts","Dividing cells", "Undetermined","Fibroblasts-2","Neural-like",
                     "MYH7+ CMs","Fibroblasts","Epithelial-like","NPPA+ CMs","Fibroblasts","Fibroblasts","Endoderm", "Epithelial-like","Fibroblasts",
                     "Endoderm-2","Fibroblasts")
names(new.cluster.ids) <- levels(atrial)
atrial <- RenameIdents(atrial, new.cluster.ids)
atrial@meta.data$Celltype2 <- Idents(atrial)
DimPlot(atrial, reduction = "umap", label = FALSE, group.by = "Celltype2",pt.size = 0.4, label.size = 5, repel = T) + NoLegend()

Idents(atrial) <- "SCT_snn_res.0.5"
new.cluster.ids <- c("CMs", "CMEs","Fibroblasts", "Epithelial-like", "Fibroblasts","Dividing cells", "Undetermined","Fibroblasts-2","Neural-like",
                     "MYH7+ CMs","Fibroblasts","Epithelial-like","CMs","Fibroblasts","Fibroblasts","Endoderm", "Epithelial-like","Fibroblasts",
                     "Endoderm-2","Fibroblasts")
names(new.cluster.ids) <- levels(atrial)
atrial <- RenameIdents(atrial, new.cluster.ids)
atrial@meta.data$Celltype3 <- Idents(atrial)
DimPlot(atrial, reduction = "umap", label = T, group.by = "Celltype3",pt.size = 0.4, label.size = 5, repel = T) + NoLegend()


##Dotplot
library(ggplot2)
Idents(atrial) <- "Celltype"
DimPlot(atrial, reduction = "umap", label = T, label.size = 6, pt.size = 0.4, raster = F) + NoLegend()
genes <- c("TNNT2","MYH6","POSTN","COL3A1","EPCAM","CLDN4","TOP2A","CENPF","NEFM","MSX1","TAGLN","DCN","AFP","SERPINA1","TTR","CCNO")
DotPlot(atrial, features = genes,dot.min = 0, dot.scale = 8) + scale_y_discrete(expand = c(0.5,0))+ scale_x_discrete(expand = c(0.1,0)) + theme(axis.text.x = element_text(angle = 90,size =10,face = 'italic'))

##save object
save(atrial, file = "atrialdeep_rep1and2_nonormsctharmony.Robj")
