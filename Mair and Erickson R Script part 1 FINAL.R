######################################
###### Mair and Erickson et al #######
######################################

##-- Load or install required packages
##-- NOTE: script is using Seurat v2.3.4
library(Seurat)
library(tidyverse)
library(MAST)
library(RColorBrewer)
library(doMC)
library(data.table)
library(reticulate)

##-- Point towards the correct Python code for running UMAP
use_python("/anaconda3/bin/python")

##-- Overview of loaded packages
sessionInfo()

##-- Set working directory
setwd('~/your_working_directory/')

##-- Data loading from CSV files derived from Seven Bridges
Abseq_1 <- read.csv(file = 'Cartridge1/Combined_PBMC_AbSeq_1_DBEC_MolsPerCell_with_SampleTag.csv', sep = ',', header = TRUE, row.names = 1, check.names = FALSE)
Abseq_2 <- read.csv(file = 'Cartridge2/Combined_PBMC_AbSeq_2_DBEC_MolsPerCell_with_SampleTag.csv', sep = ',', header = TRUE, row.names = 1, check.names = FALSE)


##-- Data pre-processing
#- RNA and protein matrices
Abseq_1P <- t(Abseq_1[, str_detect(string = colnames(Abseq_1), pattern = 'pAbO')]) # 
Abseq_1RNA <- (Abseq_1[, !str_detect(string = colnames(Abseq_1), pattern = 'pAbO')]) #
Abseq_1RNA <- t(Abseq_1RNA[, !str_detect(string = colnames(Abseq_1RNA), pattern = 'Sample')]) # 
dim(Abseq_1P)
dim(Abseq_1RNA)

Abseq_2P <- t(Abseq_2[, str_detect(string = colnames(Abseq_2), pattern = 'pAbO')]) # 
Abseq_2RNA <- (Abseq_2[, !str_detect(string = colnames(Abseq_2), pattern = 'pAbO')]) #
Abseq_2RNA <- t(Abseq_2RNA[, !str_detect(string = colnames(Abseq_2RNA), pattern = 'Sample')]) # 
dim(Abseq_2P)
dim(Abseq_2RNA)

#- Sample Tag information
Sample_Tag1 <- t(Abseq_1[, str_detect(string = colnames(Abseq_1), pattern = 'Sample_Name')]) # 
dim(Sample_Tag1)
Sample_Tag2 <- t(Abseq_2[, str_detect(string = colnames(Abseq_1), pattern = 'Sample_Name')]) # 
dim(Sample_Tag2)

#- Features - Proteins
P_names <- sapply(X = str_split(string = rownames(Abseq_1P), pattern = '\\|'), 
                  FUN = function(x) ifelse(x[1] == 'CD197', 
                                           paste(paste(x[1], x[2], sep = '|'), x[3], sep = '_'),
                                           paste(x[1], x[2], sep = '|'))) # Antigen|Gene (and clone for CD197 - CD197 is not unique (two different clones))
rownames(Abseq_1P) <- paste0(P_names, '_P')
rownames(Abseq_2P) <- paste0(P_names, '_P')

#- Features - Genes
RNA_names <- str_replace(string = rownames(Abseq_1RNA), pattern = '\\|[^(PolyA)]*',  replacement = '_')
RNA_names <- str_replace(string = RNA_names, pattern = '_?$', replacement = '_RNA')
rownames(Abseq_1RNA) <- RNA_names
rownames(Abseq_2RNA) <- RNA_names

#- Have a look at dimensions and the rownames of the RNA and antibody matrix
dim(Abseq_1RNA)
dim(Abseq_1P)
dim(Abseq_2RNA)
dim(Abseq_2P)


##-- Create Seurate Objects from RNA matrices, add Sample Tag information
Cartridge1 <- CreateSeuratObject(raw.data = Abseq_1RNA, min.cells = 1) # RNA as input
Cartridge2 <- CreateSeuratObject(raw.data = Abseq_2RNA, min.cells = 1) # RNA as input
Cartridge1@meta.data[, "sample.origin"] <- as.vector(Sample_Tag1)
Cartridge2@meta.data[, "sample.origin"] <- as.vector(Sample_Tag2)
PBMC_all <- MergeSeurat(Cartridge1, Cartridge2, add.cell.id1 = "Cartridge1", add.cell.id2 = "Cartridge2")

##-- Overview of number of cells for each donor/Sample Tag
table(PBMC_all@meta.data$sample.origin)

##-- QC Violin plot for number of genes and UMI
VlnPlot(object = PBMC_all, 
        features.plot = c('nGene', 'nUMI'), 
        nCol = 2, point.size.use=0.05, group.by = "sample.origin", x.lab.rot = T)



##-- Add protein expression data to the Seurat object
#- Make sure that the Cell IDs have the same names
colnames(Abseq_1P) <- paste("Cartridge1", colnames(Abseq_1P), sep="_")
colnames(Abseq_2P) <- paste("Cartridge2", colnames(Abseq_2P), sep="_")
Combined_protein <- cbind(Abseq_1P, Abseq_2P)
dim(Combined_protein)

#- Copy the object and add protein expression
PBMC_all <- SetAssayData(PBMC_all, assay.type = 'CITE', slot = 'raw.data', new.data = Combined_protein) # Even if it is not CITE-seq, use CITE as assay.type

#- Normalization and standard scaling, per Satija-lab recommendation using centered log-ratio (CLR)
PBMC_all <- NormalizeData(PBMC_all, assay.type = 'CITE', normalization.method = 'genesCLR')
PBMC_all <- ScaleData(object = PBMC_all, assay.type = 'CITE', display.progress = TRUE)


##-- Subset data to remove multiplets and undetermined cells
PBMC_all <- SetAllIdent(PBMC_all, id="sample.origin")
PBMC_clean <- SubsetData(PBMC_all, ident.remove=c("Multiplet", "Undetermined"))
table(PBMC_clean@meta.data$sample.origin)

##-- Copy object
PBMC <- PBMC_clean

##-- Processing of RNA information
#- Normalization - Standard log-normalization
PBMC <- NormalizeData(object = PBMC, assay.type = 'RNA', normalization.method = 'LogNormalize', scale.factor = 10000)

#- Calculate highly variable genes, not used, only for information
PBMC <- FindVariableGenes(object = PBMC, do.plot = TRUE, y.cutoff = 0.2)
length(PBMC@var.genes) # 122 genes 

#- Standard scaling
PBMC <- ScaleData(object = PBMC, vars.to.regress="nUMI", display.progress = TRUE, do.par=T, num.cores = 7)

#- PCA calculation, use all genes as input
PBMC <- RunPCA(object = PBMC, pc.genes = RNA_names, pcs.compute = 30, do.print = TRUE)
PCElbowPlot(object = PBMC, num.pc = 30) # 20-25 PCs seem enough as input for visualization and clustering
PCHeatmap(object = PBMC, pc.use = 1:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE) 
PCHeatmap(object = PBMC, pc.use = 11:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

#- Calculate UMAP dimensionality reduction
PBMC <- RunUMAP(object = PBMC, dims.use = 1:25, min_dist = 0.2)

#- Clustering based on RNA information only, resolution 0.8
PBMC <- FindClusters(object = PBMC, 
                     dims.use = 1:25,
                     print.output = 0, 
                     save.SNN = TRUE,
                     force.recalc = T,
                     resolution = 0.8)

DimPlot(PBMC, reduction.use = "umap", do.label=T, pt.size=0.1)
DimPlot(PBMC, reduction.use = "umap", do.label=F, pt.size=0.1, group.by="sample.origin", cols.use = brewer.pal(6, "Paired"))

#- Visualize all clusters (resolution 0.8)
DimPlot(PBMC, reduction.use = "umap", do.label=T, pt.size=0.1, group.by="res.0.8", 
        cols.use = c(brewer.pal(10, "Paired"), brewer.pal(8, "Dark2")), vector.friendly = TRUE)

#- Explore the data a little bit, look at Tetramer+ cells vs bulk cells, plot key transcripts and proteins
PBMC <- SetAllIdent(PBMC, id="sample.origin")
DimPlot(SubsetData(PBMC, ident.use = c("Donor1_EBV", "Donor3_EBV", "Donor2_EBV")), 
        reduction.use = "umap", do.label=F, pt.size=0.2, group.by="sample.origin", cols.use = brewer.pal(3, "Set1"))

DimPlot(SubsetData(PBMC, ident.use = c("Donor1_Bulk", "Donor3_Bulk", "Donor2_Bulk")), 
        reduction.use = "umap", do.label=F, pt.size=0.2, group.by="sample.origin", cols.use = brewer.pal(3, "Set1"))

PBMC <- SetAllIdent(PBMC, id="res.0.8")
markersRNA <- FindAllMarkers(object = PBMC, 
                             only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.5, # log-FC
                             return.thresh = 0.01, # Only return markers that have a p-value < 0.01
                             test.use = "MAST") # nUMI as proxy for CDR

top10PBMC_RNA <- markersRNA %>% 
  group_by(cluster) %>% 
  top_n(10, avg_logFC)

DoHeatmap(SubsetData(object = PBMC, max.cells.per.ident=500), 
          genes.use = top10PBMC_RNA$gene, 
          slim.col.label = TRUE, 
          remove.key = T, rotate.key = T, cex.row=10, cex.col=2, group.label.rot = T, group.cex = 10, group.spacing = 0.5,
          col.low = "#bdbbff", col.mid = "#FFFFFF", col.high = "#c20000", disp.min = -2.0, disp.max = 2.0)

FeaturePlot(PBMC, features.plot = c("CD3E_RNA", "CD4_RNA", "CD8A_RNA", "NCAM1_RNA", "CD14_RNA", "CD79A_RNA"),
            cols.use = c("#DCDCDC", "blue"), 
            reduction.use = 'umap')

FeaturePlot(PBMC, features.plot = c("CD4|CD4_P", "CD8|CD8A_P", "CD56|NCAM1_P", 
                                    "CD19|CD19_P", "CD16|FCGR3A_P", "CD14|CD14_P"),
            cols.use = brewer.pal(9, "Blues"), max.cutoff = "q99", pt.size = 0.1,
            reduction.use = 'umap', no.legend = F)


#- UMAP plots for FIGURE 2A, only donor Donor3 bulk cells, use resolution 0.8 for figure
PBMC <- SetAllIdent(PBMC, id="sample.origin")
PBMC_Donor3 <- SubsetData(PBMC, ident.use = "Donor3_Bulk")
PBMC_Donor3 <- RunUMAP(object = PBMC_Donor3, dims.use = 1:25, min_dist = 0.2)
PBMC_Donor3 <- FindClusters(object = PBMC_Donor3, 
                     dims.use = 1:25,
                     print.output = 0, 
                     save.SNN = TRUE,
                     force.recalc = T,
                     resolution = 0.8)
DimPlot(PBMC_Donor3, reduction.use = "umap", do.label=T, pt.size=0.5, cols.use = brewer.pal(12, "Paired"))

#- Explore the data from this one donor, plot key transcripts and proteins
FeaturePlot(PBMC_Donor3, features.plot = c("CD3E_RNA", "CD4_RNA", "CD8A_RNA", "FOXP3_RNA", "SELL_RNA", "CCR7_RNA", "CD79A_RNA",
                                            "KLRB1_RNA", "NCAM1_RNA", "CD14_RNA", "FCGR3A_RNA", "ITGAX_RNA"),
            cols.use = brewer.pal(9, "Blues"), max.cutoff = "q99",  pt.size = 0.1, 
            reduction.use = 'umap')
FeaturePlot(PBMC_Donor3, features.plot = c("CD3|CD3E_P", "CD4|CD4_P", "CD8|CD8A_P", "CD56|NCAM1_P", "CD11c|ITGAX_P", 
                                            "CD19|CD19_P", "CD45RA|PTPRC_P", "CD16|FCGR3A_P"),
            cols.use = brewer.pal(9, "Blues"), max.cutoff = "q99",  pt.size = 0.1, 
            reduction.use = 'umap')

#- Find differentially expressed genes and generate heatmap with top10 differentially expressed genes for FIGURE 2B
markersPBMC_Donor3 <- FindAllMarkers(object = PBMC_Donor3, 
                             only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.5, # log-FC
                             return.thresh = 0.01, # Only return markers that have a p-value < 0.01
                             test.use = "MAST") # nUMI as proxy for CDR

top10PBMC_PBMC_Donor3 <- markersPBMC_Donor3 %>% 
  group_by(cluster) %>% 
  top_n(10, avg_logFC)

DoHeatmap(SubsetData(object = PBMC_Donor3, max.cells.per.ident=100), 
          genes.use = top10PBMC_PBMC_Donor3$gene, 
          slim.col.label = TRUE, 
          remove.key = T, rotate.key = T, cex.row=10, cex.col=2, group.label.rot = T, group.cex = 10, group.spacing = 0.5,
          col.low = "#bdbbff", col.mid = "#FFFFFF", col.high = "#c20000", disp.min = -2.0, disp.max = 2.0)

#- Generate Violin and Ridge plots with Treg features for FIGURE 2C
VlnPlot(SubsetData(PBMC_Donor3, ident.use = c(0,1,7)), features.plot = c("CD4_RNA", "FOXP3_RNA", "CTLA4_RNA"),
        cols.use = c("#a6cee3", "#1f78b4", "#ff7f00"), point.size.use = 0, nCol = 3, single.legend = F)

RidgePlot(PBMC_Donor3, ident.include = c(0,1,7),
          features.plot = c("CD4|CD4_P", "CD127|IL7R_P", "CD25|IL2RA_P"),
          cols.use = c("#a6cee3", "#1f78b4", "#ff7f00"))





#- UMAP plot based on transcript for FIGURE 3A, only Bulk cells (right panel)
PBMC <- SetAllIdent(PBMC, id="sample.origin")
DimPlot(SubsetData(PBMC, ident.use = c("Donor1_Bulk", "Donor3_Bulk", "Donor2_Bulk")), 
        reduction.use = "umap", do.label=F, pt.size=0.1, group.by="sample.origin", cols.use = brewer.pal(3, "Set1"))



#- Copy object, generate UMAP plot based only on protein data using the adt distance matrix
PBMC_RNA_P <- PBMC
PBMC_RNA_P <- RunPCA(PBMC_RNA_P, pc.genes = rownames(Abseq_1P), assay.type = "CITE", pcs.compute = 30)
adt.data <- GetAssayData(PBMC_RNA_P, assay.type = 'CITE', slot = 'data')
adt.dist <- as.matrix(dist(t(adt.data)))
PBMC_RNA_P <- SetAllIdent(PBMC_RNA_P, id="res.0.8")
PBMC_RNA_P <- StashIdent(PBMC_RNA_P, 'rnaClusterID') #  Stash the RNA cluster IDs for later, in case we want to use them
PBMC_RNA_P <- RunUMAP(object = PBMC_RNA_P, distance.matrix = adt.dist, dims.use = 1:20, min_dist = 0.2)
DimPlot(PBMC_RNA_P, reduction.use = "umap", do.label=T, pt.size=0.1)

#- Clustering, based on protein information only, resolution 0.4
PBMC_RNA_P <- FindClusters(object = PBMC_RNA_P, 
                           dims.use = 1:20,
                           distance.matrix = adt.dist,
                           print.output = 0, 
                           resolution = 0.4,
                           save.SNN = TRUE)
DimPlot(PBMC_RNA_P, reduction.use = "umap", do.label=T, pt.size=0.2, group.by="res.0.4")

#- UMAP plot based on PROTEIN for FIGURE 3A, only Bulk cells (left panel)
PBMC_RNA_P <- SetAllIdent(PBMC_RNA_P, id = "sample.origin")
DimPlot(SubsetData(PBMC_RNA_P, ident.use = c("Donor1_Bulk", "Donor3_Bulk", "Donor2_Bulk")), 
        reduction.use = "umap", do.label=F, pt.size=0.1, group.by="sample.origin", cols.use = brewer.pal(3, "Set1"))


#- Explore the data a little bit, plot some key proteins
DimPlot(PBMC_RNA_P, reduction.use = "umap", do.label=T, pt.size=0.2, group.by="res.0.4",
        cols.use = c(brewer.pal(10, "Paired"), brewer.pal(8, "Dark2")))
DimPlot(PBMC_RNA_P, reduction.use = "umap", do.label=T, pt.size=0.2, group.by="rnaClusterID",
        cols.use = c(brewer.pal(10, "Paired"), brewer.pal(8, "Dark2")))
FeaturePlot(PBMC_RNA_P, features.plot = c("CD4|CD4_P", "CD8|CD8A_P", "CD56|NCAM1_P", 
                                            "CD19|CD19_P", "CD16|FCGR3A_P", "CD14|CD14_P"),
            cols.use = brewer.pal(9, "Blues"), max.cutoff = "q99",  pt.size = 0.1, nCol=3, reduction.use = 'umap')








#- Subsetting for myeloid data and T cell data
#- Clusters are selected based on expression of key proteins
PBMC <- SetAllIdent(PBMC, id="res.0.8")
PBMC.myeloid.all <- SubsetData(PBMC, ident.use=c(4,14,13, 17))
PBMC.CD3 <- SubsetData(PBMC, ident.use=c(0,12,6,15,1,8,9,3,5))

#- Run UMAP and plot data
PBMC.myeloid.all <- RunUMAP(object = PBMC.myeloid.all, dims.use = 1:20, min_dist = 0.2)
PBMC.myeloid.all <- FindClusters(object = PBMC.myeloid.all, 
                     dims.use = 1:20,
                     print.output = 0,
                     resolution = 0.3,
                     save.SNN = TRUE,
                     force.recalc = T)
DimPlot(PBMC.myeloid.all, reduction.use = "umap", do.label=T, pt.size=2, cols.use = brewer.pal(6, "Paired"))


#- Explore clusters, remove contaminating T cell cluster 4 from the myeloid object
VlnPlot(PBMC.myeloid.all, features.plot = c("CD3E_RNA", "CD4_RNA", "HLA-DRA_RNA"),
        cols.use = brewer.pal(6, "Paired"), point.size.use = 0)
RidgePlot(PBMC.myeloid.all, features.plot = c("CD8|CD8A_P", "CD4|CD4_P", "HLA-DR|CD74_P"),
          cols.use = brewer.pal(6, "Paired"))
PBMC.myeloid <- SubsetData(PBMC.myeloid.all, ident.use=c(0:3, 5), subset.raw = T)


#- Rerun UMAP with cleaned object, clustering based on transcript only, resolution 0.3 to avoid overclustering
PBMC.myeloid <- RunUMAP(object = PBMC.myeloid, dims.use = 1:20, min_dist = 0.2)
PBMC.myeloid <- FindClusters(object = PBMC.myeloid, 
                             dims.use = 1:20,
                             print.output = 0,
                             resolution = 0.3,
                             save.SNN = TRUE,
                             force.recalc = T)

#- UMAP plots for FIGURE 5A, based on transcript
DimPlot(PBMC.myeloid, reduction.use = "umap", do.label=T, pt.size=2, cols.use = brewer.pal(6, "Paired"), group.by = "res.0.3")

#- Violin plots for FIGURE 5B
FeaturePlot(PBMC.myeloid, features.plot = c("CD14|CD14_P", "CD16|FCGR3A_P"),
            cols.use = brewer.pal(9, "Blues"), max.cutoff = "q99", pt.size = 2,
            reduction.use = 'umap')

#- Find differentially expressed genes and generate heatmap with top10 differentially expressed genes
markers.myeloid <- FindAllMarkers(object = PBMC.myeloid, 
                             only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.5, # log-FC
                             return.thresh = 0.01, # Only return markers that have a p-value < 0.01
                             test.use = "MAST") # nUMI as proxy for CDR

top10myeloid <- markers.myeloid %>% 
  group_by(cluster) %>% 
  top_n(10, avg_logFC)

#- Heatmap for FIGURE 5C
DoHeatmap(SubsetData(PBMC.myeloid, ident.use = c(0:4), max.cells.per.ident = 150), 
          genes.use = top10myeloid$gene, 
          slim.col.label = TRUE, 
          remove.key = T, rotate.key = T, cex.row=10, cex.col=2, group.label.rot = T, group.cex = 10, group.spacing = 0.5,
          col.low = "#bdbbff", col.mid = "#FFFFFF", col.high = "#c20000", disp.min = -4.0, disp.max = 4.0)

#- Ridgeplots for FIGURE 5D
RidgePlot(SubsetData(PBMC.myeloid, ident.use = c(0:4)), 
          features.plot = c("CD11c|ITGAX_P", "CD123|IL3RA_P", "CD14|CD14_P", "CD16|FCGR3A_P", 
          "CD1c|CD1C_P", "CD11b|ITGAM_P", "HLA-DR|CD74_P", "CD86|CD86_P", "CD163|CD163_P"),
          cols.use = brewer.pal(5, "Paired"))

#- Violin plots for FIGURE 5F
VlnPlot(PBMC.myeloid, features.plot = c("CD14_RNA", "FCGR3A_RNA", "CD1C_PolyA_1_RNA", "GZMB_RNA",
                                        "IL1B_RNA", "TNF_RNA","CXCL3_RNA", "CCL4_RNA"),
        cols.use = brewer.pal(5, "Paired"), point.size.use = 0,  nCol = 4, single.legend = F)






#- Exploring the T cell object - Figures 3D through H were generated in Seurat 3, see separate script
PBMC.CD3 <- RunUMAP(object = PBMC.CD3, dims.use = 1:20, min_dist = 0.2)
PBMC.CD3 <- FindClusters(object = PBMC.CD3, 
                             dims.use = 1:20,
                             print.output = 0,
                             resolution = 0.4,
                             save.SNN = TRUE,
                             force.recalc = T, k.param = 30)

DimPlot(PBMC.CD3, reduction.use = "umap", do.label=T, pt.size=1, group.by = "res.0.4", cols.use = brewer.pal(8, "Paired"))
PBMC.CD3 <- SubsetData(PBMC.CD3, ident.remove = 7)
DimPlot(PBMC.CD3, reduction.use = "umap", do.label=T, pt.size=1, group.by = "res.0.4", cols.use = brewer.pal(8, "Paired"))

DimPlot(PBMC.CD3, reduction.use = "umap", do.label=F, pt.size=1, group.by = "sample.origin", 
        cols.use = c("#A9A9A9", "#FF0000", "#A9A9A9", 
        "#FF0000", "#A9A9A9", "#FF0000"))

PBMC.CD3 <- SetAllIdent(PBMC.CD3, id="sample.origin")
DimPlot(SubsetData(PBMC.CD3, ident.use = c("Donor1_EBV", "Donor3_EBV", "Donor2_EBV")), 
        reduction.use = "umap", do.label=F, pt.size=2, group.by = "sample.origin", cols.use = brewer.pal(3, "Set1"))
DimPlot(SubsetData(PBMC.CD3, ident.use = c("Donor1_Bulk", "Donor3_Bulk", "Donor2_Bulk")), 
        reduction.use = "umap", do.label=T, pt.size=1, group.by = "sample.origin", cols.use = brewer.pal(3, "Set1"))

PBMC.CD3 <- SetAllIdent(PBMC.CD3, id="res.0.4")
markers.CD3 <- FindAllMarkers(object = PBMC.CD3, 
                                  only.pos = F,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.5, # log-FC
                                  return.thresh = 0.01, # Only return markers that have a p-value < 0.01
                                  test.use = "MAST") # nUMI as proxy for CDR

top10CD3 <- markers.CD3 %>% 
  group_by(cluster) %>% 
  top_n(10, avg_logFC)

DoHeatmap(SubsetData(PBMC.CD3, max.cells.per.ident = 500), 
          genes.use = top10CD3$gene,
          slim.col.label = TRUE,
          remove.key = T, rotate.key = T, cex.row=10, cex.col=2, group.label.rot = T, group.cex = 10, group.spacing = 0.5,
          col.low = "#e30000", col.mid = "#FFFFFF", col.high = "#194eb0", disp.min = -2.0, disp.max = 2.0)

RidgePlot(PBMC.CD3, 
          features.plot = c("CD69|CD69_P", "CD27|CD27_P"),
          cols.use = brewer.pal(7, "Paired"))

VlnPlot(PBMC.CD3, features.plot = c("CD69_RNA", "CD27_RNA"),
        cols.use = brewer.pal(7, "Paired"), point.size.use = 0, nCol = 2, single.legend = T)




#- Working with the 10x data of Donor3, directly import Chromium 10x object and annotated Rhapsody object
Chromium10x <- read_rds("Chromium10x.rds")
Chromium10x
Rhapsody <- read_rds("Rhapsody.rds")
Rhapsody

DimPlot(Chromium10x, reduction.use = "umap", do.label=F, pt.size=0.5, cols.use = brewer.pal(12, "Paired"))
DimPlot(Rhapsody, reduction.use = "umap", do.label=F, pt.size=0.5, cols.use = brewer.pal(12, "Paired"))


#- Generate Heatmaps for SUPPLEMENTARY FIGURE 2B, genes are based on top diff genes from Rhapsody clusters
DoHeatmap(SubsetData(object = Rhapsody, max.cells.per.ident=500, 
                     ident.use = c("Naive CD4 T cells", "Memory CD8 T cells", "NK cells", "CD14 Monocytes")),
          genes.use = c("IL32_RNA", "IL7R_RNA", "TRAT1_RNA", "TIMP1_RNA", "LEF1_RNA", "SELL_RNA",
                        "PIK3IP1_RNA", "CCR7_RNA", "CD3D_RNA", "KLRB1_RNA", 
                        "GZMH_RNA", "GNLY_RNA","CCL5_RNA", "ZNF683_RNA",
                        "KLRF1_RNA", "GZMB_RNA", "PRF1_RNA", "FCGR3A_RNA", "CD10_PolyA_1_RNA",
                        "S100A9_RNA", "S100A12_RNA", "CD14_RNA", "LYZ_RNA", "CXCL8_RNA"),
          slim.col.label = TRUE,
          remove.key = T, rotate.key = T, cex.row=7, cex.col=2, group.label.rot = T, group.cex = 10, use.scaled = T, 
          col.low = "#bdbbff", col.mid = "#FFFFFF", col.high = "#c20000", disp.min = -2.0, disp.max = 2.0)

DoHeatmap(SubsetData(object = Chromium10x, max.cells.per.ident=500, 
                     ident.use = c("CD4 T cells", "CD8 T cells", "NK cells", "Monocytes")),
          genes.use = c("IL32", "IL7R", "TRAT1", "TIMP1", "LEF1", "SELL",
                        "PIK3IP1", "CCR7", "CD3D", "KLRB1", 
                        "GZMH", "GNLY","CCL5", "ZNF683",
                        "KLRF1", "GZMB", "PRF1", "FCGR3A", "CD10_PolyA_1",
                        "S100A9", "S100A12", "CD14", "LYZ", "CXCL8"),
          slim.col.label = TRUE,
          remove.key = T, rotate.key = T, cex.row=7, cex.col=2, group.label.rot = T, group.cex = 10, use.scaled = T, 
          col.low = "#bdbbff", col.mid = "#FFFFFF", col.high = "#c20000", disp.min = -2.0, disp.max = 2.0)





# Generate FCS-Files using premessa and FlowCore, FCS file downstream analysis in FlowJo v10.5
library(premessa)
library(flowCore)

table(PBMC_RNA_P@meta.data$sample.origin)
PBMC_RNA_P <- SetAllIdent(PBMC_RNA_P, id = "sample.origin")
Donor1_Bulk <- SubsetData(PBMC_RNA_P, ident.use = "Donor1_Bulk")
Donor3_Bulk <- SubsetData(PBMC_RNA_P, ident.use = "Donor3_Bulk")
Donor2_Bulk <- SubsetData(PBMC_RNA_P, ident.use = "Donor2_Bulk")
EBV <- SubsetData(PBMC_RNA_P, ident.use = c("Donor1_EBV", "Donor3_EBV", "Donor2_EBV"))
All.clean <- PBMC_RNA_P

Donor1_Bulk.ab <- t(Donor1_Bulk@assay$CITE@raw.data)
Donor3_Bulk.ab <- t(Donor3_Bulk@assay$CITE@raw.data)
Donor2_Bulk.ab <- t(Donor2_Bulk@assay$CITE@raw.data)

Donor1_Bulk.fcs <- as_flowFrame(Donor1_Bulk.ab)
Donor3_Bulk.fcs <- as_flowFrame(Donor3_Bulk.ab)
Donor2_Bulk.fcs <- as_flowFrame(Donor2_Bulk.ab)
EBV.fcs <- as_flowFrame(t(EBV@assay$CITE@raw.data))
All.fcs <- as_flowFrame(t(All.clean@assay$CITE@raw.data))

write.FCS(Donor1_Bulk.fcs, filename="Donor1_Bulk.fcs")
write.FCS(Donor3_Bulk.fcs, filename="Donor3_Bulk.fcs")
write.FCS(Donor2_Bulk.fcs, filename="Donor2_Bulk.fcs")
write.FCS(EBV.fcs, filename="EBV_all.fcs")
write.FCS(All.fcs, filename="Abseq_all.fcs")





