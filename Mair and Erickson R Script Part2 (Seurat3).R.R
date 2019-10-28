###--- Load packages
library(Seurat)
library(tidyverse)
library(data.table)
library(MAST)
library(circlize)
library(cowplot)
library(SingleCellExperiment)
library(zinbwave)
library(BiocParallel)
library(doMC)
library(scamp)
library(reticulate)
library(RColorBrewer)

###--- Data Import
Abseq.cartridge_1 <- fread(input = '~/Desktop/R work/Combined_PBMC_AbSeq_1_DBEC_MolsPerCell_EDIT.csv')
Abseq.cartridge_2 <- fread(input = '~/Desktop/R work/Combined_PBMC_AbSeq_2_DBEC_MolsPerCell_EDIT.csv')

###--- Data processing
#- RNA matrices
Abseq.cartridge.RNA_1 <- Abseq.cartridge_1[, str_detect(string = colnames(Abseq.cartridge_1), pattern = 'Reference_end|PolyA_1|PolyA_2'), with = FALSE] # 499 genes
Abseq.cartridge.RNA_1[, Cell_ID := Abseq.cartridge_1$Cell_Index]
Abseq.cartridge.RNA_1[, Sample_Name := Abseq.cartridge_1$Sample_Name]
Abseq.cartridge.RNA_1[, Sample_Tag := Abseq.cartridge_1$Sample_Tag]
Abseq.cartridge.RNA_1[, Batch := 1]
Abseq.cartridge.RNA_2 <- Abseq.cartridge_2[, str_detect(string = colnames(Abseq.cartridge_2), pattern = 'Reference_end|PolyA_1|PolyA_2'), with = FALSE] # 499 genes
Abseq.cartridge.RNA_2[, Cell_ID := Abseq.cartridge_2$Cell_Index]
Abseq.cartridge.RNA_2[, Sample_Name := Abseq.cartridge_2$Sample_Name]
Abseq.cartridge.RNA_2[, Sample_Tag := Abseq.cartridge_2$Sample_Tag]
Abseq.cartridge.RNA_2[, Batch := 2]
#- Ab matrices
Abseq.cartridge.Ab_1 <- Abseq.cartridge_1[, str_detect(string = colnames(Abseq.cartridge_1), pattern = 'pAbO'), with = FALSE] # 42 markers
Abseq.cartridge.Ab_1[, Cell_ID := Abseq.cartridge_1$Cell_Index]
Abseq.cartridge.Ab_1[, Sample_Name := Abseq.cartridge_1$Sample_Name]
Abseq.cartridge.Ab_1[, Sample_Tag := Abseq.cartridge_1$Sample_Tag]
Abseq.cartridge.Ab_1[, Batch := 1]
Abseq.cartridge.Ab_2 <- Abseq.cartridge_2[, str_detect(string = colnames(Abseq.cartridge_2), pattern = 'pAbO'), with = FALSE] # 42 markers
Abseq.cartridge.Ab_2[, Cell_ID := Abseq.cartridge_2$Cell_Index]
Abseq.cartridge.Ab_2[, Sample_Name := Abseq.cartridge_2$Sample_Name]
Abseq.cartridge.Ab_2[, Sample_Tag := Abseq.cartridge_2$Sample_Tag]
Abseq.cartridge.Ab_2[, Batch := 2]
#- Features - RNA
identical(names(Abseq.cartridge.RNA_1), names(Abseq.cartridge.RNA_2)) # TRUE
RNA_names <- str_replace(string = names(Abseq.cartridge.RNA_1)[1:499], pattern = '\\|[^(PolyA)]*',  replacement = '_')
RNA_names <- str_replace(string = RNA_names, pattern = '_?$',  replacement = '.rna')
setnames(x = Abseq.cartridge.RNA_1, old = names(Abseq.cartridge.RNA_1), new = c(RNA_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
setnames(x = Abseq.cartridge.RNA_2, old = names(Abseq.cartridge.RNA_2), new = c(RNA_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
#- Features - Ab
identical(names(Abseq.cartridge.Ab_1), names(Abseq.cartridge.Ab_2)) # TRUE
Ab_names <- sapply(X = str_split(string = names(Abseq.cartridge.Ab_1)[1:42], pattern = '\\|'), 
                   FUN = function(x) paste(x[1], x[2], sep = '|'))
Ab_names <- str_replace(string = Ab_names, pattern = '[|]', replacement = '.')
Ab_names <- paste(Ab_names, 'ab', sep = '.')
setnames(x = Abseq.cartridge.Ab_1, old = names(Abseq.cartridge.Ab_1), new = c(Ab_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
setnames(x = Abseq.cartridge.Ab_2, old = names(Abseq.cartridge.Ab_2), new = c(Ab_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))

###--- Seurat v3
#- Setup a Seurat object
Abseq.RNA <- bind_rows(Abseq.cartridge.RNA_1, Abseq.cartridge.RNA_2)
Abseq.RNA.rawData <- t(Abseq.RNA[, -c('Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch')])
colnames(Abseq.RNA.rawData) <- paste('cell', 1:dim(Abseq.RNA.rawData)[2], sep = '-')
PBMC.RNA <- CreateSeuratObject(counts = Abseq.RNA.rawData, min.cells = 1, min.features = 1, assay = 'RNA', project = 'Ab-seq - Mair et al., 2019')
PBMC.RNA # 476 genes across 29,033 cells (18 genes have been removed)

#- Add features - meta.data
PBMC.RNA@meta.data$cell.id <- rownames(PBMC.RNA@meta.data)
PBMC.RNA@meta.data$sample.tag <- Abseq.RNA$Sample_Tag
PBMC.RNA@meta.data$sample.name <- Abseq.RNA$Sample_Name
PBMC.RNA@meta.data$batch <- Abseq.RNA$Batch
PBMC.RNA@meta.data$donor <- plyr::mapvalues(x = PBMC.RNA@meta.data$sample.tag, 
                                            from = unique(PBMC.RNA@meta.data$sample.tag), 
                                            to = c('Donor1', 'Multiplet', 'Donor2', 'Donor3', 'Donor3', 'Donor2', 'Donor1', 'Undetermined'))
PBMC.RNA@meta.data$celltype <- plyr::mapvalues(x = PBMC.RNA@meta.data$sample.tag, 
                                               from = unique(PBMC.RNA@meta.data$sample.tag), 
                                               to = c('Bulk', 'Multiplet', 'Bulk', 'Bulk', 'EBV', 'EBV', 'EBV', 'Undetermined'))

#- Add Ab expression to Seurat object
Abseq.Ab <- bind_rows(Abseq.cartridge.Ab_1, Abseq.cartridge.Ab_2)
Abseq.Ab.rawData <- t(Abseq.Ab[, -c('Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch')])
colnames(Abseq.Ab.rawData) <- paste('cell', 1:dim(Abseq.Ab.rawData)[2], sep = '-')
Abseq.Ab.rawData <- Abseq.Ab.rawData[, colnames(PBMC.RNA)]
PBMC.RNA[['Ab']] <- CreateAssayObject(counts = Abseq.Ab.rawData)


#- Subsetting and QC
PBMC.RNA <- SetIdent(object = PBMC.RNA, value = 'celltype')
PBMC.RNA <- SubsetData(object = PBMC.RNA, ident.use = c('Bulk', 'EBV')) # 481 genes across 27,258 cells (1,775 cells have been removed)
FeatureScatter(object = PBMC.RNA, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', cex.use = 0.5) # remove cells having nUMI > 2,000 and nGene < 115; and nUMI > 4,000
PBMC.RNA <- SubsetData(object = PBMC.RNA, cells = setdiff(colnames(PBMC.RNA), rownames(PBMC.RNA@meta.data)[which(PBMC.RNA@meta.data$nCount_RNA > 2000 & PBMC.RNA@meta.data$nFeature_RNA < 115)]))
PBMC.RNA <- SubsetData(object = PBMC.RNA, cells = setdiff(colnames(PBMC.RNA), rownames(PBMC.RNA@meta.data)[which(PBMC.RNA@meta.data$nCount_RNA > 4000)])) # 481 genes across 27,242 cells 

#- RNA Normalization
PBMC.RNA <- NormalizeData(object = PBMC.RNA, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = 10000)

#- Ab Normalization
PBMC.RNA <- NormalizeData(object = PBMC.RNA, assay = 'Ab', normalization.method = 'CLR')


#- Highly variable genes
PBMC.RNA <- FindVariableFeatures(object = PBMC.RNA, verbose = FALSE)
PBMC.RNA@assays$RNA@var.features # 480 genes - CXCL13.rna has zero variance

#- Remove uninteresting sources of variation (nUMI as proxy for CDR)
PBMC.RNA <- ScaleData(object = PBMC.RNA, assay = 'RNA', vars.to.regress = 'nCount_RNA')

#- PCA
PBMC.RNA <- RunPCA(object = PBMC.RNA, assay = 'RNA')
ElbowPlot(object = PBMC.RNA, ndims = 50) # 30 PCs
nPC <- 25 #different than valentin's analysis. he set nPC to 30

#- Clustering
PBMC.RNA <- FindNeighbors(object = PBMC.RNA, reduction = 'pca', dims = 1:nPC, k.param = 30, force.recalc = T)
PBMC.RNA <- FindClusters(object = PBMC.RNA, dims.use = 1:nPC, verbose = TRUE, n.start = 100)

#- Visualization (UMAP)
PBMC.RNA <- RunUMAP(object = PBMC.RNA, reduction = 'pca', dims = 1:nPC, min_dist = 0.2, seed.use = 42, n_neighbors = 30, metric = 'correlation')
DimPlot(object = PBMC.RNA, reduction = 'umap', group.by = 'donor', pt.size = 0.5)
DimPlot(object = PBMC.RNA, reduction = 'umap', group.by = 'RNA_snn_res.0.8', pt.size = 0.5) # 19 clusters

###--- UMAP
data_ggplot <- data.table(PBMC.RNA@meta.data, Embeddings(object = PBMC.RNA, reduction = 'umap'))
plot_UMAP <- ggplot(data_ggplot %>% arrange(sample(x = cell.id, replace = FALSE)), aes(x = UMAP_1, y = UMAP_2)) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

#- Batch
plot_UMAP + geom_point(aes(color = as.character(batch)), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE)
#There is no batch (cartridges) effect

#- Sequencing Depth
plot_UMAP + geom_point(aes(color = nCount_RNA), size = 0.2) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + 
  scale_color_continuous(guide = FALSE, low = 'royalblue', high = 'red')
#There is no sequencing effect (except for MPS)

#- Donor
plot_UMAP + geom_point(aes(color = donor), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE)

#- Donor
plot_UMAP + geom_point(aes(color = donor), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ donor)

#- Cell type
plot_UMAP + geom_point(aes(color = celltype), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE)

##-- Marker genes
#- Find markers for every cluster compared to all remaining cells, report only the positive ones
markers.RNA <- FindAllMarkers(object = PBMC.RNA, 
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25, # log-FC
                              test.use = 'MAST',
                              verbose = FALSE,
                              assay = 'RNA',
                              latent.vars = 'nCount_RNA') # nUMI as proxy for CDR
top10RNA <- markers.RNA %>% 
  group_by(cluster) %>%
  top_n(n = 10)

#DoHeatmap(object = PBMC.RNA, features = top10RNA$gene, size = 1, label = F)


# 0 - Naive CD4+ T cells
# 1 - Memory CD4+ T cells
# 2 - NK cells CD56 dim
# 3 - Naive CD8+ T cells
# 4 - Tregs
# 5 - Effector CD8+ and CD4+ T cells
# 6 - Naive B cells
# 7 - CD14+ Monocytes
# 8 - Central Memory CD8+ T cells
# 9 - Effector Memory CD8+ T cells
# 10 - Memory B cells
# 11 - CD14+ Monocytes
# 12 - CD16+ Monocytes
# 13 - Naive CD4+ T cells
# 14 - NK T cells
# 15 - DC (MPS)
# 16 - NK cells CD56 Bright
# 17 - CD4+ T/ B cell doublets? 
# 18 - pDCs

#- chordDiagram
data_ggplot <- data.table(Cluster = paste('Cluster', PBMC.RNA@active.ident, sep = ' '), PBMC.RNA@meta.data)
data_plot <- data.table(table(data_ggplot$donor,
                              data_ggplot$Cluster))
colnames(data_plot) <- c('to', 'from', 'value')
# scales::hue_pal()(3)
# scales::hue_pal()(19)
order.tmp <- c('Donor1', 'Donor3', 'Donor2', paste('Cluster', 0:18))
grid.tmp <- c('Donor1' = '#F8766D', 'Donor3' = '#00BA38', 'Donor2' = '#619CFF', 
              'Cluster 0' = '#F8766D', 'Cluster 1' = '#E9842C', 'Cluster 2' = '#D69100', 
              'Cluster 3' = '#BC9D00', 'Cluster 4' = '#9CA700', 'Cluster 5' = '#6FB000', 
              'Cluster 6' = '#00B813', 'Cluster 7' = '#00BD61', 'Cluster 8' = '#00C08E', 
              'Cluster 9' = '#00C0B4', 'Cluster 10' = '#00BDD4', 'Cluster 11' = '#00B5EE',
              'Cluster 12' = '#00A7FF', 'Cluster 13' = '#7F96FF', 'Cluster 14' = '#BC81FF',
              'Cluster 15' = '#E26EF7', 'Cluster 16' = '#F863DF', 'Cluster 17' = '#FF62BF',
              'Cluster 18' = '#FF6A9A')

chordDiagram(data_plot, 
             order = order.tmp, 
             grid.col = grid.tmp,
             annotationTrack = 'grid',
             annotationTrackHeight = c(0.03, 0.04),
             preAllocateTracks = list(track.height = 0.1))
circos.track(track.index = 1,
             panel.fun = function(x, y){
               circos.text(CELL_META$xcenter,
                           CELL_META$ylim[1],
                           CELL_META$sector.index,
                           facing = 'clockwise',
                           niceFacing = TRUE,
                           cex = 0.5,  
                           adj = c(0, 0.5))}, 
             bg.border = NA)

###--- Markers
PBMC.RNA@active.ident <- plyr::mapvalues(x = PBMC.RNA@active.ident, 
                                         from = 0:18,  
                                         to = c('0 - Naive CD4+ T cells', '1 - Memory CD4+ T cells', '2 - NK cells CD56 dim', '3 - Naive CD8+ T cells', '4 - Tregs', '5 - Effector CD8+ and CD4+ T cells', '6 - Naive B cells', '7 - CD14+ Monocytes', '8 - Central Memory CD8+ T cells', '9 - Effector Memory CD8+ T cells', '10 - Memory B cells', '11 - CD14+ Monocytes', '12 - CD16+ Monocytes', '13 - Naive CD4+ T cells', '14 - NK T cells', '15 - DCs', '16 - NK cells CD56 Bright', '17 - CD4+ T/ B cell doublets?', '18 - pDCs'))

markers.to.plot <- c('CD3D.rna', 'CD3E.rna', 'CD4.rna', 'CD8A.rna', 'CD8B.rna', 'CCR7.rna', 'SELL.rna', 'IL32.rna',  'CCR4.rna', 'GNLY.rna', 'NKG7.rna', 'CD14.rna', 'CTLA4.rna', 'FOXP3.rna', 'CD79A.rna', 'FCER1A.rna', 'LYZ.rna', 'FCGR3A.rna', 'IGHG1-secreted.rna', 'PCNA.rna', 'MCM4-PolyA-1.rna', 'MCM2.rna', 'IRF8.rna', 'IRF7.rna')

DotPlot(object = PBMC.RNA, 
        features = rev(markers.to.plot), 
        dot.scale = 8,
        split.by = 'donor',
        cols = c('red', 'forestgreen', 'blue')) + RotatedAxis() + scale_fill_continuous()

###- Ab expression
#- Setup seurat object
DefaultAssay(object = PBMC.RNA) <- 'Ab'
PBMC.AB <- PBMC.RNA

#- Normalization
PBMC.AB <- NormalizeData(object = PBMC.AB, assay = "Ab", normalization.method = 'CLR')

#- Highly variABle genes
PBMC.AB <- FindVariableFeatures(object = PBMC.AB, verbose = FALSE)
PBMC.AB@assays$Ab@var.features # 42 features

#- Scale data
PBMC.AB <- ScaleData(object = PBMC.AB, assay = 'Ab')

#- PCA
PBMC.AB <- RunPCA(object = PBMC.AB, assay = 'Ab', npcs = 30)
ElbowPlot(object = PBMC.AB, ndims = 30) # 20 PCs
nPC <- 20

#- Clustering
PBMC.AB <- FindNeighbors(object = PBMC.AB, reduction = 'pca', dims = 1:nPC)
PBMC.AB <- FindClusters(object = PBMC.AB, dims.use = 1:nPC, verbose = TRUE)

#- Visualization (UMAP)
PBMC.AB <- RunUMAP(object = PBMC.AB, reduction = 'pca', dims = 1:nPC, min_dist = 0.2, seed.use = 42, n_neighbors = 30, metric = 'correlation')
DimPlot(object = PBMC.AB, reduction = 'umap', group.by = 'donor', pt.size = 0.5)
DimPlot(object = PBMC.AB, reduction = 'umap', group.by = 'Ab_snn_res.0.8', pt.size = 0.5) # 13 clusters

###--- UMAP
data_ggplot <- data.table(PBMC.AB@meta.data, Embeddings(object = PBMC.AB, reduction = 'umap'))
plot_UMAP <- ggplot(data_ggplot %>% arrange(sample(x = cell.id, replace = FALSE)), aes(x = UMAP_1, y = UMAP_2)) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

#- Batch
plot_UMAP + geom_point(aes(color = as.character(batch)), size = 0.2) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE)
#There is no batch (cartridges) effect. 

#- Sequencing depth
plot_UMAP + geom_point(aes(color = nCount_Ab), size = 0.2) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + 
  scale_color_continuous(guide = FALSE, low = 'royalblue', high = 'red')
#There is no sequencing effect (except for MPS).  

#- Donor
plot_UMAP + geom_point(aes(color = donor), size = 0.2) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE)

#- Donor
plot_UMAP + geom_point(aes(color = donor), size = 0.2) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ donor)
#We can observe some clusters composed/driven by donors.  

#- Cell type
plot_UMAP + geom_point(aes(color = celltype), size = 0.2) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE)

##-- Marker genes
#- Find markers for every cluster compared to all remaining cells, report only the positive ones
markers.AB <- FindAllMarkers(object = PBMC.AB, 
                             only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.25, # log-FC
                             test.use = 'MAST',
                             verbose = FALSE,
                             assay = 'Ab') 
markers.AB %>% 
  group_by(cluster) %>%
  top_n(n = 20)
# 0 - Effector Memory CD4+ T
# 1 - Naive CD4+ T
# 2 - NK
# 3 - Naive CD4+ T
# 4 - Central Memory CD4+ T
# 5 - Effector CD8+ T cells and NKT cells
# 6 - Naive CD8+ T
# 7 - Memo CD8+ T
# 8 - B
# 9 - CD14+ Mono
# 10 - CD16+ Mono
# 11 - Naive/CM T
# 12 - Effector Memory CD4+/CD8+ T

RidgePlot(object = PBMC.AB, features = c('CD3.CD3E.ab', 'CD4.CD4.ab', 'CD8.CD8A.ab', 'CD45RA.PTPRC.ab', 'CD45RO.PTPRC.ab', 'CD197.CCR7.ab'), assay = 'Ab')

FeaturePlot(object = PBMC.AB, 
            features = c("CD3.CD3E.ab", "CD4.CD4.ab", "CD8.CD8A.ab", "CD56.NCAM1.ab", "CD14.CD14.ab", "CD19.CD19.ab", "CD197.CCR7.ab", "CD16.FCGR3A.ab", 'CD45RA.PTPRC.ab', 'CD45RO.PTPRC.ab', 'CD11c.ITGAX.ab', "CD27.CD27.ab", "CD28.CD28.ab"))


data_ggplot$cluster.id <- plyr::mapvalues(x = data_ggplot$Ab_snn_res.0.8, 
                                          from = 0:12,  
                                          to = c('Effector Memory CD4+ T', 'Naive CD4+ T', 'NK', 'Naive CD4+ T', 'Central Memory CD4+ T', 'Effector CD8+ T cells and NKT cells', 'Naive CD8+ T', 'Naive CD8+ T', 'B', 'CD14+ Mono', 'CD16+ Mono', 'Naive/CM T', 'Effector Memory CD4+/CD8+ T'))

#- chordDiagram
data_ggplot <- data.table(Cluster = paste('Cluster', PBMC.AB@active.ident, sep = ' '), PBMC.AB@meta.data)
data_plot <- data.table(table(data_ggplot$donor,
                              data_ggplot$Cluster))
colnames(data_plot) <- c('to', 'from', 'value')
# scales::hue_pal()(3)
# scales::hue_pal()(13)
order.tmp <- c('Donor1', 'Donor3', 'Donor2', paste('Cluster', 0:12))
grid.tmp <- c('Donor1' = '#F8766D', 'Donor3' = '#00BA38', 'Donor2' = '#619CFF', 
              'Cluster 0' = '#F8766D', 'Cluster 1' = '#E18A00', 'Cluster 2' = '#BE9C00', 
              'Cluster 3' = '#8CAB00', 'Cluster 4' = '#24B700', 'Cluster 5' = '#00BE70', 
              'Cluster 6' = '#00C1AB', 'Cluster 7' = '#00BBDA', 'Cluster 8' = '#00ACFC', 
              'Cluster 9' = '#8B93FF', 'Cluster 10' = '#D575FE', 'Cluster 11' = '#F962DD',
              'Cluster 12' = '#FF65AC')

chordDiagram(data_plot, 
             order = order.tmp, 
             grid.col = grid.tmp,
             annotationTrack = 'grid',
             annotationTrackHeight = c(0.03, 0.04),
             preAllocateTracks = list(track.height = 0.1))
circos.track(track.index = 1,
             panel.fun = function(x, y){
               circos.text(CELL_META$xcenter,
                           CELL_META$ylim[1],
                           CELL_META$sector.index,
                           facing = 'clockwise',
                           niceFacing = TRUE,
                           cex = 0.5,  
                           adj = c(0, 0.5))}, 
             bg.border = NA)

#- Protein markers
#plot <- FeaturePlot(object = PBMC.AB, 
#                    features = rownames(PBMC.AB@assays$Ab), 
#                    min.cutoff = 'q05', 
#                    max.cutoff = 'q95', 
#                    ncol = 6,
#                    cols = c('lightgrey', 'blue'),
#                    pt.size = 0.5, 
#                    reduction = 'umap')
#plot + theme(axis.text = element_blank(), axis.ticks = element_blank()) 

####- Analysis of T cell subsets

markers <- rownames(PBMC.RNA)
#RidgePlot(object = PBMC.AB, features = markers, group.by = 'donor', assay = 'Ab', ncol = 6)

#manually use `scamp` to get our cut-offs for T cell markers (such as CD3+CD8+CD4-, etc). 

###--- SCAMP - get cut-offs within each donor
PBMC.RNA@meta.data$CD3.CD3E_Ab <- PBMC.RNA@assays$Ab@data['CD3.CD3E.ab', ]
PBMC.RNA@meta.data$CD4.CD4_Ab <- PBMC.RNA@assays$Ab@data['CD4.CD4.ab', ]
PBMC.RNA@meta.data$CD8.CD8A_Ab <- PBMC.RNA@assays$Ab@data['CD8.CD8A.ab', ]

PBMC.RNA.list <- SplitObject(object = PBMC.RNA, split.by = 'donor') #I had to change to PBMC.AB, PBMC.RNA had error that 'Ab' was not in object

##- Donor1
PBMC.RNA.tmp <- PBMC.RNA.list$Donor1
#- CD3+
x <- PBMC.RNA.tmp@meta.data$CD3.CD3E_Ab
cutoff.Donor1.CD3 <- tsGates(xVec = sort(addNoiseToDataVector(sort(x), 10, 123)), modePrior = 0)[2]
plot(density(x)); abline(v = cutoff.Donor1.CD3)
PBMC.RNA.tmp@meta.data$phenotype.CD3 <- ifelse(PBMC.RNA.tmp@meta.data$CD3.CD3E_Ab > cutoff.Donor1.CD3, 'CD3+', 'CD3-')

#- CD4+
x <- PBMC.RNA.tmp@meta.data$CD4.CD4_Ab[which(PBMC.RNA.tmp@meta.data$phenotype.CD3 == 'CD3+')]
cutoff.Donor1.CD4 <- tsGates(xVec = sort(addNoiseToDataVector(sort(x), 10, 123)), modePrior = 0)[2]
plot(density(x)); abline(v = cutoff.Donor1.CD4)
PBMC.RNA.tmp@meta.data$phenotype.CD4 <- ifelse(PBMC.RNA.tmp@meta.data$CD4.CD4_Ab > cutoff.Donor1.CD4, 'CD4+', 'CD4-')
PBMC.RNA.tmp@meta.data$phenotype.CD4 <- paste0(PBMC.RNA.tmp@meta.data$phenotype.CD3, PBMC.RNA.tmp@meta.data$phenotype.CD4)

#- CD8+
x <- PBMC.RNA.tmp@meta.data$CD8.CD8A_Ab[which(PBMC.RNA.tmp@meta.data$phenotype.CD3 == 'CD3+')]
cutoff.Donor1.CD8 <- tsGates(xVec = sort(addNoiseToDataVector(sort(x), 10, 123)), modePrior = 0)[2]
plot(density(x)); abline(v = cutoff.Donor1.CD8)
PBMC.RNA.tmp@meta.data$phenotype.CD8 <- ifelse(PBMC.RNA.tmp@meta.data$CD8.CD8A_Ab > cutoff.Donor1.CD8, 'CD8+', 'CD8-')
PBMC.RNA.tmp@meta.data$phenotype.CD8 <- paste0(PBMC.RNA.tmp@meta.data$phenotype.CD3, PBMC.RNA.tmp@meta.data$phenotype.CD8)

output.Donor1 <- PBMC.RNA.tmp@meta.data

##- Donor3
PBMC.RNA.tmp <- PBMC.RNA.list$Donor3
#- CD3+
x <- PBMC.RNA.tmp@meta.data$CD3.CD3E_Ab
cutoff.Donor3.CD3 <- tsGates(xVec = sort(addNoiseToDataVector(sort(x), 10, 123)), modePrior = 0)[2]
plot(density(x)); abline(v = cutoff.Donor3.CD3)
PBMC.RNA.tmp@meta.data$phenotype.CD3 <- ifelse(PBMC.RNA.tmp@meta.data$CD3.CD3E_Ab > cutoff.Donor3.CD3, 'CD3+', 'CD3-')

#- CD4+
x <- PBMC.RNA.tmp@meta.data$CD4.CD4_Ab[which(PBMC.RNA.tmp@meta.data$phenotype.CD3 == 'CD3+')]
cutoff.Donor3.CD4 <- tsGates(xVec = sort(addNoiseToDataVector(sort(x), 10, 123)), modePrior = 0)[2]
plot(density(x)); abline(v = cutoff.Donor3.CD4)
PBMC.RNA.tmp@meta.data$phenotype.CD4 <- ifelse(PBMC.RNA.tmp@meta.data$CD4.CD4_Ab > cutoff.Donor3.CD4, 'CD4+', 'CD4-')
PBMC.RNA.tmp@meta.data$phenotype.CD4 <- paste0(PBMC.RNA.tmp@meta.data$phenotype.CD3, PBMC.RNA.tmp@meta.data$phenotype.CD4)

#- CD8+
x <- PBMC.RNA.tmp@meta.data$CD8.CD8A_Ab[which(PBMC.RNA.tmp@meta.data$phenotype.CD3 == 'CD3+')]
cutoff.Donor3.CD8 <- tsGates(xVec = sort(addNoiseToDataVector(sort(x), 10, 123)), modePrior = 0)[2]
plot(density(x)); abline(v = cutoff.Donor3.CD8)
PBMC.RNA.tmp@meta.data$phenotype.CD8 <- ifelse(PBMC.RNA.tmp@meta.data$CD8.CD8A_Ab > cutoff.Donor3.CD8, 'CD8+', 'CD8-')
PBMC.RNA.tmp@meta.data$phenotype.CD8 <- paste0(PBMC.RNA.tmp@meta.data$phenotype.CD3, PBMC.RNA.tmp@meta.data$phenotype.CD8)

output.Donor3 <- PBMC.RNA.tmp@meta.data

##- Donor2
PBMC.RNA.tmp <- PBMC.RNA.list$Donor2
#- CD3+
x <- PBMC.RNA.tmp@meta.data$CD3.CD3E_Ab
cutoff.Donor2.CD3 <- tsGates(xVec = sort(addNoiseToDataVector(sort(x), 10, 123)), modePrior = 0)[2]
plot(density(x)); abline(v = cutoff.Donor2.CD3)
PBMC.RNA.tmp@meta.data$phenotype.CD3 <- ifelse(PBMC.RNA.tmp@meta.data$CD3.CD3E_Ab > cutoff.Donor2.CD3, 'CD3+', 'CD3-')

#- CD4+
x <- PBMC.RNA.tmp@meta.data$CD4.CD4_Ab[which(PBMC.RNA.tmp@meta.data$phenotype.CD3 == 'CD3+')]
cutoff.Donor2.CD4 <- tsGates(xVec = sort(addNoiseToDataVector(sort(x), 10, 123)), modePrior = 0)[2]
plot(density(x)); abline(v = cutoff.Donor2.CD4)
PBMC.RNA.tmp@meta.data$phenotype.CD4 <- ifelse(PBMC.RNA.tmp@meta.data$CD4.CD4_Ab > cutoff.Donor2.CD4, 'CD4+', 'CD4-')
PBMC.RNA.tmp@meta.data$phenotype.CD4 <- paste0(PBMC.RNA.tmp@meta.data$phenotype.CD3, PBMC.RNA.tmp@meta.data$phenotype.CD4)

#- CD8+
x <- PBMC.RNA.tmp@meta.data$CD8.CD8A_Ab[which(PBMC.RNA.tmp@meta.data$phenotype.CD3 == 'CD3+')]
cutoff.Donor2.CD8 <- tsGates(xVec = sort(addNoiseToDataVector(sort(x), 10, 123)), modePrior = 0)[2]
plot(density(x)); abline(v = cutoff.Donor2.CD8)
PBMC.RNA.tmp@meta.data$phenotype.CD8 <- ifelse(PBMC.RNA.tmp@meta.data$CD8.CD8A_Ab > cutoff.Donor2.CD8, 'CD8+', 'CD8-')
PBMC.RNA.tmp@meta.data$phenotype.CD8 <- paste0(PBMC.RNA.tmp@meta.data$phenotype.CD3, PBMC.RNA.tmp@meta.data$phenotype.CD8)

output.Donor2 <- PBMC.RNA.tmp@meta.data

###--- cutoff
CD3.cutoff <- c(cutoff.Donor1.CD3, cutoff.Donor3.CD3, cutoff.Donor2.CD3)
CD4.cutoff <- c(cutoff.Donor1.CD4, cutoff.Donor3.CD4, cutoff.Donor2.CD4)
CD8.cutoff <- c(cutoff.Donor1.CD8, cutoff.Donor3.CD8, cutoff.Donor2.CD8)

###--- meta.data
meta.data <- rbind(output.Donor1, output.Donor3, output.Donor2)
meta.data <- meta.data[rownames(PBMC.RNA@meta.data), ]
PBMC.RNA@meta.data$phenotype.CD3 <- meta.data$phenotype.CD3
PBMC.RNA@meta.data$phenotype.CD4 <- meta.data$phenotype.CD4
PBMC.RNA@meta.data$phenotype.CD8 <- meta.data$phenotype.CD8
PBMC.RNA@meta.data$phenotype.CD4 <- str_remove(string = PBMC.RNA@meta.data$phenotype.CD4, pattern = 'CD3[+/-]')
PBMC.RNA@meta.data$phenotype.CD8 <- str_remove(string = PBMC.RNA@meta.data$phenotype.CD8, pattern = 'CD3[+/-]')
PBMC.RNA@meta.data$T.phenotype <- paste0(PBMC.RNA@meta.data$phenotype.CD3, 
                                         PBMC.RNA@meta.data$phenotype.CD4, 
                                         PBMC.RNA@meta.data$phenotype.CD8)
PBMC.RNA@meta.data$T.phenotype[str_detect(string = PBMC.RNA@meta.data$T.phenotype, pattern = 'CD3[-]')] <- 'Non-T'
PBMC.RNA@meta.data$T.phenotype[str_detect(string = PBMC.RNA@meta.data$T.phenotype, pattern = 'CD3[+]CD4[+]CD8[+]')] <- 'CD3+CD4+CD8+'
PBMC.RNA@meta.data$T.phenotype[str_detect(string = PBMC.RNA@meta.data$T.phenotype, pattern = 'CD3[+]CD4[-]CD8[-]')] <- 'CD3+CD4-CD8-'
View(PBMC.RNA@meta.data)



###--- RidgePlot

#- CD3+
plot.CD3 <- RidgePlot(object = PBMC.RNA,  features = 'CD3.CD3E_Ab', group.by = 'donor') + 
  NoLegend() +
  geom_vline(xintercept = CD3.cutoff, color = c('red', 'forestgreen', 'royalblue')) +
  theme(axis.title.y = element_blank(), title = element_blank())

#- CD4+
PBMC.RNA.tmp <- PBMC.RNA
PBMC.RNA.tmp <- SetIdent(object = PBMC.RNA.tmp, value = 'phenotype.CD3')
PBMC.RNA.tmp <- SubsetData(object = PBMC.RNA.tmp,  ident.use = 'CD3+')
plot.CD4 <- RidgePlot(object = PBMC.RNA.tmp,  features = 'CD4.CD4_Ab', group.by = 'donor') + 
  NoLegend() +
  geom_vline(xintercept = CD4.cutoff, color = c('red', 'forestgreen', 'royalblue')) +
  theme(axis.title.y = element_blank(), title = element_blank())

#- CD8+
PBMC.RNA.tmp <- PBMC.RNA
PBMC.RNA.tmp <- SetIdent(object = PBMC.RNA.tmp, value = 'phenotype.CD3')
PBMC.RNA.tmp <- SubsetData(object = PBMC.RNA.tmp,  ident.use = 'CD3+')
plot.CD8 <- RidgePlot(object = PBMC.RNA.tmp,  features = 'CD8.CD8A_Ab', group.by = 'donor') + 
  NoLegend() +
  geom_vline(xintercept = CD8.cutoff, color = c('red', 'forestgreen', 'royalblue')) +
  theme(axis.title.y = element_blank(), title = element_blank())


cowplot::plot_grid(plot.CD3, plot.CD8, 
                   plot.CD4,  
                   cols = c("#e41a1c", "#377eb8", "#4daf4a"),
                   ncol = 3, 
                   labels = c('CD3+', 'CD3+CD4-CD8+', 'CD3+CD4+CD8-'))


#- Subsetting

PBMC.RNA.AllCD8 <- subset(x = PBMC.RNA, phenotype.tmp == 'CD3+CD8+CD4-')
PBMC.RNA.AllCD8 # 3502 cells


#- Normalization
DefaultAssay(object = PBMC.RNA.AllCD8) <- 'RNA'
PBMC.RNA.AllCD8 <- NormalizeData(object = PBMC.RNA.AllCD8, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = 10000)

#- Highly variable genes
PBMC.RNA.AllCD8 <- FindVariableFeatures(object = PBMC.RNA.AllCD8, verbose = FALSE)
PBMC.RNA.AllCD8@assays$RNA@var.features # 432 genes 

#- Remove uninteresting sources of variation (nUMI as proxy for CDR)
PBMC.RNA.AllCD8 <- ScaleData(object = PBMC.RNA.AllCD8, assay = 'RNA', vars.to.regress = 'nCount_RNA')

#- PCA
PBMC.RNA.AllCD8 <- RunPCA(object = PBMC.RNA.AllCD8, assay = 'RNA')
ElbowPlot(object = PBMC.RNA.AllCD8, ndims = 50) # 30 PCs
nPC <- 30

#- Clustering
PBMC.RNA.AllCD8 <- FindNeighbors(object = PBMC.RNA.AllCD8, reduction = 'pca', dims = 1:nPC, k.param = 30, force.recalc = T)
PBMC.RNA.AllCD8 <- FindClusters(object = PBMC.RNA.AllCD8, dims.use = 1:nPC, verbose = TRUE, n.start = 100, resolution = 0.8)
PBMC.RNA.AllCD8 <- RunUMAP(object = PBMC.RNA.AllCD8, reduction = 'pca', dims = 1:nPC, min_dist = 0.2, seed.use = 42, n_neighbors = 30, metric = 'correlation')
DimPlot(object = PBMC.RNA.AllCD8, reduction = 'umap', group.by = 'RNA_snn_res.0.8', pt.size = 1, cols = brewer.pal(7, "Paired")) # 7 clusters

#remove b and mono contamination
PBMC.RNA.AllCD8 <- subset(PBMC.RNA.AllCD8, ident = c("0", "1", "2", "3", "4"))
PBMC.RNA.AllCD8 <- FindNeighbors(object = PBMC.RNA.AllCD8, reduction = 'pca', dims = 1:nPC, k.param = 30, force.recalc = T)
PBMC.RNA.AllCD8 <- FindClusters(object = PBMC.RNA.AllCD8, dims.use = 1:nPC, verbose = TRUE, n.start = 100, resolution = 0.8)

#- Visualization (UMAP)
PBMC.RNA.AllCD8 <- RunUMAP(object = PBMC.RNA.AllCD8, reduction = 'pca', dims = 1:nPC, min_dist = 0.2, seed.use = 42, n_neighbors = 30, metric = 'correlation')
DimPlot(object = PBMC.RNA.AllCD8, reduction = 'umap', group.by = 'donor', pt.size = 1, cols = c("#e41a1c", "#377eb8", "#4daf4a"))
DimPlot(object = PBMC.RNA.AllCD8, reduction = 'umap', group.by = 'RNA_snn_res.0.8', pt.size = 1, cols = brewer.pal(5, "Paired")) # 5 clusters


markers.AllCD8 <- FindAllMarkers(PBMC.RNA.AllCD8, 
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.5, # log-FC
                                 test.use = 'MAST',
                                 verbose = FALSE,
                                 assay = 'RNA',
                                 latent.vars = 'nCount_RNA') # nUMI as proxy for CDR

AllCD8top5RNA <- markers.AllCD8 %>% 
  group_by(cluster) %>%
  top_n(n = 5, avg_logFC)

VlnPlot(PBMC.RNA.AllCD8, features = markers.AllCD8$gene, ncol = 5, pt.size = 0, cols = brewer.pal(5, "Paired"))

VlnPlot(PBMC.RNA.AllCD8, features = c("SELL.rna", "CCR7.rna", "GZMK.rna", "HLA-DPA1-PolyA-1.rna", "GZMB.rna", "KLRF1.rna", "IL4R.rna", "LGALS1.rna", "IL18RAP.rna"), 
        ncol = 3, pt.size = 0, cols = brewer.pal(5, "Paired"))

#Ab

DefaultAssay(object = PBMC.RNA.AllCD8) <- 'Ab'

#- Highly variable genes
PBMC.RNA.AllCD8 <- FindVariableFeatures(object = PBMC.RNA.AllCD8, verbose = FALSE)
PBMC.RNA.AllCD8@assays$Ab@var.features # 42 genes 

#- Remove uninteresting sources of variation (nUMI as proxy for CDR)
PBMC.RNA.AllCD8 <- ScaleData(object = PBMC.RNA.AllCD8, assay = 'Ab', vars.to.regress = 'nCount_Ab')

#- PCA
PBMC.RNA.AllCD8 <- RunPCA(object = PBMC.RNA.AllCD8, assay = 'Ab')
ElbowPlot(object = PBMC.RNA.AllCD8, ndims = 50) # 30 PCs
nPC <- 40

#- Clustering
PBMC.RNA.AllCD8 <- FindNeighbors(object = PBMC.RNA.AllCD8, reduction = 'pca', dims = 1:nPC, k.param = 30, force.recalc = T)
PBMC.RNA.AllCD8 <- FindClusters(object = PBMC.RNA.AllCD8, dims.use = 1:nPC, verbose = TRUE, n.start = 100, resolution = 0.8)

#- Visualization 
# (UMAP)
PBMC.RNA.AllCD8 <- RunUMAP(object = PBMC.RNA.AllCD8, reduction = 'pca', dims = 1:nPC, min_dist = 0.2, seed.use = 42, n_neighbors = 30, metric = 'correlation')
DimPlot(object = PBMC.RNA.AllCD8, reduction = 'umap', group.by = 'donor', pt.size = 1, cols = c("#e41a1c", "#377eb8", "#4daf4a"))
DimPlot(object = PBMC.RNA.AllCD8, reduction = 'umap', group.by = 'Ab_snn_res.0.8', pt.size = 1, cols = brewer.pal(8, "Paired")) # 8 clusters

#RidgePlot

RidgePlot(PBMC.RNA.AllCD8, features = c("CD45RA.PTPRC.ab", "CD27.CD27.ab", "CD28.CD28.ab", "CD25.IL2RA.ab", "CD127.IL7R.ab", "CD161.KLRB1.ab"), cols = brewer.pal(5, "Paired"))

#EBV

DefaultAssay(object = PBMC.RNA.AllCD8) <- 'RNA'
PBMC.RNA.AllCD8 -> PBMC.RNA.EBVCD8
PBMC.RNA.EBVCD8@meta.data$phenotype.tmp <- paste0(PBMC.RNA.AllCD8@meta.data$celltype)

PBMC.RNA.EBVCD8 <- subset(x = PBMC.RNA.EBVCD8, phenotype.tmp == 'EBV')

DimPlot(object = PBMC.RNA.AllCD8, reduction = 'umap', group.by = 'celltype', pt.size = 1, cols = c('#bfbfbf', '#a50f15')) # 5 clusters

DimPlot(object = PBMC.RNA.EBVCD8, reduction = 'umap', group.by = 'donor', pt.size = 1, cols = c("#e41a1c", "#377eb8", "#4daf4a"))



PBMC.RNA.AllCD8 <- SetIdent(object = PBMC.RNA.AllCD8, value = 'RNA_snn_res.0.8')
PBMC.RNA.CD8.1 <- SubsetData(object = PBMC.RNA.AllCD8, ident.use = '1')
PBMC.RNA.CD8.1 <- SetIdent(object = PBMC.RNA.CD8.1, value = 'celltype')


markers.1 <- FindAllMarkers(PBMC.RNA.CD8.1, 
                            only.pos = TRUE,
                            min.pct = 0.25,
                            logfc.threshold = 0.5, # log-FC
                            test.use = 'MAST',
                            verbose = FALSE,
                            assay = 'Ab',
                            latent.vars = 'nCount_RNA') # nUMI as proxy for CDR
markers.1top5 <- markers.1 %>% 
  group_by(cluster) %>%
  top_n(n = 5, avg_logFC)

VlnPlot(PBMC.RNA.CD8.1, features = markers.1top5$gene, ncol = 2, pt.size = 0, cols = c("#bdbdbd", "#a50f15"))

VlnPlot(PBMC.RNA.CD8.1, features = c("GNLY.rna", "YBX3.rna"), ncol = 2, pt.size = 0, cols = c("#bdbdbd", "#a50f15"))
