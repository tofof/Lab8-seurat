library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(patchwork)
library(sctransform)

# PREAMBLE
#----------
# Following tutorials at
#   https://broadinstitute.github.io/2020_scWorkshop/data-wrangling-scrnaseq.html
#   https://broadinstitute.github.io/KrumlovSingleCellWorkshop2020/data-wrangling-scrnaseq-1.html
#   https://satijalab.org/seurat/archive/v3.0/pbmc3k_tutorial.html
# Public dataset of Non-Small Cell Lung Cancer Cells (NSCLC) freely available from 10X Genomics
#   curl -O https://cf.10xgenomics.com/samples/cell-vdj/2.2.0/vdj_v1_hs_nsclc_5gex/vdj_v1_hs_nsclc_5gex_filtered_gene_bc_matrices.tar.gz
#   tar -xvzf vdj_v1_hs_nsclc_5gex_filtered_gene_bc_matrices.tar.gz
# Tirosh et al housekeeping gene list from
#   https://github.com/Michorlab/tnbc_scrnaseq/blob/master/data/housekeepers.txt

# DATA INITIALIZATION
#---------------------

datadir <- "./NSCLC/"
counts_matrix_filename <- paste0(datadir,"/filtered_gene_bc_matrices/GRCh38/") # Seurat function reads genes.tsv, matrix.mtx, barcodes.tsv

counts <- Read10X(data.dir = counts_matrix_filename)
countsPerCell <- Matrix::colSums(counts)
countsPerGene <- Matrix::rowSums(counts)
genesPerCell <- Matrix::colSums(counts > 0) # count gene only if it has non-zero reads mapped
cellsPerGene <- Matrix::rowSums(counts > 0) # count cells only where the gene is expressed

hist(log10(countsPerCell + 1), main = 'Counts per Cell', col = 'wheat')
hist(log10(genesPerCell + 1), main = 'Genes per Cell', col = 'wheat')
hist(log10(countsPerGene + 1), main = 'Counts per Gene', col = 'wheat')
hist(log10(cellsPerGene + 1), main = 'Cells per Gene', col = 'wheat')
plot(countsPerCell, genesPerCell, log = 'xy', col = 'dodgerblue4')
plot(cellsPerGene,countsPerGene, log = 'xy', col = 'dodgerblue4')
plot(sort(genesPerCell), log = 'y', xlab = 'Cell', ylab = 'Complexity', main = 'Genes per Cell (ordered)')

# Seurat data objects contain raw and computed data in 'slots' within the object.
# The Assay class stores single cell data. Typically, scRNA-seq will only yield one Assay object: $RNA.
# That assay presents multiple transformations of the data:
#   @counts - original count matrix, unnormalized raw counts or TPMs
#   @data -
#     normalized expression data [source: https://www.rdocumentation.org/packages/Seurat/versions/3.0.1/topics/Assay-class]
#     log-normalized data [source: https://github.com/satijalab/seurat/issues/3711]
#     normalized data [source: https://satijalab.org/seurat/archive/v3.0/pbmc3k_tutorial.html "Normalized values are stored in pbmc[["RNA"]]@data."]
#     WRONG WRONG WRONG: raw non-transformed, non-log-normalized counts [source: https://broadinstitute.github.io/KrumlovSingleCellWorkshop2020/data-wrangling-scrnaseq-1.html]
#   @scale.data - scaled expression matrix, used for dimensional reduction and heatmap visualization
seuratData <- CreateSeuratObject(counts, min.cells = 3, min.features = 350, project = "10x_NSCLC")
seuratData
seuratData@assays$RNA@counts[1:10,1:10]
genesPerCell_sd <- Matrix::colSums(seuratData@assays$RNA@counts)
plot(sort(genesPerCell_sd), log = 'y', xlab = 'Cell', ylab = 'Complexity', main = 'Genes per Cell (ordered, cells > 3, genes > 350)')


# QUALITY CONTROL & NORMALIZATION
#---------------------------------

# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell. We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since <WRONG> this represents
# non-transformed and non-log-normalized counts </WRONG>.

# here's a horrible way to do it:
#   mito.genes <- grep(pattern = "^MT-", x = rownames(seuratData@assays$RNA@data), value = TRUE) # MT-genes are prefixed as such
#   percent.mito <- Matrix::colSums(seuratData@assays$RNA@data[mito.genes, ])/Matrix::colSums(seuratData@assays$RNA@data)
#   seurat <- AddMetaData(object = seurat, metadata = percent.mito, col.name = "percent.mito")
# instead just do this, which automatically(?) stashes into meta.data:
seuratData[["percent.mt"]] <- PercentageFeatureSet(seuratData, pattern = "^MT-")
VlnPlot(seuratData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# housekeeping genes should be generally expressed in all cells
hkgenes <- as.vector(read.table(paste0(datadir, "tirosh_housekeepers.txt"))$V1) # don't add `,skip=2` to read.table call, because our file has no header
hkgeneRows <- which(toupper(rownames(seuratData@assays$RNA@data)) %in% toupper(hkgenes)) # get row indices in @data that match
seuratData[["n.exp.hkgenes"]] <- Matrix::colSums(seuratData@assays$RNA@data[hkgeneRows, ] > 0)
VlnPlot(seuratData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "n.exp.hkgenes"), ncol = 4)
# VlnPlot(seuratData, features = c("nFeature_RNA"), group.by = c("orig.ident")) # grouped by identity (currently identical to prior plots, because identity information will only change once we do clustering)

plotMTpercentVsUMIs <- FeatureScatter(seuratData, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotGenesVsUMIs <- FeatureScatter(seuratData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotMTpercentVsUMIs + plotGenesVsUMIs

# subset thresholds determined by inspecting violin plots
# alternative to explicit subset is to use FilterCells
#   seuratData <- FilterCells(seuratData, subset.names = c("nFeature_RNA", "percent.mt"),
#                             low.thresholds = c(350, -Inf), high.thresholds = c(6000, 12.5))
seuratData <- subset(seuratData, subset = nFeature_RNA > 350 & nFeature_RNA < 4000 &
                                          # nCount_RNA < 60000) &
                                          percent.mt < 15 &
                                          n.exp.hkgenes >= 55)

seuratData
VlnPlot(seuratData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "n.exp.hkgenes"), ncol = 4)
plotMTpercentVsUMIs <- FeatureScatter(seuratData, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotGenesVsUMIs <- FeatureScatter(seuratData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotMTpercentVsUMIs + plotGenesVsUMIs

# Variable regression using sctransform
#   single command replaces NormalizeData, ScaleData, and FindVariableFeatures
#   see further discussion below. Despite the loss of JackStraw, this should probably be the preferred workflow.
seuratDataSC <- SCTransform(seuratData, vars.to.regress = "percent.mt", verbose = TRUE)
rm(seuratData)


# DIMENSIONAL REDUCTION / PRINCIPAL COMPONENT ANALYSIS
#------------------------------------------------------

seuratDataSC <- RunPCA(seuratDataSC, features = VariableFeatures(seuratDataSC), ndims.print = 1:5, nfeatures.print = 5) # equivalently, features = seuratData@assays$RNA@var.features but note that that means specifying assay RNA, whereas if you do SCTransform it's assay SCT. Avoid confusion with the VariableFeatures call instead.
VizDimLoadings(seuratDataSC, dims = 1:2, reduction = "pca")
DimPlot(seuratDataSC, reduction = "pca")

seuratDataSC <- RunICA(seuratDataSC, features = VariableFeatures(seuratDataSC), ndims.print = 1:5, nfeatures.print = 5)
VizDimLoadings(seuratDataSC, dims = 1:2, reduction = "ica")
DimPlot(seuratDataSC, reduction = "ica")


# INITIAL HEATMAP
#-----------------

DimHeatmap(seuratDataSC, dims = 1:30, cells = 500, reduction = "pca", balanced = TRUE)


# DIMENSIONALITY DETERMINATION
#------------------------------
# First step in clustering is to determine the 'true' dimensionality of the dataset.
# Determine dimensionality of the dataset using the JackStraw procedure from Macosko et al (http://www.cell.com/abstract/S0092-8674(15)00549-8)
# JackStraw cannot be used on SCT-transformed data - JackStraw assumes each gene has equal variance, but SCT normalization weights gene variance by biological heterogeneity
# Contrast with less-computationally-demanding ElbowPlot

ElbowPlot(seuratDataSC, ndims = 50, reduction = "pca")


# CLUSTERING
#------------

# Seurat v3 has a FindNeighbors() that constructs a SNN/KNN graph (based on euclidian dist in PCA space), with edge weights refined by Jaccard similarity (shared overlap in local neighborhood)
seuratDataSC <- FindNeighbors(seuratDataSC, dims = 1:50) # where dimensionality used (50) is determined in the step above

# Seurat v3 also has a FindClusters() that does modularity optimization techniques (Louvain by default, or SLM)
seuratDataSC <- FindClusters(seuratDataSC,resolution = 1.6) # where resolution affects granularity of downstream clusters, increasing->more clusters.  Typically use 0.4-1.2 for sc datasets of around 3K cells, increasing as necessary for larger datasets. My use of 1.6 is completely uninformed.
head(Idents(seuratData),5) # look at cluster IDs of the first 5 cells

# Seurat v3 offers UMAP and tSNE to visualize clustering
seuratDataSC <- RunUMAP(seuratDataSC, dims = 1:50)
DimPlot(seuratDataSC, reduction = "umap", label = TRUE) + NoLegend()
seuratDataSC <- RunTSNE(seuratDataSC, reduction.use = "pca", dims.use = 1:50, perplexity = 10)
DimPlot(seuratDataSC, reduction = "tsne", label = TRUE) + NoLegend()


# DATA EXPORT
#-------------

saveRDS(seuratDataSC, file = paste0(datadir,"output/seuratData_tutorialSC.rds")) # Can be loaded back in without rerunning computationally intense steps
seuratDataSC <- readRDS(file = paste0(datadir,"output/seuratData_tutorialSC.rds")) # load


# DIFFERENTIALLY EXPRESSED FEATURES / CLUSTER BIOMARKERS
#--------------------------------------------------------

cluster1.markers <- FindMarkers(seuratDataSC, ident.1 = 1, min.pct = 0.2) # cluster to ID is ident.1 (to do all clusters, use FindAllMarkers instead). Min.pct is detection threshold, thresh.test is to require differential expression amount between two groups. Finally, max.cells.per.ident can speed computation for a small loss in power
head(cluster1.markers)

cluster2.markers <- FindMarkers(seuratDataSC, ident.1 = 2, ident2 = c(8,11)) # find all markers distinguishing cluster 2 from clusters 8 and 11
head(cluster2.markers)

cluster414.markers <- FindMarkers(seuratDataSC, ident.1 = c(4,14), only.pos = TRUE)
head(cluster414.markers)

nsclc.markers <- FindAllMarkers(seuratDataSC, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
saveRDS(nsclc.markers, file = paste0(datadir,"output/nsclc.markersSC.rds"))
nsclc.markers <- readRDS(file = paste0(datadir,"output/nsclc.markersSC.rds"))

nsclc.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
nsclc.markers %>% filter(cluster == 1 | cluster == 2 | cluster == 3) %>% filter(n() > 2 & avg_log2FC > 1.4)
top20 <- nsclc.markers %>% top_n(20, avg_log2FC)
top20pvalue <- nsclc.markers %>% top_n(-20, p_val)
top20pvalue
top5bycluster <- nsclc.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

VlnPlot(seuratDataSC, features = c("MS4A1", "JCHAIN"))
VlnPlot(seuratDataSC, features = c("NKG7", "TYROBP"), log = TRUE)
VlnPlot(seuratDataSC, features = c("SERINC2", "CD3E"))
FeaturePlot(seuratDataSC, features = c("MS4A1", "TYROBP", "SERINC2", "CD3E"))
DoHeatmap(seuratDataSC, features = top20$gene, label = TRUE, size = 3)
DoHeatmap(seuratDataSC, features = top20pvalue$gene, label = TRUE, size = 3)
DoHeatmap(seuratDataSC, features = top5bycluster$gene, label = TRUE, size = 3)

# RENAMING CLUSTERS
#-------------------

# if you know from looking cluster membership and recognize cell types from known markers, you can manually specify:
#   new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
#   names(new.cluster.ids) <- levels(seuratDataSC)
#   seuratDataSC <- RenameIdents(seuratDataSC, new.cluster.ids)
#   DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# otherwise, a you can pull from the "active.ident" field if it's been usefully populated:
new.cluster.ids <- paste0("CellType", levels(seuratDataSC@active.ident))
names(x = new.cluster.ids) <- levels(x = seuratDataSC)
seuratDataSC <- RenameIdents(object = seuratDataSC, new.cluster.ids)
DimPlot(seuratDataSC, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()

# SUBDIVIDING CLUSTERS WORKFLOW
#------------------------------

seuratDataSC[["ClusterNames_1.6"]] <- Idents(seuratDataSC) # stash current cluster names
seuratDataSC <- FindClusters(seuratDataSC, resolution = 2.0) # rerun clustering
plot1 <- DimPlot(seuratDataSC, reduction = "umap", label = TRUE) + NoLegend()
plot2 <- DimPlot(seuratDataSC, reduction = "umap", group.by = "ClusterNames_1.6", label = TRUE) + NoLegend()

# patchwork system
plot1 + plot2

# observe which group split, and can then find markers differentiating it. e.g. assuming group 0 split into new groups 0 and 1:
cell.markers <- FindMarkers(seuratDataSC, ident.1 = 0, ident.2 = 1)
FeaturePlot(seuratDataSC, features = c("S100A4", "CCR7"), cols = c("green", "blue"))

# DIFFERENTIAL EXPRESSION
#------------------------

# Alternative analyses
# Differential expression using t-test
FindMarkers(seuratDataSC, ident.1 = 0, ident.2 = 1, test.use = "t") # might have to use "CellType0" if you've done the new cluster id step

# Differential expression using ROC
FindMarkers(seuratDataSC, ident.1 = 3, test.use = "roc", only.pos = TRUE) # roc returns classification power from 0-1 for each individual marker

