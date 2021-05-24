library(dplyr)
library(Matrix)
library(gdata)
library(patchwork)
library(sctransform)
library(Seurat)

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
#counts <- counts[,1:1000] # use only first 1000 cells, for tutorial runspeed purposes only
countsPerCell <- Matrix::colSums(counts)
countsPerGene <- Matrix::rowSums(counts)
genesPerCell <- Matrix::colSums(counts > 0) # count gene only if it has non-zero reads mapped
cellsPerGene <- Matrix::rowSums(counts > 0) # count cells only where the gene is expressed

hist(log10(countsPerCell + 1), main = 'Counts per Cell', col = 'wheat')
hist(log10(genesPerCell + 1), main = 'Genes per Cell', col = 'wheat')
hist(log10(countsPerGene + 1), main = 'Counts per Gene', col = 'wheat')
hist(log10(cellsPerGene + 1), main = 'Cells per Gene', col = 'wheat')
plot(countsPerCell, genesPerCell, log = 'xy', col = 'dodgerblue4') + title('counts vs genes per cell')
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
#VlnPlot(seuratData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

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
                                          n.exp.hkgenes > 55)

seuratData
VlnPlot(seuratData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "n.exp.hkgenes"), ncol = 4)
plotMTpercentVsUMIs <- FeatureScatter(seuratData, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotGenesVsUMIs <- FeatureScatter(seuratData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotMTpercentVsUMIs + plotGenesVsUMIs

# Normalization
seuratData <- NormalizeData(seuratData, normalization.method = "LogNormalize", scale.factor = 10000) # defaults


# VARIABLE FEATURE SELECTION
#----------------------------
# sets identity information

seuratData <- FindVariableFeatures(seuratData, selection.method = "vst",
                                   nfeatures = 2000) # defaults
                                   # nfeatures = dim(seuratData@assays$RNA@counts)[1]*.003) # 3 SDs out

topVariance <- head(VariableFeatures(seuratData), 10) # 10 most highly-variable genes
plotVariableFeatures <- VariableFeaturePlot(seuratData)
plotLabeledFeatures <- LabelPoints(plotVariableFeatures, points = topVariance, repel = TRUE, xnudge = 0, ynudge = 0)
plotLabeledFeatures # is 'standardized variance' z-score, i.e. number of SDs from the mean? No. When nfeatures excludes the bottom 99.7%, i.e. 3 SDs, standardized variance threshold is about 12. It's just a particular calculation, as far as I can determine from a brief look at https://www.biorxiv.org/content/biorxiv/early/2018/11/02/460147.full.pdf


# GENE SET SCORING
#------------------

# Markers of G2/M phase and markers of S phase.
s.genes <- Seurat::cc.genes$s.genes
s.genes <- s.genes[s.genes %in% rownames(seuratData)] # genes in dataset
g2m.genes <- Seurat::cc.genes$g2m.genes
g2m.genes <- g2m.genes[g2m.genes %in% rownames(seuratData)] # genes in dataset
seuratData <- CellCycleScoring(seuratData, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
FeatureScatter(seuratData, "S.Score", "nFeature_RNA")

# Genes upregulated during dissociation of tissue into single cells, used to calculate a disassociation score
genes.dissoc <- list(c("ATF3", "BTG2", "CEBPB", "CEBPD", "CXCL3", "CXCL2", "CXCL1", "DNAJA1", "DNAJB1", "DUSP1", "EGR1", "FOS", "FOSB", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA1A", "HSPA1B", "HSPA8", "HSPB1", "HSPE1", "HSPH1", "ID3", "IER2", "JUN", "JUNB", "JUND", "MT1X", "NFKBIA", "NR4A1", "PPP1R15A", "SOCS3", "ZFP36"))
seuratData <- AddModuleScore(seuratData, features = genes.dissoc, ctrl.size = 20, name = "genes_dissoc")
FeatureScatter(seuratData, "genes_dissoc1", "nFeature_RNA")

# SCALING
#---------
# Important to scale all genes that will be used as input to PCA or heatmaps.
# The Seurat ScaleData function
#   - shifts expression so that mean expression across cells is 0 for each gene
#   - scales expression so that the variance across cells is 1 for each gene
#       this ensures genes are weighted evenly in downstream analyses; otherwise highly-expressed genes would dominate
#   - stores this data to seuratData[["RNA"]]@scale.data (equivalently: seuratData@assays$RNA@scale.data)
#   By default, only operates on previously identified highly-variable features
#   Optionally, use vars.to.regress = "percent.mt" etc, though better to use sctransform()'s vars.to.regress for this; see https://www.biorxiv.org/content/10.1101/576827v2 and https://satijalab.org/seurat/v3.0/sctransform_vignette.html

# Variable regression using sctransform - for advanced users
#   single command replaces NormalizeData, ScaleData, and FindVariableFeatures
#   if so, insert above, prior to manual normalization, but breaks JackStraw dimensional analysis later
#   SCTransform substantially mitigates the confounder of variation in sequencing depth compared to standard log-normalization (source: https://satijalab.org/seurat/archive/v3.0/sctransform_vignette.html 'Why can we choose more PCs...')
#   It also an entire end-to-end workflow more stable in a single, pithy function:
#     pbmc <- CreateSeuratObject(pbmc_data) %>% PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
#     SCTransform(vars.to.regress = "percent.mt") %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>%
#     RunUMAP(dims = 1:30) %>% FindClusters()
# seuratData <- SCTransform(seuratData, vars.to.regress = "percent.mt", verbose = FALSE)

all.genes <- rownames(seuratData) # equivalently: seuratData@assays$RNA@counts@Dimnames[1]
seuratData <- ScaleData(seuratData,  vars.to.regress = c("percent.mt")) # override default and scale all genes by adding `features = all.genes,`, not just highly-variable features


# DIMENSIONAL REDUCTION / PRINCIPAL COMPONENT ANALYSIS
#------------------------------------------------------

set.seed(2020) ## used for reporducibility for tutorial

seuratData <- RunPCA(seuratData, features = VariableFeatures(seuratData), ndims.print = 1:5, nfeatures.print = 5) # equivalently, features = seuratData@assays$RNA@var.features but note that that means specifying assay RNA, whereas if you do SCTransform it's assay SCT. Avoid confusion with the VariableFeatures call instead.
VizDimLoadings(seuratData, dims = 1:2, reduction = "pca")
DimPlot(seuratData, reduction = "pca")

seuratData <- RunICA(seuratData, features = VariableFeatures(seuratData), ndims.print = 1:5, nfeatures.print = 5)
VizDimLoadings(seuratData, dims = 1:2, reduction = "ica")
DimPlot(seuratData, reduction = "ica")

# ProjectDim scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
seuratData <- ProjectDim(seuratData, reduction = "pca")


# GENES BY PCs HEATMAP
#----------------------

DimHeatmap(seuratData, dims = 1:6, cells = 300, reduction = "pca", balanced = TRUE)
DimHeatmap(seuratData, dims = 1:6, cells = 300, reduction = "ica", balanced = TRUE)


# DIMENSIONALITY DETERMINATION
#------------------------------
# First step in clustering is to determine the 'true' dimensionality of the dataset.
# Determine dimensionality of the dataset using the JackStraw procedure from Macosko et al (http://www.cell.com/abstract/S0092-8674(15)00549-8)
# JackStraw cannot be used on SCT-transformed data - JackStraw assumes each gene has equal variance, but SCT normalization weights gene variance by biological heterogeneity
# Contrast with less-computationally-demanding ElbowPlot

seuratData <- JackStraw(seuratData, reduction = "pca")
seuratData <- ScoreJackStraw(seuratData, dims = 1:20)
JackStrawPlot(seuratData, dims = 1:20)
PCASigGenes(seuratData, pcs.use = 1, pval.cut = 0.001)[1:20] # No idea what this is supposed to illustrate
ElbowPlot(seuratData, ndims = 50, reduction = "pca")


# CLUSTERING
#------------

set.seed(2020)
# Seurat v3 has a FindNeighbors() that constructs a SNN/KNN graph (based on euclidian dist in PCA space), with edge weights refined by Jaccard similarity (shared overlap in local neighborhood)
seuratData <- FindNeighbors(seuratData, dims = 1:10) # where dimensionality used (10) is determined in the step above, in this case matching tutorial which uses a fairly severe cutoff

# Seurat v3 also has a FindClusters() that does modularity optimization techniques (Louvain by default, or SLM)
seuratData <- FindClusters(seuratData, resolution = 0.6) # where resolution affects granularity of downstream clusters, increasing->more clusters.  Typically use 0.4-1.2 for sc datasets of around 3K cells, increasing as necessary for larger datasets. My use of 0.6 is to match tutorial.
head(Idents(seuratData),5) # look at cluster IDs of the first 5 cells

# Seurat v3 offers UMAP and tSNE to visualize clustering
seuratData <- RunTSNE(seuratData, reduction.use = "pca", dims.use = 1:10, perplexity = 10)
DimPlot(seuratData, reduction = "tsne", label = TRUE) + NoLegend()
set.seed(2020)
seuratData <- RunUMAP(seuratData, dims = 1:10)
DimPlot(seuratData, reduction = "umap", label = TRUE) + NoLegend()



# DATA EXPORT
#-------------

saveRDS(seuratData, file = paste0(datadir,"output/seuratData_tutorial.rds")) # Can be loaded back in without rerunning computationally intense steps
seuratData <- readRDS(file = paste0(datadir,"output/seuratData_tutorial.rds")) # load


# DIFFERENTIALLY EXPRESSED FEATURES / CLUSTER BIOMARKERS
#--------------------------------------------------------

cluster1.markers <- FindMarkers(seuratData, ident.1 = 1, min.pct = 0.1) # cluster to ID is ident.1 (to do all clusters, use FindAllMarkers instead). Min.pct is detection threshold, thresh.test is to require differential expression amount between two groups. Finally, max.cells.per.ident can speed computation for a small loss in power
head(cluster1.markers)

cluster5.markers <- FindMarkers(seuratData, ident.1 = 5, ident2 = c(0,1)) # find all markers distinguishing cluster 5 from clusters 0 and 1
head(cluster5.markers)

cluster411.markers <- FindMarkers(seuratData, ident.1 = c(4,11), ident2 = c(2,3,6), thresh.use = 25, only.pos = TRUE)
head(cluster411.markers)

nsclc.markers <- FindAllMarkers(seuratData, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
saveRDS(nsclc.markers, file = paste0(datadir,"output/nsclc.markers.rds"))
nsclc.markers <- readRDS(file = paste0(datadir,"output/nsclc.markers.rds"))

nsclc.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
# nsclc.markers %>% filter(cluster == 1 | cluster == 2 | cluster == 3) %>% filter(n() > 2 & avg_log2FC > 1.4) # example of finding genes common to multiple clusters
top10bycluster <- nsclc.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

VlnPlot(seuratData, features = c("MS4A1", "CD79A"))
VlnPlot(seuratData, features = c("NKG7", "PF4"), log = TRUE)
FeaturePlot(seuratData, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols = c("grey", "blue"), reduction = "tsne")
DoHeatmap(seuratData, features = top10bycluster$gene, label = TRUE, size = 3)

# RENAMING CLUSTERS
#-------------------

# if you know from looking cluster membership and recognize cell types from known markers, you can manually specify:
#   new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
#   names(new.cluster.ids) <- levels(seuratData)
#   seuratData <- RenameIdents(seuratData, new.cluster.ids)
#   DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# otherwise, a you can pull from the "active.ident" field if it's been usefully populated:
new.cluster.ids <- paste0("CellType", levels(seuratData@active.ident))
names(x = new.cluster.ids) <- levels(x = seuratData)
seuratData <- RenameIdents(object = seuratData, new.cluster.ids)
DimPlot(seuratData, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()

# SUBDIVIDING CLUSTERS WORKFLOW
#------------------------------

seuratData[["ClusterNames_0.6"]] <- Idents(seuratData) # stash current cluster names
seuratData <- FindClusters(seuratData, resolution = 0.8) # rerun clustering
plot1 <- DimPlot(seuratData, reduction = "tsne", label = TRUE) + NoLegend()
plot2 <- DimPlot(seuratData, reduction = "tsne", group.by = "ClusterNames_0.6", label = TRUE) + NoLegend()

# patchwork system
plot1 + plot2

# observe which group split, and can then find markers differentiating it. e.g. assuming group 0 split into new groups 0 and 1:
cell.markers <- FindMarkers(seuratData, ident.1 = 0, ident.2 = 1)
FeaturePlot(seuratData, features = c("S100A4", "CCR7"), cols = c("green", "blue"))

# DIFFERENTIAL EXPRESSION
#------------------------

# Alternative analyses
# Differential expression using t-test
FindMarkers(seuratData, ident.1 = 0, ident.2 = 1, test.use = "t") # might have to use "CellType0" if you've done the new cluster id step

# Differential expression using ROC
FindMarkers(seuratData, ident.1 = 3, test.use = "roc", only.pos = TRUE) # roc returns classification power from 0-1 for each individual marker


# CLASSIFIER INTEGRATION TO CHECK CLUSTERING
#--------------------------------------------
# see also (might actually be unrelated?) https://satijalab.org/seurat/archive/v3.0/integration.html#standard-workflow
# must use log-normalized RNA data, not SCTransformed

# Assign the test object a three level attribute
groups <- sample(c("train", "test"), size = NROW(seuratData@meta.data), replace = TRUE, prob = c(0.8, 0.2))
names(groups) <- colnames(seuratData)
seuratData <- AddMetaData(seuratData, metadata = groups, col.name = "group")

# Find Anchors
seuratData.list <- SplitObject(seuratData, split.by = "group")
seuratData.anchors <- FindIntegrationAnchors(seuratData.list, dims = 1:30)
seuratData.integrated <- IntegrateData(seuratData.anchors, dims = 1:30)
seuratData.query <- seuratData.list[["train"]]
seuratData.anchors <- FindTransferAnchors(seuratData.integrated, query = seuratData.query, dims = 1:30)
predictions <- TransferData(seuratData.anchors, refdata = seuratData.integrated$ClusterNames_0.6, dims = 1:30)
seuratData <- AddMetaData(seuratData.query, metadata = predictions)
table(seuratData@meta.data$ClusterNames_0.6, seuratData@meta.data$predicted.id)
# I guess the table is supposed to show predictions using a classifier, vs existing clusters? So eg 00 is predicted and existing both 0, 02 is maybe existing 0 predicted 2, I think
# so a table with a lot of 0s is good because it means stuff isn't widely scattered among different clusters

# PROBABILISTIC CLUSTERING (LDA/TOPICS)
#---------------------------------------
# see https://broadinstitute.github.io/2020_scWorkshop/feature-selection-and-cluster-analysis.html#probabilistic-lda-clustering
# omitting here because tutorial uses a different dataset
