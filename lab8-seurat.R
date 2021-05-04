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
seuratDataSC <- SCTransform(seuratData, vars.to.regress = "percent.mt", verbose = TRUE)

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
plotVariableFeatures
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

# Variable regression using sctransform
#   single command replaces NormalizeData, ScaleData, and FindVariableFeatures
#   INSERTED ABOVE, PRIOR TO MANUAL NORMALIZATION
# seuratDataSC <- SCTransform(seuratData, vars.to.regress = "percent.mt", verbose = FALSE)

all.genes <- rownames(seuratData) # equivalently: seuratData@assays$RNA@counts@Dimnames[1]
seuratData <- ScaleData(seuratData,  vars.to.regress = c("percent.mt")) # override default and scale all genes by adding `features = all.genes,`, not just highly-variable features


# DIMENSIONAL REDUCTION / PRINCIPAL COMPONENT ANALYSIS
#------------------------------------------------------

seuratData <- RunPCA(seuratData, features = VariableFeatures(seuratData), ndims.print = 1:5, nfeatures.print = 5) # equivalently, features = seuratData@assays$RNA@var.features but note that that means specifying assay RNA, whereas if you do SCTransform it's assay SCT. Avoid confusion with the VariableFeatures call instead.

VizDimLoadings(seuratData, dims = 1:2, reduction = "pca")
VizDimLoadings(seuratDataSC, dims = 1:2, reduction = "pca")
DimPlot(seuratData, reduction = "pca")
DimPlot(seuratDataSC, reduction = "pca")
