library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(patchwork)

# PREAMBLE
#----------
# Following tutorials at
#   https://broadinstitute.github.io/2020_scWorkshop/data-wrangling-scrnaseq.html#examine-contents-of-seurat-object
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
genesPerCell <- Matrix::colSums(counts>0) # count gene only if it has non-zero reads mapped
cellsPerGene <- Matrix::rowSums(counts>0) # count cells only where the gene is expressed

hist(log10(countsPerCell+1),main='Counts per Cell', col='wheat')
hist(log10(genesPerCell+1),main='Genes per Cell', col='wheat')
hist(log10(countsPerGene+1),main='Counts per Gene', col='wheat')
hist(log10(cellsPerGene+1),main='Cells per Gene', col='wheat')
plot(countsPerCell, genesPerCell, log='xy', col='dodgerblue4')
plot(cellsPerGene,countsPerGene, log='xy', col='dodgerblue4')
plot(sort(genesPerCell), log='y', xlab='Cell', ylab='Complexity', main='Genes per Cell (ordered)')

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
plot(sort(genesPerCell_sd), log='y', xlab='Cell', ylab='Complexity', main='Genes per Cell (ordered, cells > 3, genes > 350)')

# QUALITY CONTROL
#-----------------

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
seuratData <- subset(seuratData, subset = nFeature_RNA > 350 & nFeature_RNA < 6000)
seuratData <- subset(seuratData, subset = nCount_RNA < 60000)
seuratData <- subset(seuratData, subset = percent.mt < 12.5 )
seuratData <- subset(seuratData, subset = n.exp.hkgenes >= 55)

seuratData
VlnPlot(seuratData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "n.exp.hkgenes"), ncol = 4)
plotMTpercentVsUMIs <- FeatureScatter(seuratData, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotGenesVsUMIs <- FeatureScatter(seuratData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotMTpercentVsUMIs + plotGenesVsUMIs
