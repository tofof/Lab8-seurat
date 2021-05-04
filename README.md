# scRNA-seq Data Processing Exploration

Following tutorials at
   * https://broadinstitute.github.io/2020_scWorkshop/data-wrangling-scrnaseq.html#examine-contents-of-seurat-object
   * https://broadinstitute.github.io/KrumlovSingleCellWorkshop2020/data-wrangling-scrnaseq-1.html
   * https://satijalab.org/seurat/archive/v3.0/pbmc3k_tutorial.html

Public dataset of Non-Small Cell Lung Cancer Cells (NSCLC) freely available from 10X Genomics
>   curl -O https://cf.10xgenomics.com/samples/cell-vdj/2.2.0/vdj_v1_hs_nsclc_5gex/vdj_v1_hs_nsclc_5gex_filtered_gene_bc_matrices.tar.gz
>    tar -xvzf vdj_v1_hs_nsclc_5gex_filtered_gene_bc_matrices.tar.gz


Tirosh et al housekeeping gene list from
   * https://github.com/Michorlab/tnbc_scrnaseq/blob/master/data/housekeepers.txt
