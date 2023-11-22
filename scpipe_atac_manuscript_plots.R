###############################
# Manuscript plots ------------
###############################

#' scPipe ATAC-Seq module
#' Authors: Shani Amarasinghe, Phil Yang, Haoyu Yang
#' Date: 07/14/2023

# Setup --------------------------------------

## Load required packages

#library(ggVennDiagram, lib.loc="/stornext/General/data/user_managed/grpu_mritchie_1/Shani/scPipe_final_modifications/R_libs")
library(cowplot)
library(Biostrings)
library(ggplot2)
library(ggpubr, lib.loc="/stornext/General/data/user_managed/grpu_mritchie_1/Shani/scPipe_final_modifications/R_libs")

library(Seurat)
library(Signac)
library(hdf5r,lib.loc="/stornext/General/data/user_managed/grpu_mritchie_1/Shani/scPipe_final_modifications/R_libs")
library(GenomeInfoDb, lib.loc="/stornext/General/data/user_managed/grpu_mritchie_1/Shani/scPipe_final_modifications/R_libs")
library(SingleCellExperiment)
library(dplyr)
library(stringr)
library(tibble)
library(tidyverse)
library(irlba)

# set up the working dir
setwd("/stornext/General/data/user_managed/grpu_mritchie_1/Shani/scPipe_final_modifications/")

# Custom functions ----------------------------------

#' @name signac_obj
#' @title create Seurat object to obtain features or counts per cell from cell ranger matrix
#' @param raw_peak_bc_matrix.h5 from outs dir of cellranger pipeline
#' @export
signac_obj <- function(file){
  counts <- Read10X_h5(filename =  list.files(file,full.names = T)[15])
  metadata <- read.csv(
    file =  list.files(file,full.names = T)[16],
    header = TRUE,
    row.names = 1
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments =  list.files(file,full.names = T)[7],
    min.cells = 0,
    min.features = 0
  )
  
  pbmc <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  return(pbmc[[]])
}


#' @name norm_rd
#' @title create a seurat object that is normalised and 
#' @param pbmc a seurat object
#' @export
norm_rd <- function(pbmc) {
  # Normalisation:  term frequency-inverse document frequency (TF-IDF)
  pbmc <- RunTFIDF(pbmc)
  # Feature Selection: All features
  pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
  # Dimension Reduction: singular value decomposition(SVD)
  pbmc <- RunSVD(pbmc)
  return(pbmc)
}


#' @name umap_pl
#' @title create a umap
#' @param pbmc seurat object
#' @param n_neighbours the number of neighbours to calculate the UMAP
#' @export
umap_pl <- function(pbmc, n_neighbours = 20){
  pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
  pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30, k.param = n_neighbours)
  pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
  return(pbmc)
}


#' @name stats_per_cell
#' @title create a df for counts and features er cell for scpipe output 
#' @details Because `filtered_stats_per_cell.csv` only have filtered data, a unfiltered version is needed.If I use unfiltered_feature_matrix.rds, the counts per cell is different from the filtered_stats_per_cell.csv file. If I use sparse_matrix.rds, the counts per cell is the same as the filtered_stats_per_cell.csv file.So some features are filtered out.unfiltered_feature_matrix.rds and features of sparse_matrix.rds 
#' @param path which has the cellranger out
#' @export
stats_per_cell = function(path){
  matrix=readRDS(paste0(path,"unfiltered_feature_matrix.rds"))
  ft=row.names(readRDS(paste0(path,"sparse_matrix.rds")))
  matrix = matrix[ft,]
  df= data.frame("counts_per_cell"=colSums(matrix), "features_per_cell" = colSums(matrix>0))
  row.names(df) = as.character(reverseComplement(DNAStringSet(row.names(df))))
  return(df)
}

#' @name cr_to_sce
#' @title creates a sce object from cellranger output and attaches cell line info
#' @param cr_dir cellranger dir
#' @export
# Note: don't need to take revcomp for cellranger output
cr_to_sce <- function(cr_dir) {
  cr_mat <- as(Matrix::readMM(file.path(cr_dir, "filtered_peak_bc_matrix", "matrix.mtx")), "CsparseMatrix")
  cr_bcs<- read.csv(file.path(cr_dir, "filtered_peak_bc_matrix", "barcodes.tsv"), header = FALSE)$V1
  cr_cell_data <- read.csv(file.path(cr_dir, "singlecell.csv"), header = TRUE)
  colnames(cr_mat) <- cr_bcs
  sce <- SingleCellExperiment(assays = list(counts = cr_mat))
  cell_data <- read.csv(file.path(cr_dir, "singlecell.csv"), header = TRUE) %>% column_to_rownames("barcode")
  colData(sce)[, colnames(cell_data)] <- cell_data[colnames(sce),]
  sce <- fix_sce_cols(sce, TRUE, FALSE)
  rownames(sce) <- 1:nrow(sce)
  return (sce)
}

#' @name fix_sce_cols
#' @title rev comp and edit the sce object for scpipe to match cellranger
#' @param cr_dir cellranger dir
#' @export
fix_sce_cols <- function(sce,
                         strip_last_two = FALSE,
                         revcomp = FALSE) {
  if (isTRUE(strip_last_two)) colnames(sce) <- gsub('.{2}$', '', colnames(sce))
  if (isTRUE(revcomp)) colnames(sce) <- Biostrings::reverseComplement(Biostrings::DNAStringSet(colnames(sce))) %>% as.character()
  
  return(sce)
}

#' @name attach_cell_line_info
#' @title Attaches cell line info to the colData of an sce object
#' @param sce The SingleCellExperiment object
#' @param cell_line_file Contains two columns, first being the barcodes, the second being the cell label
#' @export
attach_cell_line_info <- function(sce, cell_line_file) {
  message("Attaching cell line info")
  cli <- read.csv(cell_line_file, header=FALSE)
  colnames(cli) <- c("barcode", "cell_line")
  merged_cell_line_info <- base::merge(data.frame(barcode = colnames(sce)), cli, all.x = TRUE) %>% column_to_rownames("barcode")
  reordered_cli <- merged_cell_line_info[colnames(sce), , drop=FALSE]
  colData(sce) <- cbind(colData(sce), cell_line = reordered_cli$cell_line)
  return(sce)
}

#' @name attach_cell_comparison
#' @title compare the cell barcodes between two sce objects
#' @param sce first sce object
#' @param sce_to_compare the second sce object
#' @param sce_to_compare_name name of the second sce object
#' @export
# For a given sce object, add a new logical column in its colData indicating if the column name is present in `sce_to_compare`
attach_cell_comparison <- function(sce, sce_to_compare, sce_to_compare_name = "scPipe") {
  colData(sce)[, paste0("Overlap with ", sce_to_compare_name)] = colnames(sce) %in% colnames(sce_to_compare)
  return(sce)
}

#' @name run_seurat
#' @title generate a seurat object form the sce objects
#' @param sce first sce object
#' @param n_neighbours number of neighbours for clustering
#' @param clustering_comp logical to calculate values for ARI and NMI
#' @param dimplot logical to generate the UMAP
#' @export
run_seurat <- function(sce, n_neighbours = 300, clustering_comp = TRUE, dimplot = TRUE) {
  sce_seurat <- as.Seurat(sce, counts = "counts", data = "counts")
  sce_seurat_norm <- norm_rd(sce_seurat)
  sce_seurat_umap <- umap_pl(sce_seurat_norm, n_neighbours = n_neighbours)
  
  if (isTRUE(clustering_comp)) {
    ari <- mclust::adjustedRandIndex(sce_seurat_umap$cell_line, sce_seurat_umap$seurat_clusters)
    nmi <- NMI::NMI(data.frame(colnames(sce_seurat_umap), sce_seurat_umap$cell_line), 
                    data.frame(sce_seurat_umap$seurat_clusters) %>% rownames_to_column())$value
    cat("ARI:", ari, "\n")
    cat("NMI:", nmi, "\n")
  }
  
  if (isTRUE(dimplot))
    DimPlot(object = sce_seurat_umap, label = TRUE, group.by = c("cell_line", "seurat_clusters"))
}

#' @name sc_plot_umap
#' @title Generates UMAP of data from sce object
#' @description Uses feature count data from an sce object to produce a UMAP
#' @param mae The SingleCellExperiment object
#' @param by What to colour the points by. Needs to be in colData of all experiments. If is left `NULL` then will be coloured by cluster.
#' @param clustering_comp logical to calculate values for ARI and NMI
#' @param n_neighbours No. of neighbours for KNN 
#' @param output_file The path of the output file
#' @export
sc_plot_umap <- function(sce,
                         n_neighbours = 300,
                         s_neighbours = 350,
                         by = NULL,
                         clustering_comp = TRUE,
                         custom_colours = NULL,
                         output_file = NULL) {
  set.seed(123)
  
  sce_seurat <- as.Seurat(sce, counts = "counts", data = "counts")
  sce_seurat_norm <- norm_rd(sce_seurat)
  sce_seurat_umap <- umap_pl(sce_seurat_norm, n_neighbours = s_neighbours)
  
  if (isTRUE(clustering_comp)) {
    ari <- mclust::adjustedRandIndex(sce_seurat_umap$cell_line, sce_seurat_umap$seurat_clusters)
    nmi <- NMI::NMI(data.frame(colnames(sce_seurat_umap), sce_seurat_umap$cell_line), 
                    data.frame(sce_seurat_umap$seurat_clusters) %>% rownames_to_column())$value
    cat("ARI:", ari, "\n")
    cat("NMI:", nmi, "\n")
  }
  
  colData(sce) <- cbind(colData(sce), seurat_clusters = sce_seurat_umap$seurat_clusters)
  
  print(head(colData(sce)))
  
  counts <- assay(sce)
  bin_mat <- as.matrix((counts>0)+0)
  binary.mat <- TF.IDF.custom(bin_mat)
  dimnames(binary.mat) <- dimnames(bin_mat)
  n_bcs <- max(min(50, ncol(binary.mat), nrow(binary.mat))-1,0)
  mat.lsi          <- irlba(binary.mat, n_bcs)
  d_diagtsne       <- matrix(0, n_bcs, n_bcs)
  diag(d_diagtsne) <- mat.lsi$d
  mat_pcs          <- t(d_diagtsne %*% t(mat.lsi$v))
  rownames(mat_pcs)<- colnames(binary.mat)
  
  # clustering in the PCA space using KNN --------------
  knn.info<- RANN::nn2(mat_pcs, k = n_neighbours)
  
  ## convert to adjacency matrix
  knn           <- knn.info$nn.idx
  adj           <- matrix(0, nrow(mat_pcs), nrow(mat_pcs))
  rownames(adj) <- colnames(adj) <- rownames(mat_pcs)
  for(i in seq_len(nrow(mat_pcs))) {
    adj[i,rownames(mat_pcs)[knn[i,]]] <- 1
  }
  
  ## convert to graph
  g <- igraph::graph.adjacency(adj, mode="undirected")
  g <- igraph::simplify(g) ## remove self loops
  
  # identify communities, many algorithms. Use the Louvain clustering ------------
  km         <- igraph::cluster_louvain(g)
  com        <- km$membership
  names(com) <- km$names
  
  # running UMAP ------------------------------
  norm.data.umap    <- umap::umap(mat_pcs)
  
  df_umap           <- as.data.frame(norm.data.umap$layout)
  colnames(df_umap) <- c("UMAP1", "UMAP2")
  df_umap$barcode   <- rownames(mat_pcs)
  
  df_umap           <- dplyr::left_join(df_umap, enframe(com), by = c("barcode" = "name")) %>%
    dplyr::rename(cluster = value) %>%
    dplyr::mutate(cluster = as.factor(cluster))
  
  print(head(df_umap))
  
  if (!is.null(by) && !is.null(colData(sce)[[by]])) {
    sce_coldata <- colData(sce)[, c(by), drop=FALSE]
    df_umap[, colnames(sce_coldata)] = data.frame(sce_coldata)
  } else if (!is.null(by)) {
    df_umap[[by]] <- NA
    cat("Couldn't locate", by, "in the column data of the provided SCE object\n")
  }
  if (!is.null(by) && !is.null(df_umap[[by]])) {
    g <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(col = .data[[by]]), size = 2, alpha = 0.5) +
      theme_bw(base_size = 14)
    if (is.numeric(df_umap[[by]])) {
      g <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(col = .data[[by]]), size = 2, alpha = 0.5) +
        scale_colour_gradientn(colours=c("green","black")) +
        theme_bw(base_size = 14)
    }
  } else {
    # if (by =="cell_line") {
    #   cluster_col <- as.factor(sce@seurat_clusters[match(df_umap$barcode, colnames(sce)]) 
    # }
    g <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(col = cluster), size = 2, alpha = 0.5) +
      theme_bw(base_size = 14)
  }
  if (!is.null(custom_colours)) {
    g <- g + scale_color_manual(values = custom_colours)
  }
  if (!is.null(output_file)) ggsave(output_file)
  return(g)
}

#' @name TF.IDF.custom
#' @title Custom UMAP creation function
#' @description Uses feature count data from an sce object to produce a UMAP
#' @param binary.mat the binary matrix from the SC object
#' @export
TF.IDF.custom <- function(binary.mat, verbose = TRUE) {
  
  if (verbose) {
    message("Performing TF-IDF normalization")
  }
  object <- binary.mat
  npeaks       <- Matrix::colSums(x = object)
  tf           <- Matrix::tcrossprod(x = as.matrix(object), y = Matrix::Diagonal(x = 1 / npeaks))
  rsums        <- Matrix::rowSums(x = object)
  idf          <- ncol(x = object) / rsums
  norm.data    <- Matrix::Diagonal(n = length(x = idf), x = idf) %*% tf
  scale.factor <- 1e4
  slot(object = norm.data, name = "x") <- log1p(x = slot(object = norm.data, name = "x") * scale.factor)
  norm.data[which(x = is.na(x = norm.data))] <- 0
  return(norm.data)
}

#' @name cell_overlap_venn
#' @title compare the cell barcode overlap from called cells from cellranger to scPipe
#' @param sce_list sces to draw the venn from
#' @param output_file out filename
#' @export
# Draw a venn diagram comparing the cells of two sce objects
cell_overlap_venn <- function(sce_list, output_file = NULL) {
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  v <- VennDiagram::venn.diagram(x = lapply(sce_list, colnames),
                                 category.names = names(sce_list),
                                 filename = output_file,
                                 output = FALSE,
                                 fill = c("#a3fb98", "#ff8d8d"),
                                 cat.cex = 1.7,
                                 cat.default.pos = "outer",
                                 cex = 2.0)
  if (is.null(output_file)) {
    grid.newpage()
    grid.draw(v)
  }
}

#' @name plot_mito_comparison
#' @title compare the mito reads from cellranger to scPipe
#' @param sce_cr cellranger output into sce
#' @param sce_scPipe scPipe output into sce
#' @param output_file out filename
#' @export
# Require CellRanger sce object first then scPipe sce object
plot_mito_comparison <- function(sce_cr, sce_scPipe, output_file = NULL) {
  sce <- attach_cell_comparison(sce_cr, sce_scPipe)
  sce$frac_mito <- sce$mitochondrial/sce$total
  sce_coldata <- data.frame(colData(sce))
  colnames(sce_coldata) <- colnames(colData(sce))
  g <- ggplot(sce_coldata, aes(y = frac_mito, x=`Overlap with scPipe`, fill=`Overlap with scPipe`)) +
    geom_boxplot() + 
    ylab("Mitochondrial counts percentage") +
    theme(axis.text.x=element_blank(),
          text = element_text(size = 18)) +
    guides(fill=guide_legend(title="Overlap \nwith scPipe"))
  if (!is.null(output_file)) ggsave(output_file)
  return(g)
}

#' @name plot_mito_comparison_reverse
#' @title comppare the mito reads from scPipe to cellranger
#' @param sce_scPipe scPipe output into sce
#' @param sce_cr cellranger output into sce
#' @param output_file out filename
#' @export
# Comparing scPipe to CellRanger
plot_mito_comparison_reverse <- function(sce_scPipe, sce_cr, output_file = NULL) {
  sce <- attach_cell_comparison(sce_scPipe, sce_cr, "CellRanger")
  sce_coldata <- data.frame(colData(sce))
  colnames(sce_coldata) <- colnames(colData(sce))
  g <- ggplot(sce_coldata, aes(y = frac_mito, x=`Overlap with CellRanger`, fill=`Overlap with CellRanger`)) +
    geom_boxplot() + 
    ylab("Mitochondrial counts percentage") +
    theme(axis.text.x=element_blank(),
          text = element_text(size = 18)) +
    guides(fill=guide_legend(title="Overlap with \nCellRanger"))
  if (!is.null(output_file)) ggsave(output_file)
  return(g)
}

# Load the outputs from scPipe and Cell Ranger ATAC for Fig 2.B---------------

# scpipe
scpipe_80_dir <- "/stornext/General/data/user_managed/grpu_mritchie_1/PhilYang/scPipe_testing/haoyu_scMixology2_80_output/scPipe_atac_stats/"

# cell ranger
cr_80_dir <- "/stornext/Projects/promethion/promethion_access/lab_ritchie/scM_multiome_1/short_term/scMixology2_80_atac/outs/"

## Convert the stats file for each barcode to a dataframe

# scpipe output has filtered barcodes
# scpipe output need to convert to reverse complement
scpipe_80_sc <- read.csv(paste0(scpipe_80_dir,"filtered_stats_per_cell.csv"),header = T,stringsAsFactors = F)
scpipe_80_sc$cell = as.character(reverseComplement(DNAStringSet(scpipe_80_sc$cell))) 
rownames(scpipe_80_sc) = scpipe_80_sc$cell

# cell ranger output has all barcodes, so will filter later
cr_80_sc <- read.csv(paste0(cr_80_dir,"singlecell.csv"),header = T,stringsAsFactors = F)

# Select cell barcodes from cell ranger output, which has `is__cell_barcode` being 1.
# cell ranger atac output barcodes are atac barcodes, not rna
cr_80_sc_cell <- cr_80_sc[cr_80_sc$is__cell_barcode==1,"barcode"] 
# remove -1 at the end
cr_80_sc_cell <- gsub("-1","",cr_80_sc_cell) 

## Retrieve a DF for the stats

# create signac object and retrieve a stat DF for cellranger
cr_80_stats <- signac_obj(cr_80_dir)

# create a stat DF for scPipe
scpipe_80_stats <- stats_per_cell("/stornext/General/data/user_managed/grpu_mritchie_1/PhilYang/scPipe_testing/haoyu_scMixology2_80_output/")

# Load the outputs from scPipe and Cell Ranger ATAC for Fig 2.D - Fig 2.J ---------------

## Load data 

# Cell line info
atac_80_cli <- "/stornext/General/data/user_managed/grpu_mritchie_1/PhilYang/fastq_data/scMixology2_80_fastqs/ATAC_cell_line_info.csv"

# data
# CellRanger
cr_atac_80 <- "/stornext/Projects/promethion/promethion_data/short_term/lab_ritchie/scM_multiome_1/scMixology2_80_atac/outs/"
# scPipe
scPipe_atac_80 <- "/stornext/General/data/user_managed/grpu_mritchie_1/PhilYang/scPipe_testing/haoyu_scMixology2_80_output"

# extra stats foile for cellranger 
stats_80 <- t(read.csv(file.path(cr_atac_80, "summary.csv"), header = TRUE))

## Generate sce objects
sce_80_cr     <- attach_cell_line_info(cr_to_sce(cr_atac_80), atac_80_cli)
sce_80_scPipe <- attach_cell_line_info(fix_sce_cols(readRDS(file.path(scPipe_atac_80, "scPipe_atac_SCEobject.rds")), revcomp = TRUE), atac_80_cli)


# 1. Figure 2B. Correlation plot | reads/cell -------------

union_80      <- union(cr_80_sc_cell,scpipe_80_sc$cell)
counts_80_all <- data.frame("CellRanger"=cr_80_stats[union_80,"nCount_peaks"],"scPipe"=scpipe_80_stats[union_80,"counts_per_cell"])

# Add source column: in CellRanger will be 2, and in scPipe will be 1, in Both will be 3. Change them to the source name.
counts_80_all$source <- factor((union_80 %in% cr_80_sc_cell*2) + union_80 %in% scpipe_80_sc$cell)
levels(counts_80_all$source)[levels(counts_80_all$source)=='1'] <- 'scPipe'
levels(counts_80_all$source)[levels(counts_80_all$source)=='2'] <- 'CellRanger'
levels(counts_80_all$source)[levels(counts_80_all$source)=='3'] <- 'Common'

counts_80_all$CellRanger[which(is.na(counts_80_all$CellRanger))]         <- 0
counts_80_all$scPipe[which(is.na(counts_80_all$scPipe))]                 <- 0

p5 <- ggscatter(data    = counts_80_all, 
                x       = "scPipe", 
                y       = "CellRanger",
                color   ="source",
                palette ="Dark2",
                size    = 4,
                font.label = c(20, "plain"),
                add     = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
                cor.coef.size = 12)+
  geom_rug(col="steelblue",alpha=0.1, size=1.5) + 
  theme(text = element_text(size = 18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

pdf("./plots/scartter_80.pdf")
  p5
dev.off()

pdf("./plots/scatter_hist_80.pdf")
ggscatterhist(data    = counts_80_all, 
                    x       = "scPipe", 
                    y       = "CellRanger",
                    color   ="source",
                    palette ="Dark2",
                    size    = 4,
                    font.label = c(20, "plain"),
                    add     = "reg.line",  # Add regression line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
                    cor.coef.size = 12,
                    margin.plot = c("density"),
                    margin.params = list(fill = "source", color = "source", alpha=.25),
                    margin.ggtheme = theme_pubr(),
                    margin.space = FALSE,
                    main.plot.size = 2,
                    margin.plot.size = 1,
                    legend = "bottom") +
  geom_rug(col="steelblue",alpha=0.1, size=1.5) + 
  theme(text = element_text(size = 18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

# pdf("./plots/scatter_hist_80.pdf")
# p6
dev.off()

### 2. Figure 2C. Correlation plot | features/cell -------------

features_80_all = data.frame("CellRanger"=cr_80_stats[union_80,"nFeature_peaks"],"scPipe"=scpipe_80_stats[union_80,"features_per_cell"])

# Add source column: in CellRanger will be 2, and in scPipe will be 1, in Both will be 3. Change them to the source name.
features_80_all$source = factor((union_80 %in% cr_80_sc_cell*2) + union_80 %in% scpipe_80_sc$cell)
levels(features_80_all$source)[levels(features_80_all$source)=='1'] <- 'scPipe'
levels(features_80_all$source)[levels(features_80_all$source)=='2'] <- 'CellRanger'
levels(features_80_all$source)[levels(features_80_all$source)=='3'] <- 'Common'

# convert NA to 0
features_80_all$CellRanger[which(is.na(features_80_all$CellRanger))] <- 0
features_80_all$scPipe[which(is.na(features_80_all$scPipe))]         <- 0

p9 <- ggscatter(data       = features_80_all, 
                x          = "scPipe", 
                y          = "CellRanger",
                color      = "source",
                palette    = "Dark2",
                size       = 4,
                font.label = c(20, "plain"),
                add        = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int   = TRUE, # Add confidence interval
                cor.coef   = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
                cor.coef.size = 12)+
  geom_rug(col="steelblue",alpha=0.1, size=1.5) + 
  theme(text = element_text(size = 18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))


pdf("./plots/scatter_features_80.pdf")
  p9
dev.off()

pdf("./plots/scatter_features_hist_80.pdf")
ggscatterhist(data       = features_80_all, 
                     x          = "scPipe", 
                     y          = "CellRanger",
                     color      = "source",
                     palette    = "Dark2",
                     size       = 4,
                     font.label = c(20, "plain"),
                     add        = "reg.line",  # Add regressin line
                     add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                     conf.int   = TRUE, # Add confidence interval
                     cor.coef   = TRUE, # Add correlation coefficient. see ?stat_cor
                     cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
                     cor.coef.size = 12,
                     margin.plot = c("density"),
                     margin.params = list(fill = "source", color = "source", alpha=.25),
                     margin.ggtheme = theme_pubr(),
                     margin.space = FALSE,
                     main.plot.size = 2,
                     margin.plot.size = 1,
                     legend = "bottom") +
  geom_rug(col="steelblue",alpha=0.1, size=1.5) + 
  theme(text = element_text(size = 20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

# pdf("./plots/scatter_features_hist_80.pdf")
# p10
dev.off()

# 2. Figure 2D. Venn diagram -------------

# first a bit of more setting up to do

output_folder <- "./plots"

cell_overlap_venn(list(CellRanger = sce_80_cr,
                       scPipe = sce_80_scPipe), output_file = file.path(output_folder, "venn_80_2.png"))


# x_80 <- list("CellRanger"=cr_80_sc_cell,
#              "scPipe"=scpipe_80_sc$cell) 
# 
# 
# ggVennDiagram(x_80)

# Most of the barcodes (\>80%) overlap between cell ranger and scPipe. Cell ranger keep more barcodes.

# 3. Figure 2D & 2E: Mitochondria percentages ----------

output_folder <- "./plots"

plot_mito_comparison(sce_80_cr, sce_80_scPipe, output_file = file.path(output_folder, "compare_mito_80.png"))
plot_mito_comparison_reverse(sce_80_scPipe, sce_80_cr, output_file = file.path(output_folder, "compare_mito_80_rev.png"))


# 4. Figure 2E: Seurat map and ARI,NMI ------------------

output_folder <- "./plots/"

# run seurat for the both data sets

umap_cr     <- run_seurat(sce_80_cr, n_neighbours = 300)
# ARI: 0.5423382 
# NMI: 0.6203453
umap_scPipe <- run_seurat(sce_80_scPipe, n_neighbours = 300)
# ARI: 0.6847724 
# NMI: 0.7536758 

pdf(paste0(output_folder, "umap_cr.pdf"))
  print(umap_cr)
dev.off()

pdf(paste0(output_folder, "umap_scPipe.pdf"))
  print(umap_scPipe)
dev.off()

sce_80_comp <- attach_cell_comparison(sce_80_cr, sce_80_scPipe)
umap_comp   <- run_seurat(sce_80_comp[, sce_80_comp$`Overlap with scPipe` == TRUE], n_neighbours = 300)

pdf(paste0(output_folder, "umap_comp.pdf"))
  print(umap_comp)
dev.off()

# # Comparison of cells on UMAP
# sce_comp_cr <- attach_cell_comparison(sce_80_cr, sce_80_scPipe, "scPipe")
# sc_plot_umap(sce_comp_cr, by = "Overlap with scPipe", output_file = file.path(output_folder, "compare_umap_80.png"), custom_colours = c("red", "dimgrey"), n_neighbours = 300)
# sc_plot_umap(sce_comp_cr, output_file = file.path(output_folder, "cluster_umap_80.png"), n_neighbours = 300)
# sc_plot_umap(sce_comp_cr, by = "cell_line", output_file = file.path(output_folder, "cluster_umap_80_cellline.png"), n_neighbours = 300)
# 
# sce_comp_scPipe <- attach_cell_comparison(sce_80_scPipe, sce_80_cr, "CellRanger")
# sc_plot_umap(sce_comp_scPipe, by = "Overlap with CellRanger", output_file = file.path(output_folder, "compare_umap_80_rev.png"), custom_colours = c("red", "dimgrey"), n_neighbours = 300)
# sc_plot_umap(sce_comp_scPipe, output_file = file.path(output_folder, "cluster_umap_80_rev.png"), n_neighbours = 300)
# sc_plot_umap(sce_comp_scPipe, by = "cell_line", output_file = file.path(output_folder, "cluster_umap_80__rev_cellline.png"), n_neighbours = 300)

sce_comp_cr <- attach_cell_comparison(sce_80_cr, sce_80_scPipe, "scPipe")
plot1 <- sc_plot_umap(sce_comp_cr, by = "Overlap with scPipe", output_file = file.path(output_folder, "compare_umap_80.png"), custom_colours = c("red", "dimgrey"), n_neighbours = 350)
sc_plot_umap(sce_comp_cr, output_file = file.path(output_folder, "cluster_umap_80.png"), n_neighbours = 350, custom_colours = c("#EF8860", "#57D635", "#D6E8D9", "#FF7994", "#1F3D51"))
sc_plot_umap(sce_comp_cr, by = "cell_line", output_file = file.path(output_folder, "cluster_umap_80_cellline.png"), custom_colours = c("#EF8860", "#57D635", "#D6E8D9", "#FF7994", "#1F3D51"), n_neighbours = 350)
sc_plot_umap(sce_comp_scPipe, by = "counts_per_cell", output_file = file.path(output_folder, "cluster_umap_80_counts.png"), n_neighbours = 350)

sce_comp_scPipe <- attach_cell_comparison(sce_80_scPipe, sce_80_cr, "CellRanger")
plot2 <- sc_plot_umap(sce_comp_scPipe, by = "Overlap with CellRanger", output_file = file.path(output_folder, "compare_umap_80_rev.png"), custom_colours = c("red", "dimgrey"), n_neighbours = 300)
sc_plot_umap(sce_comp_scPipe, output_file = file.path(output_folder, "cluster_umap_80_rev.png"), n_neighbours = 300, custom_colours = c("#EF8860", "#57D635", "#D6E8D9", "#FF7994", "#1F3D51"))
sc_plot_umap(sce_comp_scPipe, by = "cell_line", output_file = file.path(output_folder, "cluster_umap_80_rev_cellline.png"), custom_colours = c("#EF8860", "#57D635", "#D6E8D9", "#FF7994", "#1F3D51"), n_neighbours = 300)
# sc_plot_umap(sce_comp_scPipe, by = "cell_line", output_file = file.path(output_folder, "cluster_umap_80_rev_cellline.png"), n_neighbours = 350)
sc_plot_umap(sce_comp_scPipe, by = "counts_per_cell", output_file = file.path(output_folder, "cluster_umap_80_rev_counts.png"), n_neighbours = 300)

saveRDS(plot1, "./cr_comp_ggplot.RDS")
saveRDS(plot2, "./scpipe_comp_ggplot.RDS")