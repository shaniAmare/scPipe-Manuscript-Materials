#' scPipe ATAC-Seq manuacript figures - custom functions
#' Authors: Phil Yang, Haoyu Yang, Shani Amarasinghe
#' Date: 07/14/2023

## Helper scripts --------------------------------------


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
umap_pl <- function(pbmc, n_neighbours = 300){
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
run_seurat <- function(sce, n_neighbours = 300, clustering_comp = TRUE, dimplot = TRUE, output_file = NULL, custom_colours = NULL, by = NULL) {
  set.seed(123)
  sce_seurat      <- as.Seurat(sce, counts = "counts", data = "counts")
  sce_seurat_norm <- norm_rd(sce_seurat)
  sce_seurat_umap <- umap_pl(sce_seurat_norm, n_neighbours = n_neighbours)
  
  if (isTRUE(clustering_comp)) {
    ari <- mclust::adjustedRandIndex(sce_seurat_umap$cell_line, sce_seurat_umap$seurat_clusters)
    nmi <- NMI::NMI(data.frame(colnames(sce_seurat_umap), sce_seurat_umap$cell_line), 
                    data.frame(sce_seurat_umap$seurat_clusters) %>% rownames_to_column())$value
    cat("ARI:", ari, "\n")
    cat("NMI:", nmi, "\n")
  }
  
  
  if (isTRUE(dimplot)) {
    if (!is.null(by) && !is.null(sce_seurat_umap@meta.data[[by]])) {
      if (is.numeric(sce_seurat_umap@meta.data[[by]])) {
          sce_seurat_umap$log_total <- log(sce_seurat_umap@meta.data[[by]])
          dplot <- FeaturePlot(object = sce_seurat_umap, label = FALSE, features = "log_total") +
          scale_colour_gradientn(colours=c("green","black")) +
          theme_minimal(base_size = 16)
      } else{
        dplot <- DimPlot(object = sce_seurat_umap, label = FALSE, group.by = by) +
          theme_minimal(base_size = 16)
      }
    } else{
      dplot <- DimPlot(object = sce_seurat_umap, label = FALSE, group.by = c("cell_line", "seurat_clusters")) +
        theme_minimal(base_size = 16)
      }
    if (!is.null(custom_colours)) {
      dplot <- dplot + scale_color_manual(values = custom_colours) +
        theme_minimal(base_size = 16)
      }
    if (!is.null(output_file)) ggsave(output_file)
  }
}

#' @name sc_generate_umap_df
#' @title Generates a data frame for UMAP from sce object
#' @description Uses feature count data from an sce object to produce a UMAP
#' @param sce The SingleCellExperiment object
#' @param by What to colour the points by. Needs to be in colData of all experiments. If is left `NULL` then will be coloured by cluster.
#' @param clustering_comp logical to calculate values for ARI and NMI
#' @param n_neighbours No. of neighbours for KNN 
#' @param output_file The path of the output file
#' @export
sc_generate_umap_df <- function(sce,
                         n_neighbours = 300,
                         s_neighbours = 350,
                         by           = NULL,
                         clustering_comp = TRUE,
                         custom_colours = NULL,
                         output_file = NULL) {
  set.seed(123)
  
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
  
  colData(sce) <- cbind(colData(sce), seurat_clusters = sce_seurat_umap$seurat_clusters)
  
  print(head(colData(sce)))
  
  #saveRDS(sce, "./sce_comp_scPipe_halfway.rds")
  
  counts <- assay(sce)
  bin_mat <- as.matrix((counts>0)+0)
  binary.mat <- TF.IDF.custom(bin_mat)
  dimnames(binary.mat) <- dimnames(bin_mat)
  
  #saveRDS(binary.mat, "./sce_comp_binary_mat.rds")
  
  n_bcs <- max(min(50, ncol(binary.mat), nrow(binary.mat))-1,0)
  mat.lsi          <- irlba(binary.mat, n_bcs)
  d_diagtsne       <- matrix(0, n_bcs, n_bcs)
  diag(d_diagtsne) <- mat.lsi$d
  mat_pcs          <- t(d_diagtsne %*% t(mat.lsi$v))
  rownames(mat_pcs)<- colnames(binary.mat)
  
  #saveRDS(mat_pcs, "./sce_comp_mat_pcs.rds")
  
  # clustering in the PCA space using KNN --------------
  knn.info<- RANN::nn2(mat_pcs, k = n_neighbours)
  
  #saveRDS(knn.info, "./sce_comp_knn_info.rds")
  
  ## convert to adjacency matrix
  knn           <- knn.info$nn.idx
  adj           <- matrix(0, nrow(mat_pcs), nrow(mat_pcs))
  rownames(adj) <- colnames(adj) <- rownames(mat_pcs)
  for(i in seq_len(nrow(mat_pcs))) {
    adj[i,rownames(mat_pcs)[knn[i,]]] <- 1
  }
  
  #saveRDS(adj, "./sce_comp_adj.rds")
  
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
  
  saveRDS(df_umap, paste0(output_file, ".rds"))
  
  return(df_umap)
}

#' @name sc_plot_umap
#' @title Generates UMAP of data from sce object
#' @description Uses the output from sc_generate_umap_df() to produce a UMAP
#' @param df_umap the data frame generated from sc_gnerate_umap_df()
#' @param sce The SingleCellExperiment object
#' @param by What to colour the points by. Needs to be in colData of all experiments. If is left `NULL` then will be coloured by cluster.
#' @param output_file The path of the output file
#' @export
sc_plot_umap <- function(df_umap, sce,
                         by              = NULL,
                         custom_colours  = NULL,
                         output_file     = NULL){  
  if (!is.null(by) && !is.null(colData(sce)[[by]])) {
    #print("1 is true")
    sce_coldata <- colData(sce)[, c(by), drop=FALSE]
    df_umap[, colnames(sce_coldata)] <- data.frame(sce_coldata)
  } else if (!is.null(by)) {
    #print("2 is true")
    df_umap[[by]] <- NA
    cat("Couldn't locate", by, "in the column data of the provided SCE object\n")
  }
  if (!is.null(by) && !is.null(df_umap[[by]])) {
    #print("3 is true")
    g <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(col = .data[[by]]), size = 2, alpha = 0.5) +
      theme_bw(base_size = 14)
    if (is.numeric(df_umap[[by]])) {
      #print("4 is true")
      g <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(col = log(.data[[by]])), size = 2, alpha = 0.5) +
        scale_colour_gradientn(colours=c("green","black")) +
        theme_bw(base_size = 14)
    }
  } else {
    #print("5 is true")
    g <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(col = cluster), size = 2, alpha = 0.5) +
      theme_bw(base_size = 14)
  }
  if (!is.null(custom_colours)) {
    #print("6 is true")
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

#' @name plot_read_count_distr
#' @title compare the total reads from cellranger to scPipe
#' @param sce_cr cellranger output into sce
#' @param sce_scPipe scPipe output into sce
#' @param output_file out filename
#' @export
# Require CellRanger sce object first then scPipe sce object
plot_read_count_distr <- function(sce_cr, sce_scPipe, output_file = NULL) {
  sce <- attach_cell_comparison(sce_cr, sce_scPipe)
  sce_coldata <- data.frame(colData(sce))
  colnames(sce_coldata) <- colnames(colData(sce))
  g <- ggplot(sce_coldata, aes(y = log(total), x=`Overlap with scPipe`, fill=`Overlap with scPipe`)) +
    geom_boxplot() + 
    ylab("log(total read count)") +
    theme(axis.text.x=element_blank(),
          text = element_text(size = 18)) +
    guides(fill=guide_legend(title="Overlap \nwith scPipe"))
  if (!is.null(output_file)) ggsave(output_file)
  return(g)
}

#' @name plot_read_count_distr_rev
#' @title compare the total reads from cellranger to scPipe
#' @param sce_cr cellranger output into sce
#' @param sce_scPipe scPipe output into sce
#' @param output_file out filename
#' @export
# Require CellRanger sce object first then scPipe sce object
plot_read_count_distr_rev <- function(sce_scPipe, sce_cr, output_file = NULL) {
  sce <- attach_cell_comparison(sce_scPipe, sce_cr, "CellRanger")
  sce_coldata <- data.frame(colData(sce))
  colnames(sce_coldata) <- colnames(colData(sce))
  g <- ggplot(sce_coldata, aes(y = log(counts_per_cell), x=`Overlap with CellRanger`, fill=`Overlap with CellRanger`)) +
    geom_boxplot() + 
    ylab("log(total read count)") +
    theme(axis.text.x=element_blank(),
          text = element_text(size = 18)) +
    guides(fill=guide_legend(title="Overlap with \nCellRanger"))
  if (!is.null(output_file)) ggsave(output_file)
  return(g)
}
