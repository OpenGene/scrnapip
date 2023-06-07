#' @title integrated CytoTRACE (iCytoTRACE)
#'
#' @description This function generates single-cell predictions of differentiation status
#' across multiple, heterogeneous scRNA-seq batches/datasets. This implementation leverages
#' mutual nearest neighbor and Gaussian kernel normalization techniques from Scanorama to merge datasets.
#' It takes in a list of gene expression matrices where columns are cells and rows are genes and outputs
#' integrated gene counts, gene counts signature (GCS), CytoTRACE, and the correct gene expression matrix.
#'
#'
#' @param datasets list of gene expression matrices where columns are cells and rows are genes
#'
#' @return a list containing
#' \itemize{
#' \item exprMatrix: a matrix (genes by single cells) of the corrected gene expression values from Scanorama
#' \item CytoTRACE: a numeric vector of the predicted ordering of single cells from 1.0 (least differentiated) to 0.0 (most differentiated)
#' \item CytoTRACErank: a numeric vector of the ranked predicted ordering of single cells by differentiation status. High ranks correspond to less differentiated cells, while low ranks correspond to more differentiated cells.
#' \item GCS: a numeric vector of the merged gene counts signature (geometric mean of the top 200 genes associated with gene counts)
#' \item Counts: a numeric vector of the corrected number of genes expressed per single cell (gene counts)
#' \item coord: a matrix containing the coordinates for the merged low-dimensional embedding (number of cells x 100 t-SNE components)
#' \item filteredCells = a character vector of the names of single cells (columns) that were filtered due to poor quality.
#' }
#'
#' @author Gunsagar Gulati <cytotrace@gmail.com>
#'
#' @seealso https://cytotrace.stanford.edu
#'
#' @references https://doi.org/10.1101/649848
#'
#' @examples
#'
#' #Create a list containing two bone marrow scRNA-seq datasets profiled on different platforms, 10x and Smart-seq2
#' datasets <- list(marrow_10x_expr, marrow_plate_expr)
#'
#' #Run iCytoTRACE
#' results <- iCytoTRACE(datasets)
#' @export



iCytoTRACE <- function(datasets, enableFast = TRUE,
                       ncores = 1,subsamplesize = 1000) {
  if(!have_scanoramaCT | !have_numpy){
    stop("The necessary python modules are not accessible. This function is disabled. Please follow the instructions in https://github.com/gunsagargulati/CytoTRACE to install the Python packages for this application.")
  }
  #import python libraries
  suppressWarnings(scanoramaCT<- reticulate::import('scanoramaCT'))
  np <- reticulate::import('numpy')
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}

  if(class(datasets) != "list" | length(datasets) < 2){
    print("Please input a list of at least 2 gene by cell matrices")
  }

  processedMats <- lapply(datasets, function(mat){
    #Checkpoint: log2-normalization
    if(max(mat)<50){
      mat <- 2^mat - 1
    }

    #Checkpoint: ERCC standards
    if(length(grep("ERCC-", rownames(mat)))>0){
      mat <- mat[-grep("ERCC-", rownames(mat)),]
    }

    #Checkpoint: Sequencing depth normalization
    mat <- t(t(mat)/apply(mat, 2, sum))*1000000

    #Checkpoint: NAs and poor quality cells
    pqcells <- is.na(apply(mat>0, 2, sum)) | apply(mat>0, 2, sum) <= 10
    num_pqcells <- length(which(pqcells == TRUE))
    mat <- mat[,!pqcells]

    #Checkpoint: NAs and poor quality genes
    pqgenes <- is.na(rowSums(mat>0)) | apply(mat, 1, var) == 0
    num_pqgenes <- length(which(pqgenes == TRUE))
    mat <- mat[!pqgenes,]

    #Checkpoint: log2-normalize
    mat <- log(mat+1,2)
    mat <- data.matrix(mat)
    return(mat)
  }
  )
  #Calculate and min-max normalize gene counts
  countsOrig <- lapply(processedMats, function(x) colSums(x>0))
  countsNorm <- lapply(countsOrig, range01)
  countsNorm <- lapply(countsNorm, function(x) cbind(x, x))

  #Run modified scanoramaCT to batch correct matrix and gene counts
  #Create genes_list
  genes_list <- lapply(processedMats, rownames)
  cell_list <- unlist(lapply(processedMats, colnames))
  #transpose matrices
  processedMats <- lapply(processedMats, t)

  merge <- scanoramaCT$merge_datasets(processedMats, genes_list)
  process <- scanoramaCT$process_data(merge[[1]], merge[[2]], hvg = 0L, dimred = 100L)
  gcalign <- scanoramaCT$find_alignments_gc(process[[1]], countsNorm)
  integrated.corrected.data <- scanoramaCT$correct(processedMats, genes_list, return_dimred=TRUE, return_dense = T)
  countsNormCorrected <- unlist(lapply(gcalign, function(x) x[,1]))
  matCorrected <- t(do.call(rbind, integrated.corrected.data[[2]]))
  rownames(matCorrected) <- integrated.corrected.data[[3]]
  mat2 <- matCorrected

  #Getting plotting coordinates
  m <- do.call(rbind, integrated.corrected.data[[1]])
  colnames(mat2) <- rownames(m) <- cell_list

  #Function to identify the most variable genes
  mvg <- function(matn) {
    A <- matn
    n_expr <- rowSums(A > 0);
    A_filt <- A[n_expr >= 0.05 * ncol(A),];
    vars <- apply(A_filt, 1, var);
    means <- apply(A_filt, 1, mean);
    disp <- vars / means;
    last_disp <- tail(sort(disp), 1000)[1];
    A_filt <- A_filt[disp >= last_disp,];

    return(A_filt)
  }


  #Calculate similarity matrix
  similarity_matrix_cleaned <- function(similarity_matrix){
    D <- similarity_matrix
    cutoff <- mean(as.vector(D))
    diag(D) <- 0;
    D[which(D < 0)] <- 0;
    D[which(D <= cutoff)] <- 0;
    Ds <- D
    D <- D / rowSums(D);
    D[which(rowSums(Ds)==0),] <- 0
    return(D)
  }

  if(ncol(mat2) < 10000){
    enableFast = FALSE
    message("The number of cells in your integrated dataset is less than 10,000. Fast mode has been disabled.")
  } else {
    message("The number of cells in your dataset exceeds 10,000. CytoTRACE will now be run in fast mode (see documentation). You can multi-thread this run using the 'ncores' flag. To disable fast mode, please indicate 'enableFast = FALSE'.")
  }

  #Subsample routine
  if(enableFast == FALSE){
    size <- ncol(mat2)
  } else if (enableFast == TRUE & subsamplesize < ncol(mat2)){
    size <- subsamplesize
  } else if (enableFast == TRUE & subsamplesize >= ncol(mat2)){
    stop("Please choose a subsample size less than the number of cells in dataset.")
  }

  chunk <- round(ncol(mat2)/size)
  subsamples <- split(1:ncol(mat2), sample(factor(1:ncol(mat2) %% chunk)))
  message(paste("CytoTRACE will be run on", chunk, "sub-sample(s) of approximately",
                round(mean(unlist(lapply(subsamples, length)))), "cells each using", min(chunk, ncores),"/", ncores, "core(s)"))

  batches <- parallel::mclapply(subsamples, mc.cores = min(chunk, ncores), function(subsample){
    #Filter out cells not expressing any of the 1000 most variable genes
    mat2 <- mat2[,subsample]
    mat2.mvg <- mvg(mat2)
    rm1 <- colSums(mat2) == 0
    mat2 <- mat2[, !rm1]
    countsNormCorrected <- countsNormCorrected[subsample]
    countsNormCorrected <- countsNormCorrected[!rm1]
    D <- similarity_matrix_cleaned(HiClimR::fastCor(mvg(mat2)))
  return(list(mat2 = mat2,countsNormCorrected = countsNormCorrected, D = D))
  }
  )
  mat2 <- do.call(cbind, lapply(batches, function(x) x$mat2))
  countsNormCorrected <- do.call(c, lapply(batches, function(x) x$countsNormCorrected))
  D2 <- lapply(batches, function(x) x$D)

  #Calculate gene counts signature (GCS) or the genes most correlated with gene counts
  ds2 <- sapply(1:nrow(mat2), function(x) ccaPP::corPearson(mat2[x,],countsNormCorrected))
  names(ds2) <- rownames(mat2)
  gcs <- apply(mat2[which(rownames(mat2) %in% names(rev(sort(ds2))[1:200])),],2,mean)

  samplesize <- unlist(lapply(D2, ncol))
  gcs2 <- split(gcs, as.numeric(rep(names(samplesize), samplesize)))

  #Regress gene counts signature (GCS) onto similarity matrix
  regressed <- function(similarity_matrix_cleaned, score){
    out <- nnls::nnls(similarity_matrix_cleaned,score)
    score_regressed <- similarity_matrix_cleaned %*% out$x
    return(score_regressed)
  }

  #Apply diffusion to regressed GCS using similarity matrix
  diffused <- function(similarity_matrix_cleaned, score, ALPHA = 0.9){
    vals <- score
    v_prev <- rep(vals);
    v_curr <- rep(vals);

    for(i in 1:10000) {
      v_prev <- rep(v_curr);
      v_curr <- ALPHA * (similarity_matrix_cleaned %*% v_curr) + (1 - ALPHA) * vals;

      diff <- mean(abs(v_curr - v_prev));
      if(diff <= 1e-6) {
        break;
      }
    }
    return(v_curr)
  }

  cytotrace <- parallel::mclapply(1:length(D2), mc.cores = ncores, function(i) {
    gcs_regressed <- regressed(D2[[i]], gcs2[[i]])
    gcs_diffused <- diffused(D2[[i]], gcs_regressed)
    cytotrace <- rank(gcs_diffused)
  }
  )
  cytotrace <- cytotrace_ranked <- unlist(cytotrace)
  cytotrace <- range01(cytotrace)
  names(cytotrace) <- names(cytotrace_ranked) <- names(gcs) <- names(countsNormCorrected) <- colnames(mat2)

  #Calculate genes associated with iCytoTRACE
  cytogenes <- sapply(1:nrow(mat2),
                      function(x) ccaPP::corPearson(mat2[x,], cytotrace))
  names(cytogenes) <- rownames(mat2)
  message("Calculating genes associated with iCytoTRACE...")


  #Subset plotting coordinates
  m <- m[colnames(mat2),]

  #filter
  filter <- unlist(lapply(datasets, colnames))[-which(unlist(lapply(datasets, colnames)) %in% colnames(mat2))]

  #Final steps
  mat2 <- t(data.frame(t(mat2))[unlist(lapply(datasets, colnames)),])
  cytotrace <- cytotrace[unlist(lapply(datasets, colnames))]
  cytotrace_ranked <- cytotrace_ranked[unlist(lapply(datasets, colnames))]
  gcs <- gcs[unlist(lapply(datasets, colnames))]
  countsNormCorrected <- countsNormCorrected[unlist(lapply(datasets, colnames))]

  colnames(mat2) <- names(cytotrace) <- names(cytotrace_ranked) <- names(gcs) <- names(countsNormCorrected) <- unlist(lapply(datasets, colnames))
  return(list(exprMatrix = mat2, CytoTRACE = cytotrace, CytoTRACErank = cytotrace_ranked,
              GCS= gcs, Counts = countsNormCorrected, cytoGenes = cytogenes,
              coord = m, filteredCells = filter))
}

