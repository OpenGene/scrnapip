#' Mouse bone marrow differentiation scRNA-seq (10x) expression dataset
#'
#' Single-cell RNA-sequencing data from Tabula Muris of mouse bone marrow cells
#' profiled using the 10x platform.
#'
#' @docType data
#'
#' @usage marrow_10x_expr
#'
#' @format A data frame where single cells are columns and genes are rows.
#' Values represent unique molecular identifiers (UMIs).
#'
#' @keywords bone marrow, single-cell RNA-sequencing, 10x, gene expression
#'
#' @references Schaum, Nicholas, et al. "Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris: The Tabula Muris Consortium." Nature 562.7727 (2018): 367.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/30283141}{PubMed})
#'
#' @source GSE109774 and \href{https://github.com/czbiohub/tabula-muris}{Tabula Muris GitHub}
#'
"marrow_10x_expr"

#' Mouse bone marrow differentiation scRNA-seq (10x) expression dataset
#'
#' Single-cell RNA-sequencing data from Tabula Muris of mouse bone marrow cells
#' profiled using the Smart-seq2 platform.
#'
#' @docType data
#'
#' @usage marrow_plate_expr
#'
#' @format A data frame where single cells are columns and genes are rows.
#' Values represent raw counts.
#'
#' @keywords bone marrow, single-cell RNA-sequencing, Smart-seq2, gene expression
#'
#' @references Schaum, Nicholas, et al. "Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris: The Tabula Muris Consortium." Nature 562.7727 (2018): 367.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/30283141}{PubMed})
#'
#' @source GSE109774 and \href{https://github.com/czbiohub/tabula-muris}{Tabula Muris GitHub}
#'
"marrow_plate_expr"

#' Mouse bone marrow differentiation scRNA-seq (10x) phenotype labels
#'
#' Phenotype labels of mouse bone marrow cells from Tabula Muris corresponding to the data file 'marrow_10x_expr'.
#'
#' @docType data
#'
#' @usage marrow_10x_pheno
#'
#' @format A character vector of the phenotype labels of mouse bone marrow cells with names corresponding to the single cell IDs
#'
#' @keywords bone marrow, single-cell RNA-sequencing, 10x, phenotype
#'
#' @references Schaum, Nicholas, et al. "Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris: The Tabula Muris Consortium." Nature 562.7727 (2018): 367.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/30283141}{PubMed})
#'
#' @source GSE109774 and \href{https://github.com/czbiohub/tabula-muris}{Tabula Muris GitHub}
#'
"marrow_10x_pheno"


#' Mouse bone marrow differentiation scRNA-seq (Smart-seq2) phenotype labels
#'
#' Phenotype labels of mouse bone marrow cells from Tabula Muris corresponding to the data file 'marrow_10x_expr'.
#'
#' @docType data
#'
#' @usage marrow_plate_pheno
#'
#' @format A character vector of the phenotype labels of mouse bone marrow cells with names corresponding to the single cell IDs
#'
#' @keywords bone marrow, single-cell RNA-sequencing, Smart-seq2, phenotype
#'
#' @references Schaum, Nicholas, et al. "Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris: The Tabula Muris Consortium." Nature 562.7727 (2018): 367.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/30283141}{PubMed})
#'
#' @source GSE109774 and \href{https://github.com/czbiohub/tabula-muris}{Tabula Muris GitHub}
#'
"marrow_plate_pheno"
