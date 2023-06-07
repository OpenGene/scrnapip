<a href="https://cytotrace.stanford.edu"><p align="center"><img src="https://github.com/gunsagargulati/CytoTRACE/blob/master/man/figures/cytotrace_logo.png" alt="CytotRACE_logo" width=400></a></p>

## Overview
CytoTRACE (Cellular (**Cyto**) **T**rajectory **R**econstruction **A**nalysis using gene **C**ounts and **E**xpression) is a computational method that predicts the differentiation state of cells from single-cell RNA-sequencing data. 
CytoTRACE leverages a simple, yet robust, determinant of developmental potential—the number of detectably expressed genes per cell, or gene counts. 
CytoTRACE has been validated on ~150K single-cell transcriptomes spanning 315 cell phenotypes, 52 lineages, 14 tissue types, 9 scRNA-seq platforms, and 5 species.

### Manuscript pre-print
A pre-print of our manuscript describing the development and application of CytoTRACE, entitled **Single-cell transcriptional diversity is a hallmark of developmental potential** is now available on <a href="https://doi.org/10.1101/649848"><i>bioRxiv</i></a>.

### CytoTRACE webtool
Additional details and general purpose use of CytoTRACE can be accessed on our website, https://cytotrace.stanford.edu. Key features of the website include the ability to:
- Analyze 42 publically available, annotated scRNA-seq datasets pre-computed with CytoTRACE
- Predict differentiation states in a custom scRNA-seq dataset
- Predict differentiation states across multiple batches/datasets from different platforms and developmental stages
- Visualize predicted differentiation with an interative 3D graph (includes t-SNE, force directed layout, and UMAP)
- Summarize results by known phenotypes
- Identify predicted stemness- and differentiation-associated genes

**This GitHub page provides instructions for the installation and use of the CytoTRACE R package.**

---

## Installation

The latest development version of the CytoTRACE R package is maintained on GitHub. To install from GitHub, run:

```r
install.packages("devtools")
devtools::install_github("gunsagargulati/CytoTRACE")
```

The `iCytoTRACE()` function in the R package requires 2 Python packages, `scanoramaCT`, an adapted version of the original <a href="https://github.com/brianhie/scanorama">Scanorama code</a> for application to CytoTRACE, and `numpy`. The `cytoTRACE()` function will run without these dependencies, but to enable application of CytoTRACE across multiple batches/datasets, install the Python dependencies by running:

```shell
$ pip install scanoramaCT update
$ pip install numpy
```
The development version of `scanoramaCT` is also available on GitHub and can be installed by running the following:

```shell
$ git clone https://github.com/gunsagargulati/scanoramaCT.git
$ cd scanorama/
$ python setup.py install --user
```

**Troubleshooting**: 
- We use the `reticulate` package to call Python functions in R. If you have multiple versions of Python, make sure the Python used to install the above packages is used by `reticulate`. You can do this by running the following in R:

```python
Sys.setenv(RETICULATE_PYTHON="/PATH/TO/PYTHON/BIN"). 
```
- CytoTRACE was developed and tested on R >= v3.5.1 using Python 2/3. Other versions may be incompatible. 

---

## Running CytoTRACE

Load CytoTRACE in R with `library(CytoTRACE)`. The package contains the following contents:

- `cytotrace()`: function to run CytoTRACE on a custom scRNA-seq dataset
- `iCytoTRACE`: function to run CytoTRACE across multiple, heterogeneous scRNA-seq batches/dataset
- Two bone marrow differentiation scRNA-seq datasets (`marrow_10x_expr` and `marrow_plate_expr`) with corresponding phenotype labels (`marrow_10x_pheno` and `marrow_plate_pheno`)

### Example I: Run CytoTRACE on a custom scRNA-seq dataset

Use the bone marrow 10x scRNA-seq dataset to run CytoTRACE
```r
results <- cytoTRACE(marrow_10x_expr)
```

Run this dataset on fast mode using 8 cores
```r
results <- cytoTRACE(marrow_10x_expr, enableFast = TRUE, ncores = 8)
```

The ouput is a list object containing numeric values for CytoTRACE, GCS, and gene counts, a numeric vector of the Pearson correlation between each gene and gene counts, and the IDs of filtered cells (see package documentation for more details).  

### Example II: Run iCytoTRACE on multiple scRNA-seq batches/datasets

Run `iCytoTRACE` on a list containing two bone marrow scRNA-seq datasets profiled on different platforms, 10x and Smart-seq2

```r
datasets <- list(marrow_10x_expr, marrow_plate_expr)
results <- iCytoTRACE(datasets)
```
The ouput is a list object containing numeric values for the merged CytoTRACE, GCS, and gene counts, the Scanorama-corrected gene expression matrix, the merged low dimensional embedding, and the IDs of filtered cells (see package documentation for more details).  

---

## Questions, comments, and concerns:

We are constantly improving the software and welcome all feedback from the community. For all questions, comments, or concerns regarding the CytoTRACE code and manuscript, please contact us at cytotrace@gmail.com. 


