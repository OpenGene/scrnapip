# scrnapip
A Systematic and Dynamic Pipeline for Single-Cell RNA Sequencing Analysis


**Single cell biological information analysis process**

##### Haplox (2023/02/27)

### Introduction

The types, states, and interactions of cells in human tissues vary greatly. Single-cell transcriptome sequencing (scRNA-SEQ) is a new technique for high-throughput sequencing analysis of the transcriptome in single cell. Single-cell transcriptome sequencing can complement conventional transcriptome sequencing (mRNA-seq: Batch RNA sequencing, comparing the average expression values of genes in all cells of the cell population), revealing the expression situation of all genes in the all-cause group in single cell, including the identified tissue cell types, reflecting the cell heterogeneity between different samples and the tissue microenvironment, so that we can better understand the real state and correlation of each cell in a Bulk tissue. Presents a real and comprehensive cellular world. Currently, single-cell transcriptome sequencing is mostly used in complex multicellular systems such as tumor, developmental, neural, and immune microenvironments.

The purpose of this tool is connect the analysis of single-cell data into a complete process to accelerate the speed of analysis and contribute to the progress in this field.

## A. Environment set up

### 1. Required software and environment

```bash
#Download fastp
conda install -c bioconda fastp
#Download and decompress cellranger
wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1677511262&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2Nzc1MTEyNjJ9fX1dfQ__&Signature=iP4gsUveuUB0-WDVjtfGo8zqxW9vAa~RX9WLeHwjeH9~PsV3YAPsibucKGUvd9LbRhLZgFIkH9YbdMt0lTFXOL77DuLASDEJxdsYw1a1JbGGfr5ovy-xMDyWhsVuhz-zd~hmjRmTDvck0FwefwZWyUQXNGhcpBn-CVJQrAa5RMAId91NmFk8lKA1-ZBPat-fJ0Ivyq3gnXjkHpsNb3oJw2lgUDtTLVDtnuQo~LRuFLz026RaKgTAhI1YQfw7ug7N4shTtbxo275AWIi870fpizNd-K6npOrGyDAkuB87x8AigJBIpfs2Loi33g9hRXxgs4irjOov5C7c-F3wV7iTjA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar -xzvf cellranger-7.1.0.tar.gz
```

### 2. Install R packages

```R
#Enter the R environment and install the required R packages,If this fails, try to install through conda
install.packages(c('rmarkdown', 'BiocManager','stringi','circlize','logging','optparse','colourpicker', 'formattable', 'msigdbr','readr', 'shinydashboard','shinyWidgets','DT','configr','Seurat'))
BiocManager::install(c('SingleR','monocle','biomaRt','ggtree','GSVA','qvalue','clusterProfiler','org.Hs.eg.db','org.Mm.eg.db','ReactomePA'))
```



## B. Start Workflow

### 1. Set config file

All input files and parameters are set in this configuration file. The main Settings that need to be changed are the following：
``` bash
#####[fastp_cellrange]:RAW data path. The pair end data must be split into two files
S1.R1=["/usr/data/SAMPLE1_S1_L001_R1_001.fastq.gz"]
S1.R2=["/usr/data/SAMPLE1_S1_L001_R2_001.fastq.gz"]
#If a sample has more than one raw data,you can merge them before or add path split by ",":
S1.R1=["/usr/data/SAMPLE1.1_S1_L001_R1_001.fastq.gz","/usr/data/SAMPLE1.2_S1_L001_R1_001.fastq.gz"]
S1.R2=["/usr/data/SAMPLE1.1_S1_L001_R2_001.fastq.gz","/usr/data/SAMPLE1.2_S1_L001_R2_001.fastq.gz"]

#####[indata]:cellranger matrix file path
S1="/usr/workout/02.cellranger/S1/outs/filtered_feature_bc_matrix"

#####[outpath]:output path
outpath="/usr/workout"

#####[tempdata]rds file output path
tempdata="workout"

#####[run]:The analysis that needs to be done should set to true,for example:
fastp=true#run fastp

#####[fastp]:Configure the fastp path and parameters
fastppath="/usr/fastp"
longr=26#R1 length after trim
ncode=5#The maximum number of N-bases

#####[cellrangle]:Configure the cellranger path and parameters
dockerusr="1025:1025"#user id
dir="/user/name"#The folder which docker mount
ref="/user/refdata-gex-GRCh38-2020-A"#Reference genome path
cellrangerpath="/usr/cellranger-6.1.2/cellranger"#software path of cellranger
expectcell=10000#expect cell number
localcores=32#Number of threads
localmem=64#Memory size
include_introns="false"#Whether to analyze introns

#####[step1]:
nFeature_RNA=[200,5000]#The cells were filtered by feature, keeping cells that feature between 200 and 5,000 
percent_mt=[0,10]#The cells were filtered by percent of mitochondria, keeping cells that percent of mitochondria less than 10%
mttype="MT"#Mitochondrial type, MT for humans and mt for mice

#####[step2]:
kfilter=200#Minimum number of cells per sample
normethod="SCT"#The merge method, which uses SCT by default, can also use vst to simply group samples together
nFeature=3000#Genes for subsequent analysis

#####[step3]:
heatmapnumber=9#Number of heatmaps drawn for pca
elbowdims=100#The number of PCS shown in the elbow diagram
dims=30#Select the top 30 PCs for dimensionality reduction
reduction="umap"#tSNE or UMAP
clustercell=true#Whether you need to cluster cells
resolution=0.6#Set the resolution when clustering
algorithm=1#Cluster modular optimization algorithm (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm
singler="/singleRdata/test.rds"#singleR database position

#####[step4]:
clustermarkers=true#Whether marker genes of each cluster need to be found
min_pct=0.25#The minimum proportion of marker gene in the number of cells is 0.25 by default
findmarkers_testuse="wilcox"#The method of finding marker gene
difcluster.test.a=[0,1]#Find Differential gene.If you want to find differences between samples,change cluster to ident
difcluster.test.b=[5,6]#Test indicates the group name,a for case and b for control
difcluster.test.testuse="wilcox"#Inspection method

#####[step5]:
meanexpression=0.5#Select the appropriate gene to mark the state, intercept the condition, default is 0.5
genenum=50#Number of gene in differential analysis heat map
numclusters=4#The number of clusters in a cluster
pointid=1#The branching points used in BEAM analysis
BEAMnumclusters=4#Number of clusters in heat map clustering
BEAMgn=50#BEAM analyzes heat map gene count
BEAMgenelist=["S100A12", "ALOX5AP", "PAD14", "NRG1", "MCEMP1", "THBS1","testgene"]#BEAM analyzes specific gene names

#####[step6]:
circosbin="/bin/get_exp.r"#Extraction expression
circos_perl_bin="/bin/circos_plot.pl"#Plot circos
```

### 2. Filtered data by fastp and cellranger
This R script is used for data filtering and comparison quantitative analysis, and related parameters are set in the configuration file config.ini.
```bash
Rscript fastp_cellranger.r -i config.ini
```

### 3. Seurat analysis
This R script is used for all advanced analysis, and related parameters are set in the configuration file config.ini.
```bash
Rscript singlecell.r -i config.ini
```

## C. Result

### 1. fastp
Sequence statistics and reads filtering result files were performed on the original data.

``` bash
── 01.fastp/
     └── <SampleName>/                                  <- config for report
            ├── <SampleName>_fastp.html                 <- Report generated by fastp
            ├── <SampleName>_fastp.json                 <- Statistical information generated by fastp
            ├── <SampleName>_S1_L001_R1_001.fastq.gz    <- R1 clean data
            └── <SampleName>_S1_L001_R2_001.fastq.gz    <- R2 clean data
```

![fastp_summary](readme_files/fastp.png)

Quality control summary statistics by fastp.

### 2. Cellranger

Cellranger results after mapping and quantitative.
``` bash
── 02.Cellranger/
     └── <SampleName>/               		
              └── outs/                                 <- Cellranger analysis results 
                    ├── analysis/                       <- Cluster by cellranger
                    ├── raw_feature_bc_matrix/          <- Unfiltered feature-barcode matrices MEX (usually not used)
                    ├── raw_feature_bc_matrix.h5        <- Unfiltered feature-barcode matrices HDF5
                    ├── filtered_feature_bc_matrix/     <- Filtered feature-barcode matrices MEX
                    ├── filtered_feature_bc_matrix.h5   <- Filtered feature-barcode matrices HDF5
                    ├── molecule_info.h5                <- Per-molecule read information
                    ├── metrics_summary.xls             <- Run summary information
                    ├── possorted_genome_bam.bam        <- Bam file of single cell alignment
                    ├── possorted_genome_bam.bam.bai    <- Bam index
                    └── web_summary.html                <- Run summary HTML
```

![RNAseq Workflow](readme_files/cellranger_web_report.png)

Cellranger web report.

### 3. CellFilter
Filter low-quality cells according to mitochondrial proportion and gene number
``` bash
── 03.CellFilter/
     ├── <SampleName>/
     │          ├── <SampleName>_countVfeature.png/pdf  <- Scatter plot of nCounts and nFeature
     │          ├── <SampleName>_countVmt.png/pdf       <- Scatter plot of nCounts and percent.mito
     │          ├── <SampleName>_libraryVmt.png/pdf     <- Scatter plot of nCounts and percent.mito,nCounts as the color of the dot   
     │          └── <SampleName>_voilin.png/pdf         <- Violin plot with dots of nCounts,nFeature and percent.mito
     └── summary.txt/                                   <- Summary statistical information of cell filtration
```

![filter](readme_files/S1_countVfeature.png)

Scatter plot of feature and UMIs for all cells.Filter the cells outside the two red lines

![filter](readme_files/S1_countVmt.png)
Scatter plot of  percent mitochondria and UMIs for all cells.Filter the cells above the red lines.

### 4. PCA_UMAP
Dimensionality reduction and clustering to filtered cells.Annotate cells with singleR.
``` bash
── 04.PCA_UMAP/
     ├── pca/
     │    ├── bowplot.png/pdf               <- Standard Deviation of top 50 PC
     │    ├── pca.png/pdf                   <- Sample distribution of pc1 and pc2
     │    └── pcaheatmap1.png/pdf           <- Heatmap of genes in a pc
     └── umap/
          ├── barplot.png/pdf               <- The proportion of each sample in each cluster
          ├── celltype.png/pdf              <- The cell type identified by singleR
          ├── cluster0.png/pdf              <- Pie charts of the percentage of cell types in each cluster
          ├── plotall_ident.png/pdf         <- Umap plot of samle
          ├── plotby_cluster.png/pdf        <- Umap plot of cluster
          ├── plotby_ident.png/pdf          <- Umap plot splited by sample 
          ├── plotby_nCount.png/pdf         <- Umap plot of nCount
          ├── singleR_celltype.xls          <- Statistical table of cell types determined by singleR
          └── summary.xls                   <- The number of cells per cluster in each sample
```
![cluster](readme_files/plotby_cluster.png)
The single cell expression matrix was reduced and clustering to obtain the final umap plot. The dots represent cells, and clusters split by colors.
![sample](readme_files/plotall_ident.png)
The distribution of samples in umap is shown. The dots represent cells, and samples split by colors .
![singleR](readme_files/celltype.png)
UMAP plot for cell type annotation.

### 5. MarkerGene
Find marker genes and display result by violin and feature umap.
``` bash
── 05.MarkerGene/
     ├── cluster* /
     │       ├── cluster_genelist.xls                    <- Marker genes of a cluster
     │       ├── cluster_genereduction.png/pdf           <- Feature plot of top 10 marker genes
     │       ├── cluster_genevlnplot.png/pdf             <- Violin plot of top 10 marker genes
     │       └── cluster_padj0.05_logFC0.5genelist.xls   <- Filtered marker genes by padj<0.05 and logFC>0.5
     ├── custer/
     │       └── Celltype
     │              ├── *_reduction.png/pdf              <- Feature plot of the genes entered in the config file
     │              └── *_vlnplot.png/pdf                <- Violin plot of the genes entered in the config file
     ├── cluster_top10geneheatmap.png/pdf                <- Heat maps of the top 10 marker genes in all clusters
     ├── clusterall_adj0.05_logFC0.5genelist.xls         <- All filtered marker genes of all clusters
     └── clusterall_top10genelist.xls                    <- Top 10 marker genes of all clusters
```
![marker_umap](readme_files/cluster0_genereduction.png)

Umap plot of top10 marker gene  for each cluster.

![marker_violin](readme_files/cluster0_genevlnplot.png)

 Violin plot of top 10 marker gene for each cluster.

![scRNASeq](readme_files/cluster_top10geneheatmap.png)

### 6. Pseudotime

To perform pseudotime analysis of the cells,we used monocle2  to select high discrete gene and draw trajectory diagram.By default, the first branch point is used for beam analysis.If you want to analyze other branch points,setting in the configuration file.

```bash
── 06.Pseudotime/
     ├── pseudotime /
     │    ├── geneplot.png/pdf                      <- Scatter plot for selecting high discrete gene
     │    ├── pseudotime_byclusters.png/pdf         <- Cell trajectory plot,colors represent different clusters
     │    ├── pseudotime_byPseudotime.png/pdf       <- Cell trajectory plot,colors represent the level of expression
     │    ├── pseudotime_bystate.png/pdf            <- Cell trajectory plot,colors represent different state
     │    ├── pseudotime_splitbyclusters.png/pdf    <- Cell trajectory plot split by clusters
     │    ├── pseudotime_splitbyorig.ident.png/pdf  <- Cell trajectory plot split by samples
     │    └── cell_state_clusters.xls               <- Cell state statistical table
     ├── BEAM /
     │    ├── custergene.png/pdf                    <- Gene expression in different branches of cells
     │    ├── top10gene.png/pdf                     <- Top 10 gene expression in different branches of cells
     │    ├── topgeneheatmap.pdf                    <- Heatmap for genes in different branches
     │    └── diffgenelist.xls                      <- p-value and q-value of genes
     ├── diffgenes /
     │    ├── topgeneheatmap.pdf                    <- Heatmap for genes in different branches
     │    └── cell_diffgene.xls                     <- p-value and q-value of genes
     └── QC.png/pdf
```

![pseudotime](readme_files/pseudotime_byclusters.png)

Cell trajectory plot drawed by monocle.

### 7. Cerebro

Cerebro(cell report browser), which allows users to interactively visualize various parts of single cell transcriptomics data without requiring bioinformatic expertise.Cerebro can draw various graphs to display single cell results like umap/tsne for 2D/3D,bar plot,violin plot,cluster tree etc.


![cerebro](readme_files/cerebro.png)

### 8. Circos

In order to display the marker gene more intuitively, we using the R package circlize for circos. Circos shows the marker genes of all clusters or cell types, and draws bubble, heatmap and scatter diagram on different tracks. 

![scRNASeq](readme_files/cluster_circos.png)



#### Citations:

[1]R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.

[2]Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34(17):i884-i890.

[3]Zheng GX, Terry JM, Belgrader P, et al. Massively parallel digital transcriptional profiling of single cells. Nat Commun. 2017;8:14049. Published 2017 Jan 16.

[4]Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data. Cell. 2021;184(13):3573-3587.

[5]Qiu X, Hill A, Packer J, Lin D, Ma YA, Trapnell C. Single-cell mRNA quantification and differential analysis with Census. Nat Methods. 2017;14(3):309-315. 

[6]Zhao J, Zhang S, Liu Y, et al. Single-cell RNA sequencing reveals the heterogeneity of liver-resident immune cells in human. Cell Discov. 2020;6:22. Published 2020 Apr 28.

[7]Abdi H, Williams LJ. Principal component analysis. Wiley Interdisciplinary Reviews: Computational Statistics. 2010;2:433–459. 

[8]McInnes, L., Healy, J., Saul, N. & Großberger, L. UMAP: uniform manifold approximation and projection. J. Open Source Softw. 3, 861 (2018).

[9]Van Der Maaten, L. & Hinton, G. Visualizing high-dimensional data using t-SNE. journal of machine learning research. J. Mach. Learn. Res. 9, 26 (2008).

[10]Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. Fast unfolding of communities in large networks. J. Stat. Mech. 2008, P10008 (2008).

[11]Aran D, Looney AP, Liu L, et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nat Immunol. 2019;20(2):163-172.

[12]Hillje R, Pelicci PG, Luzi L. Cerebro: interactive visualization of scRNA-seq data. Bioinformatics. 2020;36(7):2311-2313.

[13]Gu Z, Gu L, Eils R, Schlesner M, Brors B. circlize Implements and enhances circular visualization in R. Bioinformatics. 2014;30(19):2811-2812.

[14]Wickham H. (2016) ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag, New York. ISBN 978-3-319-24277-4
