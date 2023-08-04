**Single cell biological information analysis process**


### Introduction

The types, states, and interactions of cells in human tissues vary greatly. Single-cell transcriptome sequencing (scRNA-seq) is a new technique for high-throughput sequencing analysis of the transcriptome in single cell. Single-cell transcriptome sequencing can complement conventional transcriptome sequencing (mRNA-seq: Batch RNA sequencing, comparing the average expression values of genes in all cells of the cell population), revealing the expression situation of all genes in the all-cause group in single cell, including the identified tissue cell types, reflecting the cell heterogeneity between different samples and the tissue microenvironment, so that we can better understand the real state and correlation of each cell in a Bulk tissue. Presents a real and comprehensive cellular world. Currently, single-cell transcriptome sequencing is mostly used in complex multicellular systems such as tumor, developmental, neural, and immune microenvironments.

The purpose of this tool is connect the analysis of single-cell data into a complete process to accelerate the speed of analysis and contribute to the progress in this field.

## A. Environment set up

### 1. Download docker

```bash
docker pull zhangjing12/scrnapip
```

### 2. Download reference genomic

```bash
#Human reference (GRCh38) dataset required for Cell Ranger.
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
#Mouse reference dataset required for Cell Ranger.
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
```

### 3. Use docker

```bash
docker run -d -p 1921:8787 -p 1882:3838 -e PASSWORD=yourpassword -e USERID=youruserid -e GROUPID=yourgroupid -v /yourdatapath:/dockerpath zhangjing12/scrnapip
```

The image is created based on Rocker (https://rocker-project.org/images/versioned/rstudio.html). You can use the above command to access rstudio through port 8787, which is more convenient for users to use the process. The userid and groupid can be queried through the id command. For the port number, please confirm whether the corresponding port is open.



## B. Start Workflow

### 1. Set config file

All input files and parameters are set in this configuration file. The main settings that need to be changed are the following：

```bash
#####[fastp_cellrange]: RAW data path. The pair end data must be split into two files
S1.R1=["/usr/data/SAMPLE1_S1_L001_R1_001.fastq.gz"]
S1.R2=["/usr/data/SAMPLE1_S1_L001_R2_001.fastq.gz"]
#If a sample has more than one raw data,you can merge them before or add path split by ",":
S1.R1=["/usr/data/SAMPLE1.1_S1_L001_R1_001.fastq.gz","/usr/data/SAMPLE1.2_S1_L001_R1_001.fastq.gz"]
S1.R2=["/usr/data/SAMPLE1.1_S1_L001_R2_001.fastq.gz","/usr/data/SAMPLE1.2_S1_L001_R2_001.fastq.gz"]

#####[indata]: cellranger matrix file path
S1="/usr/workout/02.cellranger/S1/outs/filtered_feature_bc_matrix"

#####[outpath]: output path
outpath="/usr/workout"

#####[tempdata]rds file output path
tempdata="workout"

#####[run]: The analysis that needs to be done should set to true,for example:
fastp=true #run fastp

#####[fastp]: Configure the fastp path and parameters
fastppath="/home/bin/fastp"
longr=26 #R1 length after trim
ncode=5 #The maximum number of N-bases

#####[cellrangle]: Configure the cellranger path and parameters
dockerusr="1025:1025" #user id
dir="/user/name" #The folder which docker mount
ref="/user/refdata-gex-GRCh38-2020-A" #Reference genome path
cellrangerpath="/home/bin/cellranger-7.1.0/cellranger" #software path of cellranger
expectcell=10000 #expect cell number
localcores=32 #Number of threads
localmem=64 #Memory size
include_introns="false" #Whether to analyze introns

#####[step1]:
filetype="10x" #The format of the input file,could be "10x" or "csv"
csv_sep="" #Separator of the csv file
nFeature_RNA=[200,5000] #The cells were filtered by feature, keeping cells that feature between 200 and 5000 
percent_mt=[0,10] #The cells were filtered by percent of mitochondria, keeping cells that percent of mitochondria less than 10%
mttype="MT" #Mitochondrial type, MT for humans and mt for mice

#####[step2]:
kfilter=200 #Minimum number of cells per sample
normethod="SCT" #The merge method, which uses SCT by default, can also use vst to simply group samples together
nFeature=3000 #Genes for subsequent analysis

#####[step3]:
recluster=true #Whether to perform batch correction
mode="harmony" #Method of batch correction
heatmapnumber=9 #Number of heatmaps drawn for pca
elbowdims=100 #The number of PCS shown in the elbow diagram
dims=30 #Select the top 30 PCs for dimensionality reduction
reduction="umap" #tSNE or UMAP
clustercell=true #Whether you need to cluster cells
resolution=0.6 #Set the resolution when clustering
algorithm=1 #Cluster modular optimization algorithm (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm)
singler="/home/bin/singleRdata/singleRdata/test.rds" #singleR database position

#####[step4]:
clustermarkers=true #Whether marker genes of each cluster need to be found
min_pct=0.25 #The minimum proportion of marker gene in the number of cells is 0.25 by default
findmarkers_testuse="wilcox" #The method of finding marker gene
difcluster.test.a=[0,1] #Find Differential gene.If you want to find differences between samples,change cluster to ident
difcluster.test.b=[5,6] #Test indicates the group name, a for case and b for control
difcluster.test.testuse="wilcox" #Inspection method
ClusterProfiler=["true","Rscript","/home/bin/clusterProfiler.R","-a true -s org.Hs.eg.db,hsa,human -g 6 -t SYMBOL -d KEGG,BioCyc,PID,PANTHER,BIOCARTA -C 0.05"] #Enrichment analysis of difference analysis results. -a: Whether to use all background genes; -s: species; -g: The column of the gene in the file; -t: gene name type(SYMBOL,ENTREZID); -d: database name

#####[step5]:
meanexpression=0.5 #Select the appropriate gene to mark the state, intercept the condition, default is 0.5
genenum=50 #Number of gene in differential analysis heat map
numclusters=4 #The number of clusters in a cluster
pointid=1 #The branching points used in BEAM analysis
BEAMnumclusters=4 #Number of clusters in heat map clustering
BEAMgn=50 #BEAM analyzes heat map gene count
BEAMgenelist=["S100A12","ALOX5AP","PAD14","NRG1","MCEMP1","THBS1"] #BEAM analyzes specific gene names

#####[step6]:
circosbin="/home/bin/get_exp.r" #Extraction expression
circos_perl_bin="/home/bin/circos_plot.pl" #Plot circos

#####[step7]:
copykat_bin="/home/bin/copykat_v4.r" #Identify tumor cells

#####[step8]:
cytoTRACE_bin="/home/bin/cytotrace_230508.R" #Developmental potential analysis

#####[step9]:
genomicinstably_bin="/home/bin/genomicinstably.R" #Genomic instability analysist
org="human" #species

#####[step11]:
ClusterProfiler=["true","Rscript","/home/bin/clusterProfiler.R","-a true -s org.Hs.eg.db,hsa,human -g 1 -t SYMBOL -d KEGG,BioCyc,PID,PANTHER,BIOCARTA -C 0.05"] #Enrichment analysis of marker gene
```

### 2. Filtered data by fastp and cellranger

This R script is used for data filtering and comparison quantitative analysis, and relevant parameters are set in the configuration file config_Example.ini.

```bash
Rscript /home/bin/fastp_cellranger.r -i config_Example.ini
```

### 3. Seurat analysis

This R script is used for all advanced analyses, and relevant parameters are set in the configuration file config_Example.ini.

```bash
Rscript /home/bin/singlecell.r -i config_Example.ini
```

## C. Result

### 1. Fastp

Sequence statistics and reads filtering result files were performed on the original data.

```bash
── 01.Fastp/
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

```bash
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

Filter low-quality cells according to mitochondrial proportion and gene number.

```bash
── 03.CellFilter/
     ├── <SampleName>/
     │          ├── <SampleName>_countVfeature.png/pdf  <- Scatter plot of nCounts and nFeature
     │          ├── <SampleName>_countVmt.png/pdf       <- Scatter plot of nCounts and percent.mito
     │          ├── <SampleName>_libraryVmt.png/pdf     <- Scatter plot of nCounts and percent.mito,nCounts as the color of the dot   
     │          └── <SampleName>_voilin.png/pdf         <- Violin plot with dots of nCounts,nFeature and percent.mito
     └── summary.txt/                                   <- Summary statistical information of cell filtration
```

![filter](readme_files/S1_countVfeature.png)

Scatter plot of feature and UMIs for all cells.Filter the cells outside the two red lines.


![filter](readme_files/S1_countVmt.png)
Scatter plot of  percent mitochondria and UMIs for all cells.Filter the cells above the red lines.



### 4. PCA_UMAP

Dimensionality reduction and clustering to filtered cells.Annotate cells with singleR.

```bash
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
The distribution of samples in umap is shown. The dots represent cells, and samples split by colors.

![singleR](readme_files/celltype.png)
UMAP plot for cell type annotation.


### 5. MarkerGene

Find marker genes and display result by violin and feature umap.

```bash
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

Heatmap of top 10 marker gene for each cluster.



### 6. Pseudotime

To perform pseudotime analysis of the cells, we used monocle2 to select high discrete gene and draw trajectory diagram. By default, the first branch point is used for beam analysis.If you want to analyze other branch points, setting in the configuration file.

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

Cerebro(cell report browser), which allows users to interactively visualize various parts of single cell transcriptomics data without requiring bioinformatic expertise. Cerebro can draw various graphs to display single cell results like umap/tsne for 2D/3D, bar plot, violin plot, cluster tree etc.

![cerebro](readme_files/cerebro.png)

### 8. Copykat

Copykat is used to perform copy number analysis and predict tumor cells. The umap plot is used to display the results.

Copykat (Copynumber Karyotyping of Tumors) is a computational tool using integrative Bayesian approaches to identify genome-wide aneuploidy at 5MB resolution in single cells to separate tumor cells from normal cells, and tumor subclones using high-throughput sc-RNAseq data.

```bash
── 08.Copykat/
     ├── <SampleName1> /
     │    ├── <SampleName1>_copykat_clustering_results.rds            <- RDS file containing copykat results
     │    ├── <SampleName1>_copykat_CNA_raw_results_gene_by_cell.txt  <- Expression matrix sorted by absolute position of genes
     │    ├── <SampleName1>_copykat_CNA_results.txt                   <- Copy number variation of cells at each chromosomal location
     │    ├── <SampleName1>_copykat_heatmap.jpeg                      <- Heatmap of tumor cell prediction results
     │    ├── <SampleName1>_copykat.pdf                               <- Heatmap of tumor cell prediction results
     │    ├── <SampleName1>_copykat_prediction.txt                    <- Tumor cell predicted results of copycat
     │    ├── <SampleName1>_copykat_with_genes_heatmap.pdf            <- Heatmap of tumor cell prediction results with gene names
     │    ├── <SampleName1>_rawdata.txt                               <- Raw umi expression matrix
     │    ├── <SampleName1>.tumor_subtype.pdf                         <- Heatmap of prediction results of tumor subclonal types
     │    └── <SampleName1>.tumor_subtype.txt                         <- Tumor subclonal type
     ├── <SampleName2> /
     ├── ...
     │
     ├── copykat.counts.clusterbarplot.pdf/png                        <- Barplot of the cell number by clusters, arranged by tumor cells
     ├── copykat.counts.samplebarplot.pdf/png                         <- Barplot of the cell number by samples
     ├── copykat.prop.clusterbarplot.pdf/png                          <- Barplot of the proportion of cell by clusters
     ├── copykat.prop.samplebarplot.pdf/png                           <- Barplot of the proportion of cell by samples
     ├── copykat.umap.pdf/png                                         <- Umap plot of tumor cell prediction results
     ├── copykat.umap.splitsample.pdf/png                             <- Umap plot of tumor cell prediction results splited by samples
     └── Summary_copykat_prediction.txt                               <- Predicted results for all sample cells
```

<img src="readme_files/copykat_heatmap.png" title="" alt="copykat" data-align="center">

Heatmap of copykat prediction results.


### 9. CytoTRACE

CytoTRACE (Cellular (Cyto) Trajectory Reconstruction Analysis using gene Counts and Expression) is a computational method that predicts the differentiation state of cells from single-cell RNA-sequencing data. CytoTRACE leverages a simple, yet robust, determinant of developmental potential—the number of detectably expressed genes per cell, or gene counts. CytoTRACE have been validated on ~150K single-cell transcriptomes spanning 315 cell phenotypes, 52 lineages, 14 tissue types, 9 scRNA-seq platforms, and 5 species.

```bash
── 9.CytoTRACE/ 
       ├── <prefix>.CytoTRACE.boxplot_raw.pdf/png     <- Boxplot of cytotrace scores by cluster 
       ├── <prefix>.CytoTRACE.boxplot_type.pdf/png    <- Boxplot of cytotrace scores by cell type 
       ├── <prefix>.cytovalue.FeaturePlot.pdf/png     <- Umap plot characterized by cytotrace scores 
       ├── <prefix>CytoGenes.pdf                      <- Genes that most correlated with cytotrace score 
       └── <prefix>.CytoTRACE.table.txt               <- Table containing cytotrace scores and cell names
```

![cytotrace](readme_files/cytotest.CytoTRACE.boxplot_raw.png)

Boxplots ordered by median cytotrace score.


### 10. GenomicInstability

Genomic instability analysis (GIA) uses the aREA algorithm to quantitatively estimate the association between gene expression and chromosomal location by performing enrichment analysis of contiguously coded genes (loci-blocks) on the single cell gene expression profiles.

```bash
── 10.genomicInstability/
     ├── <prefix>_genomicinstability.pdf             <- Density plot of the gene instability score
     ├── <prefix>_genomicinstability.withanno.pdf    <- Density plot of the gene instability score with annotated information
     └── <prefix>_gis.result.txt                     <- Table containing genomic instability score and cell names
```

<img src="readme_files/genomicInstability.png" title="" alt="genomicInstability" data-align="center">

The genomic Instability score density plot.

### 11. CellChat

CellChat, a tool that is able to quantitatively infer and analyze intercellular communication networks from single-cell RNA-sequencing (scRNA-seq) data. CellChat predicts major signaling inputs and outputs for cells and how those cells and signals coordinate for functions using network analysis and pattern recognition approaches. Through manifold learning and quantitative contrasts, CellChat classifies signaling pathways and delineates conserved and context-specific pathways across different datasets. 

```bash
── 11.CellChat  /
     ├── bubble /
     │      └── <CellType>_bubble.pdf        <- A bubble plot showing the strength of gene communication between cell types
     ├── pathway /
     │      ├── <GeneName>_circle.pdf        <- Circle plot of gene signaling pathway network of each gene
     │      └── <GeneName>_hier.pdf          <- Hierarchy plot of gene signaling pathway network of each gene
     ├── <CellType>_circle.pdf               <- Circle plot of interacting gene numbers of each cell type
     └── topgeneheatmap.pdf                  <- Circle plot of gene signaling pathway newtwork
```

<img src="readme_files/cellchat.png" title="" alt="cellchat" data-align="center">

### 12. ClusterProfiler

Gene pathway enrichment analysis is to find a class of overexpressed genes in a set of genes. Here, based on five databases including BIOCARTA, BioCyc, GO, KEGG, and reactome, we perform enrichment analysis on the marker genes of the cluster respectively.

```bash
── 12.clusterProfiler/ 
           ├── <prefix>.BIOCARTA_Enrich.xls          <- Gene enrichment analysis results based on BIOCARTA database 
           ├── <prefix>.BioCyc.png/pdf               <- Bubble plot of gene enrichment analysis based on BIOCARTA database 
           ├── <prefix>.BioCyc_Enrich.xls            <- Gene enrichment analysis results based on BioCyc database 
           ├── <prefix>.BP.DAG.svg                   <- Directed acyclic graph of biological processes in gene ontology 
           ├── <prefix>.CC.DAG.svg                   <- Directed acyclic graph of cellular component in gene ontology 
           ├── <prefix>.MF.DAG.svg                   <- Directed acyclic graph of molecular function in gene ontology 
           ├── <prefix>.go.pdf/png                   <- Barplot with significantly enriched GO terms 
           ├── <prefix>.GO_Enrich.xls                <- Gene enrichment analysis results based on gene ontology database 
           ├── <prefix>.KEGG.png/pdf                 <- Bubble plot of gene enrichment analysis based on KEGG database 
           ├── <prefix>.KEGG_Enrich.xls              <- Gene enrichment analysis results based on KEGG database 
           ├── <prefix>.PANTHER_Enrich.xls           <- Gene enrichment analysis results based on PANTHER database 
           ├── <prefix>.PID_Enrich.xls               <- Gene enrichment analysis results based on PID database 
           ├── <prefix>.reactome.png/pdf             <- Bubble plot of gene enrichment analysis based on reactome database 
           └── <prefix>.Reactome_Enrich.xls          <- Gene enrichment analysis results based on reactome database
```

![clusterProfiler](readme_files/patho.go.png)

Barplot with significantly enriched GO terms.



![keggbubbleplot](readme_files/patho.KEGG.png)

The bubble plot for KEGG enriched analysis.




#### Citations:
[1] Griffiths, Jonathan A., Antonio Scialdone, John C. Marioni. 2018. “Using single-cell genomics to understand developmental processes and cell fate decisions.” Molecular Systems Biology 14: e8046. https://doi.org/https://doi.org/10.15252/msb.20178046

[2] scRNA-tools. https://www.scrna -tools.org/.

[3] R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.

[4] Chen, Shifu, Yanqing Zhou, Yaru Chen, Jia Gu. 2018. “fastp: an ultra-fast all-in-one FASTQ preprocessor.” Bioinformatics 34: i884-i890. https://doi.org/10.1093/bioinformatics/bty560

[5] Chen,Shifu.2023.“Ultrafastone-passFASTQdatapreprocessing,qualitycontrol,and deduplicationusingfastp.”Imeta2:e107.https://doi.org/https://doi.org/10.1002/imt2.107

[6] Zheng, Grace X. Y., Jessica M. Terry, Phillip Belgrader, Paul Ryvkin, Zachary W. Bent, Ryan Wilson, Solongo B. Ziraldo, et al. 2017. “Massively parallel digital transcriptional profiling of single cells.” Nature Communications 8: 14049. https://doi.org/10.1038/ncomms14049

[7] Hao, Yuhan, Stephanie Hao, Erica Andersen-Nissen, William M. Mauck, III, Shiwei Zheng, Andrew Butler, Maddie J. Lee, et al. 2021. “Integrated analysis of multimodal single-cell data.” Cell 184: 3573-3587.e3529. https://doi.org/10.1016/j.cell.2021.04.048

[8] Qiu, Xiaojie, Andrew Hill, Jonathan Packer, Dejun Lin, Yi-An Ma, Cole Trapnell. 2017. “Single-cell mRNA quantification and differential analysis with Census.” Nature Methods 14: 309-315. https://doi.org/10.1038/nmeth.4150 

[9] Dobin, Alexander, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras. 2013. “STAR: ultrafast universal RNA-seq aligner.” Bioinformatics 29: 15-21. https://doi.org/10.1093/bioinformatics/bts635

[10] McGinnis, Christopher S., Lyndsay M. Murrow, Zev J. Gartner. 2019. “DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors.” Cell Systems 8: 329-337.e324. https://doi.org/10.1016/j.cels.2019.03.003

[11] Anders, Simon, Wolfgang Huber. 2010. “Differential expression analysis for sequence count data.” Genome Biology 11: R106. https://doi.org/10.1186/gb-2010-11-10-r106

[12] Love, Michael I., Wolfgang Huber, Simon Anders. 2014. “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology 15: 550. https://doi.org/10.1186/s13059-014-0550-8

[13] Korsunsky, Ilya, Nghia Millard, Jean Fan, Kamil Slowikowski, Fan Zhang, Kevin Wei, Yuriy Baglaenko, Michael Brenner, Po-ru Loh, Soumya Raychaudhuri. 2019. “Fast, sensitive and accurate integration of single-cell data with Harmony.” Nature Methods 16: 1289-1296. https://doi.org/10.1038/s41592-019-0619-0

[14] Haghverdi, Laleh, Aaron T. L. Lun, Michael D. Morgan, John C. Marioni. 2018. “Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors.” Nature Biotechnology 36: 421-427. https://doi.org/10.1038/nbt.4091

[15] Tran, Hoa Thi Nhu, Kok Siong Ang, Marion Chevrier, Xiaomeng Zhang, Nicole Yee Shin Lee, Michelle Goh, Jinmiao Chen. 2020. “A benchmark of batch-effect correction methods for single-cell RNA sequencing data.” Genome Biology 21: 12. https://doi.org/10.1186/s13059-019-1850-9

[16] Abdi, Hervé, Lynne J. Williams. 2010. “Principal component analysis.” WIREs Computational Statistics 2: 433-459. https://doi.org/https://doi.org/10.1002/wics.101 

[17] McInnes, Leland, John Healy. 2018. “UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.” The Journal of Open Source Software 3: 861. https://doi.org/10.48550/arXiv.1802.03426

[18] Becht, Etienne, Leland McInnes, John Healy, Charles-Antoine Dutertre, Immanuel W. H. Kwok, Lai Guan Ng, Florent Ginhoux, Evan W. Newell. 2019. “Dimensionality reduction for visualizing single-cell data using UMAP.” Nature Biotechnology 37: 38-44. https://doi.org/10.1038/nbt.4314

[19] Laurens, Van Der Maaten, Geoffrey Hinton. 2008. “Visualizing Data using t-SNE.” Journal of Machine Learning Research 9: 2579-2605.

[20] Blondel, Vincent D., Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre. 2008. “Fast unfolding of communities in large networks.” Journal of Statistical Mechanics: Theory and Experiment 2008: P10008. https://doi.org/10.1088/1742-5468/2008/10/P10008

[21] Aran, Dvir, Agnieszka P. Looney, Leqian Liu, Esther Wu, Valerie Fong, Austin Hsu, Suzanna Chak, et al. 2019. “Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage.” Nature Immunology 20: 163-172. https://doi.org/10.1038/s41590-018-0276-y

[22] Hillje, Roman, Pier Giuseppe Pelicci, Lucilla Luzi. 2020. “Cerebro: interactive visualization of scRNA-seq data.” Bioinformatics 36: 2311-2313. https://doi.org/10.1093/bioinformatics/btz877

[23] Gu, Zuguang, Lei Gu, Roland Eils, Matthias Schlesner, Benedikt Brors. 2014. “circlize implements and enhances circular visualization in R.” Bioinformatics 30: 2811-2812. https://doi.org/10.1093/bioinformatics/btu393

[24] Wickham, Hadley. 2009. Ggplot2: Elegant Graphics for Data Analysis. https://doi.org/10.1007/978-0-387-98141-3

[25] Gao, Ruli, Shanshan Bai, Ying C. Henderson, Yiyun Lin, Aislyn Schalck, Yun Yan, Tapsi Kumar, et al. 2021. “Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes.” Nature Biotechnology 39: 599-608. https://doi.org/10.1038/s41587-020-00795-2

[26] Gulati, Gunsagar S., Shaheen S. Sikandar, Daniel J. Wesche, Anoop Manjunath, Anjan Bharadwaj, Mark J. Berger, Francisco Ilagan, et al. 2020. “Single-cell transcriptional diversity is a hallmark of developmental potential.” Science 367: 405-411. https://doi.org/10.1126/science.aax0249

[27] Jin, Suoqin, Christian F. Guerrero-Juarez, Lihua Zhang, Ivan Chang, Raul Ramos, Chen-Hsiang Kuan, Peggy Myung, Maksim V. Plikus, Qing Nie. 2021. “Inference and analysis of cell-cell communication using CellChat.” Nature Communications 12: 1088. https://doi.org/10.1038/s41467-021-21246-9

[28] Yu, Guangchuang, Li-Gen Wang, Yanyan Han, Qing-Yu He. 2012. “clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters.” OMICS: A Journal of Integrative Biology 16: 284-287. https://doi.org/10.1089/omi.2011.0118

[29] Zhao, Juanjuan, Shuye Zhang, Yang Liu, Xiaomeng He, Mengmeng Qu, Gang Xu, Hongbo Wang, et al. 2020. “Single-cell RNA sequencing reveals the heterogeneity of liver-resident immune cells in human.” Cell Discovery 6: 22. https://doi.org/10.1038/s41421-020-0157-z
