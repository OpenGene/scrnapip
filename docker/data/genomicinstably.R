library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(SingleCellExperiment)
library(ExperimentHub)
library(genomicInstability)
library(readr)
library(dplyr)
library(ggplot2)
library(getopt)
library(clusterProfiler)
command = matrix(c("rds","r",1,"character",
                   "outpath","o",1,"character",
                   "anno","a",2,"character",
                   "species","s",1,"character",
                   "prefix","h",2,"character")
                 ,byrow = T,ncol = 4)
args = getopt(command)
options(stringsAsFactors = FALSE)
output = args$outpath
if(!dir.exists(output)){
  dir.create(output)
}
#setwd("C:\\Users\\Administrator\\Desktop\\zl")
#seuset<- readRDS("step4.RDS")
#output = "C:\\Users\\Administrator\\Desktop\\zl\\genomic"
if(!is.null(args$prefix)){prefix= args$prefix}else{prefix="data"}
if(args$species=="human"){
  database = "org.Hs.eg.db"
}else if(args$species=="mouse"){
  database = "org.Mm.eg.db"
}else{
  stop("check the species,only for human or mouse")
}

#input rds file
seuset = readRDS(args$rds)
sc_matrix = as.matrix(seuset@assays$RNA@counts)

list = bitr(rownames(sc_matrix),fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = database)
if(TRUE %in% duplicated(list$SYMBOL)){
  list = list[-which(duplicated(list$SYMBOL)),]
}

#list = AnnotationDbi::select(org.Hs.eg.db,keys=rownames(sc_matrix),columns=c("ENTREZID"),keytype="SYMBOL")
#dim(list)
#list = filter(list,ENTREZID!="NA")

sc_matrix = sc_matrix[rownames(sc_matrix)%in%list$SYMBOL,]
list = list[match(rownames(sc_matrix),list$SYMBOL),]
identical(list$SYMBOL,rownames(sc_matrix))
rownames(sc_matrix)=list$ENTREZID


#genomic Instability analysis
cnv_result <- inferCNV(sc_matrix,species = args$species)

cnv_result <- genomicInstabilityScore(cnv_result)

cnv_result <- giLikelihood(cnv_result, recompute=FALSE, normal=1, tumor=2:3)


#Cells with likelihood less than 0.25 are defined as normal cells for further analysis
cnv_norm <- inferCNV(sc_matrix, nullmat=sc_matrix[, cnv_result$gi_likelihood<0.25, drop=FALSE],species=args$species)

if(!is.null(cnv_norm$gis)){
  cnv_norm <- genomicInstabilityScore(cnv_norm, likelihood=TRUE)
}else{
  cnv_norm=cnv_result
}


# Plotting the density distribution for GIS and fitted models
pdf(paste0(output,"/",prefix,"_genomicinstability.pdf"))
par(mai=c(0.8,0.8, 0.2, 0.8))
giDensityPlot(cnv_norm, ylim=c(0, 1.1))
# Adding the likelihood data and second-axis
pos <- order(cnv_norm$gis)
lines(cnv_norm$gis[pos], cnv_norm$gi_likelihood[pos], lwd=2, col="blue")
axis(4, seq(0, 1, length=6), seq(0, 1, 0.2), col="blue", col.axis="blue")
axis(4, 0.5, "Relative likelihood", tick=FALSE, line=1.5, col.axis="blue")
pos5 <- which.min((0.5-cnv_norm$gi_likelihood)^2)
lines(c(rep(cnv_norm$gis[pos5], 2), max(cnv_norm$gis*1.05)),
      + c(0, rep(cnv_norm$gi_likelihood[pos5], 2)), lty=3, col="blue")
dev.off()

nes_result = as.data.frame(cnv_norm$nes)
gis = as.data.frame(cnv_norm$gis)
likehood = as.data.frame(cnv_norm$gi_likelihood)
result = merge(gis,likehood,by="row.names")
colnames(result)=c("Cells","gis","gi_likelihood")
write.table(result,file = paste0(output,"/",prefix,"_gis.result.txt"),col.names = T,row.names = F,sep = "\t",quote = F)




if(!is.null(args$anno)){
  
  metadata= read.table(args$anno,sep = "\t",header = T,row.names = 1)
  metadata$type = ifelse(metadata$copykat.pred=="aneuploid",1,ifelse(metadata$copykat.pred=="diploid",0,"other"))
  metadata = filter(metadata,type!="other")
  sc_matrix = sc_matrix[,which(colnames(sc_matrix) %in% rownames(metadata))]
  metadata_tumor = as.vector(metadata$type)[match(colnames(sc_matrix),rownames(metadata))]
  metadata_tumor = as.logical(as.numeric(metadata_tumor))
  
  pdf(paste0(output,"/",prefix,"_genomicinstability.withanno.pdf"))
  par(mai=c(0.8, 0.8, 0.2, 0.8))
  giDensityPlot(cnv_norm, ylim=c(0, 1.1))
  # Estimating GIS density distributions for normal and tumor cells
  gis_normal <- cnv_norm$gis[!metadata_tumor]
  gis_tumor <- cnv_norm$gis[metadata_tumor]
  den_normal <- density(gis_normal, from=min(gis_normal), to=max(gis_normal))
  den_tumor <- density(gis_tumor, from=min(gis_tumor), to=max(gis_tumor))
  
  # Scaling the densities based on the inferred proportions for each group
  den_normal$y <- den_normal$y * cnv_norm$gi_fit$lambda[1]
  den_tumor$y <- den_tumor$y * sum(cnv_norm$gi_fit$lambda[2])
  # Function to add the density plot
  addDensity <- function(x, col="grey") {
    polygon(c(x$x[1], x$x, x$x[length(x$x)], x$x[1]), c(0, x$y, 0, 0), col=col)
  }
  # Adding the density plots
  addDensity(den_normal, col=hsv(0.6, 0.8, 0.8, 0.4))
  addDensity(den_tumor, col=hsv(0.05, 0.8, 0.8, 0.4))
  # Adding the likelihood data and second-axis
  pos <- order(cnv_norm$gis)
  lines(cnv_norm$gis[pos], cnv_norm$gi_likelihood[pos], lwd=2, col="darkgreen")
  axis(4, seq(0, 1, length=6), seq(0, 1, 0.2), col="darkgreen", col.axis="darkgreen")
  axis(4, 0.5, "Relative likelihood", tick=FALSE, line=1.5, col.axis="darkgreen")
  pos5 <- which.min((0.5-cnv_norm$gi_likelihood)^2)
  lines(c(rep(cnv_norm$gis[pos5], 2), max(cnv_norm$gis*1.05)),
        + c(0, rep(cnv_norm$gi_likelihood[pos5], 2)), lty=3, col="blue")
  # Adding the legend
  legend(c(0.9, 1), c("Normal", "Tumor"), fill=hsv(c(0.6, 0.05), 0.8, 0.8, 0.4), bty="n")
  dev.off()
}