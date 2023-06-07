#devtools::install_local("./CytoTRACE_0.3.3.tar.gz")
#library(reticulate)
#conda_create("cytoTRACE",python_version = '3.7')
#use_condaenv("cytoTRACE")
#conda_install("cytoTRACE", "numpy")
#reticulate::py_install(packages = c("scanoramaCT"))
library(reticulate)
#use_python('/opt/conda/bin/python',required = T)
library(CytoTRACE)
library(Seurat)
library(dplyr)
library(ggplot2)
library(getopt)
command = matrix(c("rds","r",1,"character",#输入rds文件
                   "outpath","o",1,"character",#输出路径
                   "type","t",2,"character",
                   "barcolor","b",2,"character",
                   "featurecolor","f",2,"character",
                   "ncore","n",2,"numeric",
                   "runfast","F",0,"logical",
                   "prefix","h",2,"character")
                 ,byrow = T,ncol = 4)

print("test1")

args = getopt(command)
print("test")
#输入rds文件
seuset = readRDS(args$rds)

#输出文件目录
output <- args$outpath
if(!dir.exists(output)){
  dir.create(output)
}
if(!is.null(args$prefix)){prefix <- args$prefix}else{prefix <- "CytoTRACE"}
if(!is.null(args$ncore)){core <- args$ncore}else{core <- 4}


print("[Start]")
#开始CytoTRACE分析
matrix = as.matrix(seuset@assays$RNA@counts)

if(is.null(args$runfast)){
  results <- CytoTRACE(mat=matrix,ncores=core,enableFast = FALSE)
}else{
  results <- CytoTRACE(mat=matrix,ncores=core)
}

plotCytoGenes(results, numOfGenes = 10,outputDir = paste0(output,"/",prefix))

#提取打分结果，绘制umap图
cyto_value = as.data.frame(results$CytoTRACE)
colnames(cyto_value) = "CytoTRACE_value"
seuset = AddMetaData(object = seuset, metadata = cyto_value)
if(!is.null(args$featurecolor)){
  fcolor = strsplit(args$featurecolor,split=",")[[1]]
  p = FeaturePlot(seuset,"CytoTRACE_value",raster=F, reduction = "umap",cols = fcolor)
}else{
  p = FeaturePlot(seuset,"CytoTRACE_value",raster=F, reduction = "umap")
}
ggsave(paste0(output,"/",prefix,".cytovalue.FeaturePlot.png"),p,width = 8.5,height = 8)
ggsave(paste0(output,"/",prefix,".cytovalue.FeaturePlot.pdf"),p,width = 8.5,height = 8)

#添加注释的cluster信息
seu_ident = as.data.frame(seuset@active.ident)
colnames(seu_ident)="Cluster"
seu_ident$Cluster = as.character(seu_ident$Cluster)
      
cyto_ident = merge(cyto_value,seu_ident,by="row.names")
#确定cluster的排列顺序
cyto_order = group_by(cyto_ident,Cluster) %>% summarise(cyto_median = median(CytoTRACE_value))
cyto_order = cyto_order[rev(order(cyto_order$cyto_median)),]
cyto_ident$Cluster=factor(cyto_ident$Cluster,levels = cyto_order$Cluster)

#开始绘制基础图形
p = ggplot(cyto_ident,aes(x=Cluster,y=CytoTRACE_value,fill=Cluster,color=Cluster))+
  geom_boxplot(alpha=0.5)+
  geom_jitter(width = 0.1,aes(color=Cluster),cex=0.5)+labs(x="Cell phenotypes",y="Predicted ordering by CytoTRACE")+
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1),legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=10),axis.title = element_text(size=20,face = "bold"))
if(!is.null(args$barcolor)){
  mycolor = read.table(args$barcolor,header = T,sep = "\t",comment.char = "")
  colnames(mycolor)=c("Cluster","color")
  mycolor$Cluster = factor(mycolor$Cluster,levels = cyto_order$Cluster)
  mycolor = mycolor[order(mycolor$Cluster),]
  p = p +  scale_fill_manual(values = mycolor$color)+
    scale_color_manual(values = mycolor$color)
}
ggsave(paste0(output,"/",prefix,".CytoTRACE.boxplot_raw.png"),p,width = length(unique(cyto_ident$Cluster)),height = 8)
ggsave(paste0(output,"/",prefix,".CytoTRACE.boxplot_raw.pdf"),p,width = length(unique(cyto_ident$Cluster)),height = 8)

#如果输入了类型文件，则按类型绘图
if(!is.null(args$type)){
  clustertype = read.table(args$type,header = T,sep = "\t",comment.char = "")
  colnames(clustertype)[c(1,2)]=c("Cluster","Celltype")
  cyto_ident$Cluster = as.character(cyto_ident$Cluster)
  clustertype$Cluster = as.character(clustertype$Cluster)
  cyto_ident = left_join(cyto_ident,clustertype,by="Cluster")

  cyto_order = group_by(cyto_ident,Celltype) %>% summarise(cyto_median = median(CytoTRACE_value))
  cyto_order = cyto_order[rev(order(cyto_order$cyto_median)),]
  cyto_ident$Celltype=factor(cyto_ident$Celltype,levels = cyto_order$Celltype)
  cyto_ident = cyto_ident[order(cyto_ident$Celltype),]
  p = ggplot(cyto_ident,aes(x=Celltype,y=CytoTRACE_value,fill=Celltype,color=Celltype))+
    geom_boxplot(alpha=0.5)+
    geom_jitter(width = 0.1,aes(color=Celltype),cex=0.5)+labs(x="Cell phenotypes",y="Predicted ordering by CytoTRACE")+
    theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1),legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(size=15),axis.title = element_text(size=20,face = "bold"))
  if(ncol(clustertype)==3){
    cyto_ident=cyto_ident[order(cyto_ident$Celltype),]
    p = p +  scale_fill_manual(values = unique(cyto_ident[,5]))+
      scale_color_manual(values = unique(cyto_ident[,5]))
  }
  ggsave(paste0(output,"/",prefix,".CytoTRACE.boxplot_type.png"),p,width = length(unique(cyto_ident$Cluster)),height = 8)
  ggsave(paste0(output,"/",prefix,".CytoTRACE.boxplot_type.pdf"),p,width = length(unique(cyto_ident$Cluster)),height = 8)
}

write.table(cyto_ident,paste0(output,"/",prefix,".CytoTRACE.table.txt"),col.names = T,row.names = F,sep = "\t",quote = F)
