library(Seurat)
library(getopt)
library(stringr)
library(Matrix)
library(dplyr)
library(tidyr)
command = matrix(c('rds','r',1,"character",
                   'outpath','o',1,"character",
                   'celltype','c',2,"character",
                   'type2gene','g',2,"character",
                   'marker','m',2,"character",
                   'prefix','s',2,"character")
                 ,byrow=TRUE, ncol=4)
args=getopt(command)
options(stringsAsFactors = FALSE)

rdsFile = args$rds;
outpath = args$outpath;
if(is.null(args$prefix)){args$prefix = "RDS"}
prefix = args$prefix;
if(!dir.exists(outpath)) {
  dir.create(outpath)
}


seuset=readRDS(rdsFile)


if(!is.null(args$celltype)){
  type = read.table(args$celltype,sep="\t",header=T,comment.char = "")
  colnames(type)=c("Cell","type")
}else{
  type = data.frame(Cell=rownames(seuset@meta.data),type=seuset@meta.data$seurat_cluster)
}
sample = data.frame(Cell=rownames(seuset@meta.data),sample=seuset@meta.data$orig.ident)
celllist = merge(type,sample,by="Cell")
celllist$type_sample = paste0(celllist$type,"_",celllist$sample)
celllist = celllist[,c("Cell","type_sample")]
realcluster = celllist[is.finite(match(celllist[,1],colnames(seuset))),]
seuset = subset(seuset,cells=as.character(realcluster[,1]))
DefaultAssay(seuset) = "RNA"
data = GetAssayData(seuset, slot = "data")

if(!is.null(args$type2gene)){
  seuset = SetIdent(object = seuset, cells= as.character(type[,1]), as.character(type[,2]))
  markers <- FindAllMarkers(seuset,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
  allmarkers <- data.frame(Gene=markers[,7],markers[,1:6])
  write.table(allmarkers,paste0(outpath,"/allmarker_gene.txt"),quote = F,sep = "\t",row.names = F,col.names = T)
  
  type2gene = read.table(args$type2gene,sep="\t",header=F)
  genelist = type2gene[,2]
  genelist = genelist[is.finite(match(genelist,rownames(data)))]
}else{
  allmarkers = read.table(args$marker,header=T,sep="\t")
  genelist=NULL
  type2gene=data.frame(Type=c(),Gene=c())
  for(i in unique(allmarkers$cluster)){
    tmp_marker = allmarkers[allmarkers$cluster==i,]
    for(j in 1:nrow(tmp_marker)){
      if(!tmp_marker$Gene[j]%in%genelist && tmp_marker$Gene[j]%in%rownames(seuset@assays$RNA@data)){
        genelist = c(genelist,tmp_marker$Gene[j])
        type2gene=rbind(type2gene,data.frame(Type=i,Gene=tmp_marker$Gene[j]))
        print(paste0("cluster",i))
        genelist
        type2gene
        break
      }
    }
  }
  
  write.table(type2gene,paste0(outpath,"/type2gene.txt"),quote=F,sep="\t",row.names=F,col.names=T)
  
}
seuset = SetIdent(object = seuset, cells= as.character(realcluster[,1]), as.character(realcluster[,2]))

data = (data)[genelist,,drop=F]
nrow = nrow(data)
line = 1
meandata = c()
pctdata = c()
while( line <= nrow) {
  line2 = line+1000-1
  line2 = min(line2,nrow)
  print(line2)
  data1 = data[line:line2,,drop=F]
  newdata = aggregate(t(as.matrix(data1)), data.frame(seuset@active.ident), FUN=function(x) log(x = mean(x = expm1(x = x)) + 1))
  names = newdata[,1]
  newdata = t(newdata[,-1,drop=F])
  newdata = data.frame(Gene = rownames(newdata),newdata)
  colnames(newdata) = c("Gene",as.character(names))
  meandata = rbind(meandata,newdata)
  
  newdata = aggregate(t(as.matrix(data1)), data.frame(seuset@active.ident), FUN=function(x) {return(sum(x > 0) / length(x = x))})
  names = newdata[,1]
  newdata = t(newdata[,-1,drop=F])
  newdata = data.frame(Gene = rownames(newdata),newdata)
  colnames(newdata) = c("Gene",as.character(names))
  pctdata = rbind(pctdata,newdata)
  line = line2+1
}
meandata[mapply(is.infinite,meandata)]<-0
mean_dat = data.frame(type_sample=colnames(meandata)[-1],t(meandata[,-1]))
mean_dat <- mean_dat %>% separate(type_sample,c("type","sample"),"_")
pctdata[mapply(is.infinite,pctdata)]<-0
pct_dat = data.frame(type_sample=colnames(pctdata)[-1],t(pctdata[,-1]))
pct_dat <- pct_dat %>% separate(type_sample,c("type","sample"),"_")
colnames(mean_dat)=gsub("\\.","-",colnames(mean_dat))
colnames(pct_dat)=gsub("\\.","-",colnames(mean_dat))
write.table(mean_dat,file=paste0(outpath,"/",prefix,".ClusterMean.txt"),sep="\t",row.names=F,quote=F)
write.table(pct_dat,file=paste0(outpath,"/",prefix,".cellpercent.txt"),sep="\t",row.names=F,quote=F)
