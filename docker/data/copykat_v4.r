library(foreach)  
library(doParallel)   
library(Seurat)
library(copykat)
library(getopt)
library(devtools)
library(ggplot2)
library(ggrepel)
library(dplyr)
command = matrix(c("rds","r",1,"character",#rds文件
                   "outpath","o",1,"character",#输出路径
                   "ngenechr","n",2,"numeric",#每个染色体最小基因数量
                   "subtype","t",2,"numeric",#肿瘤亚克隆数量
                   "specie","e",2,"character",#物种,human或mouse
                   "ncore","c",2,"numeric",#核心数量
                   "ncol","l",2,"numeric",#按样本绘图时每列的个数
                   "color","C",2,"character")#柱状图的颜色
                 ,byrow = T,ncol = 4)
args = getopt(command)
options(stringsAsFactors = FALSE)

rdsfile <- args$rds
output <- args$outpath
ngene_chr <- 5
subn <- 2
ncore <- 8
species <- "hg20"
ncol <- 3
colors <- c("#c82d31,#194f97,#00686b")
if(!is.null(args$ngenechr)){ngene_chr <- args$ngenechr}
if(!is.null(args$subtype)){subn <- args$subtype}
if(!is.null(args$specie) && args$specie%in%c("MOUSE","mouse","m")){
  species<-"mm10"
}
if(!is.null(args$ncore)){ncore <- args$ncore}
if(!is.null(args$ncol)){ncol = args$ncol}
if(!is.null(args$color)){colors = args$color}
mycol = strsplit(colors,split=",")[[1]]
if(length(mycol)>=3){
  mycol = mycol[1:3]
}else{
  stop("[ERROR]This program needs at least a color vector of length 3")
}

print("Start")

raw <- readRDS(rdsfile)
samplenum <- length(unique(raw@meta.data$orig.ident))

copykat.analysis <-function(rds,output,sampleid,ngene_chr,subn,species,ncore){
  sample <- unique(rds@meta.data$orig.ident)[sampleid]
  outpath <- paste0(output,"/",sample)
  if(!dir.exists(outpath)){
    dir.create(outpath)
  }
  setwd(outpath)
  
  realcell <- rownames(rds@meta.data[rds@meta.data$orig.ident==sample,])
  exp.rawdata <- as.matrix(rds@assays$RNA@counts[,colnames(rds@assays$RNA@counts) %in% realcell])
  write.table(exp.rawdata, file=paste0(sample,"_rawdata.txt"), sep="\t", quote = FALSE, row.names = TRUE)
  
  #copykat分析
  copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=ngene_chr, win.size=25, KS.cut=0.1, sam.name= sample, distance="euclidean", norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", genome=species,n.cores=ncore)
  
  
  #
  pred.test <- data.frame(copykat.test$prediction)
  rownames(pred.test) <- gsub("\\.","-",rownames(pred.test))
  pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
  CNA.test <- data.frame(copykat.test$CNAmat)
  colnames(CNA.test) <- gsub("\\.","-",colnames(CNA.test))
  
  print("Prediction result:")
  head(pred.test)
  print("Matrix:")
  head(CNA.test[ , 1:5])
  
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
  
  chr <- as.numeric(CNA.test$chrom) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)
  
  rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
  com.preN <- pred.test$copykat.pred
  pred <- rbPal5(2)[as.numeric(factor(com.preN))]
  
  cells <- rbind(pred,pred)
  col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
  
  pdf(file = paste0(sample,"_copykat.pdf"))
  heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
  
  legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")
  dev.off()
  
  #识别恶性细胞的亚克隆
  tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]
  
  tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
  hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =4, method = "euclidean"), method = "ward.D2")
  hc.umap <- cutree(hcc,subn)
  hc <- as.data.frame(hc.umap)
  write.table(hc,paste0(sample,".tumor_subtype.txt"),row.names = T,col.names = F,sep = "\t",quote = F)
  
  rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
  subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
  cells <- rbind(subpop,subpop)
  
  pdf(paste0(sample,".tumor_subtype.pdf"))
  heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
  legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')
  dev.off()
}

#分进程进行分析
system.time({
  cl = makeCluster(samplenum)
  registerDoParallel(cl)
  mydata2<-foreach(i=1:samplenum,
                   .combine = rbind,
                   .packages = c("copykat","devtools","Seurat")
  ) %dopar% copykat.analysis(rds=raw,output=output,sampleid=i,ngene_chr=ngene_chr,subn=subn,species=species,ncore=ncore)
  stopCluster(cl)
})

#copykat绘图
datasum = NULL
for(sampleid in unique(raw@meta.data$orig.ident)){
  cnvdata = read.table(paste0(output,"/",sampleid,"/",sampleid,"_copykat_prediction.txt"),sep = "\t",header = T)
  cnvdata$samples = sampleid
  datasum = rbind(datasum,cnvdata)
}
setwd(output)
print("Merged copykat result:")
dim(datasum)
#seu = subset(seu,cells = datasum$cell.names)
write.table(datasum,file = "Summary_copykat_prediction.txt",sep = "\t",quote = F,col.names = T,row.names = F)
#绘制umap图
print("Plot umap")
dat = data.frame(row.names = datasum$cell.names,copykat.pred = datasum$copykat.pred)
seu = AddMetaData(raw,dat)
p = DimPlot(seu,reduction = "umap",group.by = "copykat.pred",cols = mycol)
ggsave("copykat.umap.pdf",p,height = 10,width = 10)
ggsave("copykat.umap.png",p,height = 10,width = 10)
p = DimPlot(seu,reduction = "umap",group.by = "copykat.pred",cols = mycol,split.by = "orig.ident",ncol = ncol)
heig = ifelse(samplenum%%ncol>0,floor(samplenum/ncol)+1,samplenum/ncol)
wid = ifelse(samplenum<ncol,samplenum,ncol)
ggsave("copykat.umap.splitsample.pdf",p,height = heig*5,width = wid*5)
ggsave("copykat.umap.splitsample.png",p,height = heig*5,width = wid*5)

#绘制柱状图
print("Plot bar plot")
databar = seu@meta.data[,c("orig.ident","copykat.pred","seurat_clusters")]

stat <- databar %>%
  group_by(orig.ident) %>%
  count(copykat.pred) %>%
  rename(counts=n)
sum_counts <- stat %>%
  group_by(orig.ident) %>%
  summarise(sum(counts)) %>%
  rename(sum = 2)
stat = merge(stat,sum_counts,by="orig.ident")
stat$prop = stat$counts/stat$sum
print("Group by samples")
stat
tumors=stat[stat$copykat.pred=="aneuploid",]
stat = stat[stat$orig.ident %in% tumors$orig.ident,]
arranges = tumors[order(-tumors$counts),]

p = ggplot(stat,aes(x=factor(orig.ident,levels = arranges$orig.ident),y=counts,fill=copykat.pred))+
  geom_col(position = position_stack(reverse = TRUE))+
  labs(x="Samples",y="Cell Counts")+theme_bw()+
  scale_fill_manual(values = mycol)+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())
ggsave("copykat.counts.samplebarplot.pdf",p,height = 8,width = length(arranges$orig.ident)*0.6+2)
ggsave("copykat.counts.samplebarplot.png",p,height = 8,width = length(arranges$orig.ident)*0.6+2)

arranges = tumors[order(-tumors$prop),]
p = ggplot(stat,aes(x=factor(orig.ident,levels = arranges$orig.ident),y=prop,fill=copykat.pred))+
  geom_col(position = position_stack(reverse = TRUE))+
  labs(x="Samples",y="Cell Proportion")+theme_bw()+
  scale_fill_manual(values = mycol)+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())
ggsave("copykat.prop.samplebarplot.pdf",p,height = 8,width = length(arranges$orig.ident)*0.6+2)
ggsave("copykat.prop.samplebarplot.png",p,height = 8,width = length(arranges$orig.ident)*0.6+2)


stat <- databar %>%
  group_by(seurat_clusters) %>%
  count(copykat.pred) %>%
  rename(counts=n)
sum_counts <- stat %>%
  group_by(seurat_clusters) %>%
  summarise(sum(counts)) %>%
  rename(sum = 2)
stat = merge(stat,sum_counts,by="seurat_clusters")
stat$prop = stat$counts/stat$sum
print("Group by clusters")
stat
tumors=stat[stat$copykat.pred=="aneuploid",]
stat = stat[stat$seurat_clusters %in% tumors$seurat_clusters,]
arranges = tumors[order(-tumors$counts),]
p = ggplot(stat,aes(x=factor(seurat_clusters,levels = arranges$seurat_clusters),y=counts,fill=copykat.pred))+
  geom_col(position = position_stack(reverse = TRUE))+
  labs(x="Clusters",y="Cell Counts")+theme_bw()+
  scale_fill_manual(values = mycol)+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())
ggsave("copykat.counts.clusterbarplot.pdf",p,height = 8,width = length(arranges$seurat_clusters)*0.6+2)
ggsave("copykat.counts.clusterbarplot.png",p,height = 8,width = length(arranges$seurat_clusters)*0.6+2)

arranges = tumors[order(-tumors$prop),]
p = ggplot(stat,aes(x=factor(seurat_clusters,levels = arranges$seurat_clusters),y=prop,fill=copykat.pred))+
  geom_col(position = position_stack(reverse = TRUE))+
  labs(x="Clusters",y="Cell Proportion")+theme_bw()+
  scale_fill_manual(values = mycol)+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())
ggsave("copykat.prop.clusterbarplot.pdf",p,height = 8,width = length(arranges$seurat_clusters)*0.6+2)
ggsave("copykat.prop.clusterbarplot.png",p,height = 8,width = length(arranges$seurat_clusters)*0.6+2)

