library(Seurat)
library(ggplot2)
library(configr)
library(RcppTOML)
library(logging)
library(plyr)
library(dplyr)
library(tidyr) #### 2020-0721
library(SingleR)### 2020-0727
library(RColorBrewer)### 2020-0727
library(stringr)

options(future.globals.maxSize = 400000 * 1024^2)
basicConfig(level='FINEST')
basicConfig(level='INFO')
suppressPackageStartupMessages(library(optparse))
savepdf<- function(outpath,outf,ispdf=1,ispng=1,notgg=0,widthle=7,heightle =7){
  if(notgg==0){
    if (ispdf){
      
      ggsave(outf,file=paste0(outpath,".pdf",collapse = ""),device="pdf",width=widthle,height =heightle)
    }
    if(ispng){
      ggsave(outf,file=paste0(outpath,".png",sep=""),device="png",width=widthle,height =heightle)
    }
  }
  else{
    if (ispdf){
      
      savePlot(filename=paste0(outpath,".pdf",collapse = ""),type = "png",device=dev.cur())
    }
    if(ispng){
      #ggsave(outf,file=paste0(outpath,".png",sep=""),device="png",width=5,height=5)
      loginfo("pass")
    }
    
  }
  loginfo("savefig at %s",outpath)
}

isstep <-function(x){
  outb=FALSE  
  c=paste0("outb=csc$run$step",x)
  tryCatch(
    {
      eval(parse(text=c))
      if(is.null(outb)){
        outb=FALSE
      }
      
    },error = function ( e ) {
      outb=FALSE
      logerror(print(outb))
      return(outb)
      
    },finally = {
      return(outb)
    }
  )
}

getrds <-function(x){
  backobsig=1
  loginfo("read rds")
  c=paste0(outtemp,"step",x,".RDS")
  loginfo(c)
  if(file.exists(c)){
    loginfo("import %s",c)
    backob=readRDS(c)
    backobsig=0
  }
  else if(!is.null( csc$tempdata$tempdata)){
    if(file.exists(csc$tempdata$tempdata)){
      loginfo("import otherdata %s",csc$tempdata$tempdata)
      backob=readRDS(csc$tempdata$tempdata)
      backobsig=0
    }
    
    
  }
  if(backobsig==1){
    logerror("can't find file exist in %s or %s",c,csc$tempdata$tempdata)
    stop()
  }
  return(backob)
}
syscomd<-function(commada,savef="",runl=T){
  if(savef!=""){
    cat(paste0(commada,"\n"),file=savef,append = T)
  }
  if(runl==T){
    loginfo(commanda)
    system(commanda)
  }
}


option_list=list(make_option(c("-i","--infile"),  action="store", help='The input file.'),
                 make_option(c("-p","--path"),  action="store", help='The input file.',default = ""))
opt<-parse_args(OptionParser(usage="%prog [options] file\n",option_list=option_list))
#opt$infile="/haplox/haprs/liyq/singlecell/test/mmu10/test/cofige.hapyun.ini"
#opt$infile="/data/haplox/users/liyq/singlecell/deadline/confige.ini"   shielded by zhangdx
csc=parseTOML(opt$infile)
outdir=csc$outpath$outpath

wdpath=getwd()
if(csc$outpath$outpath %>% str_detect("^/")){
  loginfo("outfile in %s",csc$outpath$outpath)
  outdir=csc$outpath$outpath}else{
    outdir=paste(wdpath,outdir,sep = "/")
    loginfo("outfile in %s",outdir)
    
  }

numsap=length(csc$indata)
dir.create(outdir)
dir.create(paste(outdir,"temp",sep="/"))
outtemp=paste(outdir,"temp","",sep="/")
if(isstep(1)){
  loginfo("step1 start")
  dir.create(paste(outdir,"step1",sep="/"),recursive = T)
  listobjall=list()
  cat("name\tbefore\tafter\tfiltered\n",file = paste(outdir,"step1","summary.xls",sep="/"),fill= F, labels=NULL, append=F)
  for (i in 1:length(csc$indata))
  {
    csc$step1$nFeature_RNA=as.numeric(csc$step1$nFeature_RNA)
    csc$step1$percent_mt=as.numeric(csc$step1$percent_mt)
    
    dir.create(paste(outdir,"step1",names(csc$indata[i]),sep="/"))
    if(csc$step1$filetype=="10x"){
      if(csc$indata[[i]] %>% str_detect("^/")){
        loginfo("indata from %s",csc$indata[[i]])
        
        pbmc_data= Read10X(data.dir = csc$indata[[i]])
        
        
      }else{
        csc$indata[[i]]=paste(opt$path,csc$indata[[i]],sep = "/")
        loginfo("indata from %s",csc$indata[[i]])
        
        pbmc_data= Read10X(data.dir = csc$indata[[i]])
        
      }
      
    }else if(csc$step1$filetype=="csv"){
      pbmc_data=read.table(csc$indata[[i]],sep = csc$step1$csv_sep,header = T,row.names = 1)
    }else{logerror("this is cant open ,please cheack")
      exists()}
    pbmc=CreateSeuratObject(counts = pbmc_data,project = names(csc$indata[i])) #project in orig.ident
    #pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = paste0("^",csc$step1$mttype,"-"))
    col2=length(pbmc@active.ident)
    
    mt_vln=VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=1) 
    savepdf(paste(outdir,"step1",names(csc$indata[i]),paste0(names(csc$indata[i]),"_voilin"),sep="/"),mt_vln)
    featureqc1=FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")+ geom_hline(aes(yintercept=csc$step1$percent_mt[2]),colour="#BB0000", linetype="dashed") + ggtitle("QCmt")
    featureqc2=FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ geom_hline(aes(yintercept=csc$step1$nFeature_RNA[1]),colour="#BB0000", linetype="dashed") + geom_hline(aes(yintercept=csc$step1$nFeature_RNA[2]),colour="#BB0000", linetype="dashed") + ggtitle("QCng")
    
    df <- data.frame(nc = pbmc@meta.data$nCount_RNA,pm = pbmc@meta.data$percent.mt)
    p <- ggplot(df,aes(x=nc,y=pm,colour=nc))+geom_point(size=1,shape=16)+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA),panel.border=element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+scale_y_continuous(expand=c(0,0))+theme(panel.border = element_blank(), axis.line = element_line(size = 0.7))+scale_color_gradientn(colours = c("#0066CC","#FFCC33","#FF0000"))+labs(x="nCount_RNA",y="percent.mt",colour="nCount_RNA")
    p <- p+theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),axis.title.x = element_text(size=14),axis.title.y = element_text(size = 14))+geom_hline(aes(yintercept=csc$step1$percent_mt[2]),colour="#BB0000", linetype="dashed")
    
    savepdf(paste(outdir,"step1",names(csc$indata[i]),paste0(names(csc$indata[i]),"_countVmt"),sep="/"),featureqc1)
    savepdf(paste(outdir,"step1",names(csc$indata[i]),paste0(names(csc$indata[i]),"_countVfeature"),sep="/"),featureqc2)
    savepdf(paste(outdir,"step1",names(csc$indata[i]),paste0(names(csc$indata[i]),"_libraryVmt"),sep="/"),p)    
    
    pbmc<- subset(pbmc,subset = nFeature_RNA > csc$step1$nFeature_RNA[1]& nFeature_RNA < csc$step1$nFeature_RNA[2] & percent.mt >=  csc$step1$percent_mt[1] & percent.mt < csc$step1$percent_mt[2])
    
    col3=length(pbmc@active.ident)
    
    listobjall[[i]]=pbmc
    
    cat(pbmc@project.name,
        col2,col3,col2-col3,"\n",
        sep = "\t",
        file = paste(outdir,"step1","summary.xls",sep="/"),fill= F, labels=NULL, append=T)
  }
  loginfo("save %s",paste0(outtemp,"step1.RDS"))
  saveRDS(listobjall,paste0(outtemp,"step1.RDS"))  
  loginfo("step1 end")
}
if (isstep(2)){
  loginfo("step2 start")
  if (!isstep(1)){
    loginfo("input step1 rds")
    listobjall=getrds(1)
  }
  if(length(listobjall)>1){
    
    
    if(csc$step2$normethod=="SCT"){
      filternum=NULL
      for (i in 1:length(listobjall)) {
        listobjall[[i]] <- SCTransform(listobjall[[i]])
        if(csc$step2$kfilter>length(listobjall[[i]]@active.ident)){
          logwarn("%s numcell %s less than %s,plesae check it" ,listobjall[[i]]@project.name,length(listobjall[[i]]@active.ident),csc$step2$kfilter)
          
          next
        }
        filternum=c(filternum,i)
        
      }
      listobj=listobjall[filternum]
      #listname<- unlist(lapply(listobj,levels))
      listname<- unlist(lapply(listobj,function(x){x@project.name}))
      usereflist<-unlist(lapply(csc$step2$reference,function(x){which(listname==x)}))
      if(length(usereflist)==0){
        usereflist=NULL
      }else{    logwarn("userlist in samplelist ")
      }
      pbmc_features<-SelectIntegrationFeatures(object.list = listobj, nfeatures = csc$step2$nFeature)
      listobj <- PrepSCTIntegration(object.list = listobj, anchor.features = pbmc_features)
      pbmc_anchors <- FindIntegrationAnchors(object.list = listobj, normalization.method = "SCT", k.filter = csc$step2$kfilter,anchor.features = pbmc_features,reference = usereflist)
      pbmcall <- IntegrateData(anchorset = pbmc_anchors, normalization.method = "SCT")
    }else if(csc$step2$normethod=="vst"){
      filternum=NULL
      for (i in 1:length(listobjall)) {
        logdebug("ok")
        if(csc$step2$kfilter>length(listobjall[[i]]@active.ident)){
          logwarn("%s numcell %s less than %s,plesae check it" ,listobjall[[i]]@project.name,length(listobjall[[i]]@active.ident),csc$step2$kfilter)
          
        }else{filternum=c(filternum,i)}
        listobjall[[i]] <- NormalizeData(listobjall[[i]])
        listobjall[[i]] <- FindVariableFeatures(listobjall[[i]], selection.method = "vst",nfeatures = csc$step2$nFeature)
      }
      listobj=listobjall[filternum]
      pbmc_anchors <- FindIntegrationAnchors(object.list = listobj,normalization.method = "LogNormalize",k.filter = csc$step2$kfilter)
      pbmcall <- IntegrateData(anchorset = pbmc_anchors, dims = 1:30)
      pbmcall <- ScaleData(pbmcall, verbose = FALSE)
    }else if(csc$step2$normethod=="none"){
      mergepb=merge(listobjall[[1]],listobjall[-1])
      mergepb=NormalizeData(mergepb,normalization.method = "LogNormalize")
      mergepb=FindVariableFeatures(mergepb,selection.method = "vst",nfeatures = csc$step2$nFeature)
      pbmcall <- ScaleData(mergepb, verbose = FALSE)
    }
    
  }else{
    mergepb= NormalizeData(listobjall[[1]],normalization.method = "LogNormalize")
    mergepb=FindVariableFeatures(mergepb,selection.method = "vst",nfeatures = csc$step2$nFeature)
    pbmcall <- ScaleData(mergepb, verbose = FALSE)
    
  }
  loginfo("save %s",paste0(outtemp,"step2.RDS"))
  saveRDS(pbmcall,paste0(outtemp,"step2.RDS"))
  loginfo("step2 end")
}

if (isstep(3)){
  loginfo("step3 start")
  dir.create(paste(outdir,"step3",sep="/"))
  if (!isstep(2)){
    loginfo("input step2 rds")
    pbmcall=getrds(2)
  }
  pbmcall=RunPCA(pbmcall)
  dir.create(paste(outdir,"step3","pca",sep="/"))
  fig1=DimPlot(object = pbmcall, reduction = "pca",group.by = "orig.ident")+
    theme(legend.text = element_text(size =8))+
    guides(color=guide_legend(ncol = 1,override.aes = list(size=3)))
  savepdf(paste(outdir,"step3","pca","pca",sep="/"),fig1,widthle = 10,heightle = 8)
  for(i in 1:csc$step3$heatmapnumber){
    pdf(paste(outdir,"step3","pca",paste0("pcaheatmap",i,".pdf"),sep="/"))
    DimHeatmap(object = pbmcall ,cells=1000, balanced = TRUE ,reduction = "pca",dims=i)
    dev.off()
    png(paste(outdir,"step3","pca",paste0("pcaheatmap",i,".png"),sep="/"))
    DimHeatmap(object = pbmcall ,cells=1000, balanced = TRUE ,reduction = "pca",dims=i)
    dev.off()
  }
  bow_plot=ElbowPlot(pbmcall,csc$step3$elbowdims)#    add parameters: csc$step3$elbowdims by zhangdx
  savepdf(paste(outdir,"step3","pca","bowplot",sep="/"),bow_plot)
  
  cscumap=csc$step3$reduction
  
  if(cscumap=="umap"){
    #    pbmcall=RunUMAP(pbmcall,dims = 1:csc$step3$dims,umap.method="umap-learn",metric="correlation")
    pbmcall=RunUMAP(pbmcall,dims = 1:csc$step3$dims,seed.use = 100)
    
    
  }else{
    pbmcall=RunTSNE(pbmcall,dims=1:csc$step3$dims,seed.use = 100)
  }
  dir.create(paste(outdir,"step3",cscumap,sep="/"))
  umaporig=DimPlot(pbmcall,reduction = cscumap,split.by = "orig.ident",group.by = "orig.ident",ncol= 4)+
    theme(legend.text = element_text(size =8))+
    guides(color=guide_legend(ncol = 1,override.aes = list(size=3)))
  
  savepdf(paste(outdir,"step3",cscumap,"plotby_ident",sep="/"),umaporig,widthle = min(numsap*4,16),heightle = ceiling(numsap/4)*4)
  umaporig2=DimPlot(pbmcall,reduction = cscumap,group.by = "orig.ident")+
    theme(legend.text = element_text(size =8))+
    guides(color=guide_legend(ncol = 1,override.aes = list(size=3)))
  
  savepdf(paste(outdir,"step3",cscumap,"plotall_ident",sep="/"),umaporig2,widthle = 10,heightle = 8)
  
  pbmcall <- AddMetaData(pbmcall,pbmcall@reductions$umap@cell.embeddings,col.name = colnames(pbmcall@reductions$umap@cell.embeddings))
  umap <- ggplot(pbmcall@meta.data,aes(x=UMAP_1,y=UMAP_2))+geom_point(size =1,aes(color=log2(nCount_RNA)))+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA),panel.border=element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+theme(panel.border = element_blank(), axis.line = element_line(size = 0.7))+scale_color_gradientn(colours = c("#0066CC","#FFCC33","#FF0000"))
  savepdf(paste(outdir,"step3",cscumap,"plotby_nCount",sep="/"),umap,widthle = 10,heightle = 8)
  
  if(csc$step3$clustercell){
    pbmcall <- FindNeighbors(pbmcall, dims = 1:csc$step3$dims)
    pbmcall <- FindClusters(pbmcall, resolution = csc$step3$resolution,algorithm = csc$step3$algorithm)
    umaporig3=DimPlot(pbmcall,reduction = cscumap,group.by = "seurat_clusters")+
      theme(legend.text = element_text(size =12))+
      guides(color=guide_legend(ncol = 1,override.aes = list(size=3)))
    
    savepdf(paste(outdir,"step3",cscumap,"plotby_cluster",sep="/"),umaporig3,widthle = 10,heightle = 8)
    #pbmcall.markers <- FindAllMarkers(pbmcall,only.pos = TRUE,)
    ##########################################################
    test=pbmcall@meta.data %>% group_by(seurat_clusters,orig.ident) %>% summarise(count=n()) %>% spread(key = seurat_clusters,value = count)
    test2=pbmcall@meta.data %>% group_by(orig.ident) %>% summarise(count=n())
    test3=merge(test2,test)
    test3[is.na(test3)]<-0
    
    datasum <- apply(test3[,-1],2,sum)
    sumcol <- c("all",datasum)
    finaldata <- rbind(test3,sumcol)
    
    write.table(finaldata,file = paste(outdir,"step3",cscumap,"summary.xls",sep="/"),sep = "\t",quote=FALSE,row.names = FALSE,col.names = TRUE)
    
    total3=test3 %>% gather(cluster,counts,-orig.ident)
    barpli = ggplot(total3,aes(cluster,counts,fill=orig.ident)) + geom_bar(stat = 'identity',position = 'fill') + scale_x_discrete(limits = colnames(test3)[-1])
    
    savepdf(paste(outdir,"step3",cscumap,"barplot",sep="/"),barpli,widthle = 10, heightle = 8) 
    ###############################################################
    hpca.se <- readRDS(csc$step3$singler)
    #largdata=readRDS("/data/haplox/users/liyq/singlecell/bin/singleRdata/DatabaseImmuneCellExpressionDataRDS.rds")
    
    testcd=GetAssayData(pbmcall)
    immannot=SingleR(test = testcd,ref = hpca.se,labels = hpca.se$label.fine)
    foo=immannot$labels
    names(foo)=immannot@rownames
    pbmcall$singleR=foo
    
    singRtxt= pbmcall@meta.data %>% group_by(seurat_clusters,singleR) %>% summarise(n_cells=n()) %>% arrange(seurat_clusters,desc(n_cells))
    x <- vector()
    for (i in levels(pbmcall$seurat_clusters)) {
      spart = subset(singRtxt,singRtxt$seurat_clusters == i)
      spart = subset(spart,(spart$n_cells/sum(spart$n_cells))>0.02)
      x <- union(x,spart$singleR)
    }  
    x <- sort(x)
    pbmcall1 = pbmcall
    for (j in 1:length(pbmcall$singleR)) {
      if (!pbmcall$singleR[j] %in% x){
        pbmcall1$singleR[j] = 'others'
      }
    }
    
    cols <- brewer.pal(7,"Spectral")
    pal <- colorRampPalette(cols[c(7,6,2,1)])
    mycolors <- pal(length(x)+1)
    
    anovplot = DimPlot(pbmcall1,reduction = cscumap,label = F,cols= mycolors,group.by = "singleR")
    levels(anovplot$data$singleR) <- c(x,'others')
    
    anovplot <- anovplot + theme(legend.text = element_text(size =5))+ guides(color=guide_legend(ncol = 1,keyheight = 0.5,override.aes = list(size=1.5))) 
    savepdf(paste(outdir,"step3",cscumap,"celltype",sep="/"),anovplot,widthle = 10, heightle = 8)
    
    ck <- c(x,'others')
    dzb <- cbind(ck,mycolors)
    
    write.table(singRtxt,file=paste(outdir,"step3",cscumap,"singleR_celltype.xls",sep="/"),sep='\t',row.names = F, col.names = TRUE,quote = F)
    
    s1 = pbmcall1@meta.data %>% group_by(seurat_clusters,singleR) %>% summarise(n_cells=n()) %>% arrange(seurat_clusters,desc(n_cells))
    for (i in levels(pbmcall$seurat_clusters)) {
      
      s10 = subset(s1,s1$seurat_clusters== i)
      
      s11 = s10 %>% arrange(singleR)
      s11$singleR=factor(s11$singleR,levels=levels(anovplot$data$singleR))
      s11 = s11 %>% arrange(factor(s11$singleR))
      
      thiscolor <- vector()
      for (j in 1:length(s11$singleR)) {
        for (z in 1:length(dzb[,1])) {
          if (s11$singleR[j]==dzb[z,1]) {
            thiscolor[j] = dzb[z,2]
          }
        }
      } 
      
      p <- ggplot(data = s11, mapping = aes(x='Content',y = n_cells, fill = singleR)) + geom_bar(stat = 'identity', position = 'stack',width=1,aes(fill = singleR))
      aaa <- round(s11$n_cells/sum(s11$n_cells)*100,2)
      label_value <- ifelse(aaa > 2, paste('(',aaa,'%)',sep=''),"")
      p <- p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') +theme(panel.grid.major=element_line(colour=NA),panel.border=element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) +labs(fill= paste("Cell type in cluster ",i,sep = ''))
      p <- p+geom_text(stat="identity",aes(x =sum(s11$n_cells)/(sum(s11$n_cells/1.8)),y=s11$n_cells, label = label_value), size=3, position=position_stack(vjust = 0.5))+scale_fill_manual(values = thiscolor)+
        guides(fill = guide_legend( ncol = length(aaa)%/%20+1, byrow = TRUE))
      
      savepdf(paste(outdir,"step3",cscumap,paste("cluster",i,sep = ''),sep="/"),p,widthle = ifelse(length(aaa)>30,length(aaa)%/%20*5+5,10), heightle = 8)
    }  
    
    
  }
  loginfo("save %s",paste0(outtemp,"step3.RDS"))
  saveRDS(pbmcall,paste0(outtemp,"step3.RDS"))
  loginfo("step3 end")
  
}

if (isstep(4)){
  loginfo("step4 start")
  dir.create(paste(outdir,"step4",sep="/"))
  if (!isstep(3)){
    loginfo("input step3 rds")
    pbmcall=getrds(3)
  }
  if(csc$step4$clustermarkers){
    dir.create(paste(outdir,"step4",sep="/"))
    pbmcall.markers <- FindAllMarkers(pbmcall,only.pos = TRUE,test.use = csc$step4$findmarkers_testuse,min.pct = csc$step4$min_pct)
    top10<- subset(pbmcall.markers,avg_log2FC > 0.5 & p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(10, avg_log2FC)
    sorttop10 <- top10[order(top10$cluster,-top10$avg_log2FC),]
    pbmcall.markerss=pbmcall.markers[,c(ncol(pbmcall.markers),2:ncol(pbmcall.markers)-1)]
    colnames(pbmcall.markerss)[1]="Gene"
    write.table(pbmcall.markerss,file=paste(outdir,"step4",paste0("clusterall_genelist.xls"),sep="/"),sep='\t',row.names = F,quote = F)
    write.table(sorttop10,file=paste(outdir,"step4",paste0("clusterall_top10genelist.xls"),sep="/"),sep='\t',row.names = TRUE, col.names = NA,quote = F)
    write.table(subset(pbmcall.markers,avg_log2FC > 0.5 & p_val_adj < 0.05),file=paste(outdir,"step4",paste0("clusterall_adj0.05_logFC0.5genelist.xls"),sep="/"),sep='\t',row.names = TRUE, col.names = NA,quote = F)
    #top10=read.csv("/haplox/haprs/liyq/Project/HGC20220110002_sc_guangzhoudaxue/workfile/step4/clusterall_top10genelist.xls",sep = "\t",quote = "",header = 1)
    #pdf(paste(outdir,"step4",paste0("cluster_top10geneheatmap.pdf"),sep="/"))
    #DoHeatmap(pbmcall,features = top10$gene) + scale_y_discrete(breaks=NULL)
    #dev.off()
    subobj=subset(pbmcall,downsample=1000)
    plot <- DoHeatmap(object = subobj,features = top10$gene) +NoLegend() + scale_y_discrete(breaks=NULL)
    savepdf(paste(outdir,"step4",paste0("cluster_top10geneheatmap"),sep="/"),plot)
    #png(paste(outdir,"step4",paste0("cluster_top10geneheatmap.png"),sep="/"))
    #DoHeatmap(pbmcall,features = top10$gene) + scale_y_discrete(breaks=NULL)
    #dev.off()
    
    pbmcall@active.assay <- "RNA"
    for (i in levels(pbmcall$seurat_clusters)){
      logdebug(i)
      dir.create(paste(outdir,"step4",paste0("cluster",i),sep="/"))
      write.table(subset(pbmcall.markers,cluster==i),file=paste(outdir,"step4",paste0('cluster',i,"/cluster",i,"_genelist.xls"),sep="/"),sep='\t',row.names = TRUE, col.names = NA,quote = F)
      write.table(subset(pbmcall.markers,cluster==i & avg_log2FC > 0.5 & p_val_adj < 0.05),file=paste(outdir,"step4",paste0('cluster',i,"/cluster",i,"_padj0.05_logFC0.5genelist.xls"),sep="/"),sep='\t',row.names = TRUE, col.names = NA,quote = F)
      if(length(subset(top10,cluster==i)$gene)!=0){
        logdebug("abc")
        
        fig4.1=VlnPlot(pbmcall,features=c(subset(top10,cluster==i)$gene),ncol=4,pt.size=-1)
        savepdf(paste(outdir,"step4",paste0('cluster',i,"/cluster",i,"_genevlnplot"),sep="/"),fig4.1,widthle = 20,heightle = 15)###add l
        fig4.2=FeaturePlot( pbmcall, 
                            features = c(subset(sorttop10,cluster==i)$gene), 
                            cols = c("lightgrey", "blue"),
                            ncol = 4,reduction = csc$step3$reduction)
        savepdf(paste(outdir,"step4",paste0('cluster',i,"/cluster",i,"_genereduction"),sep="/"),fig4.2,widthle = 20,heightle = 15)
      }
    } 
    if(!is.null(csc$step4$custer)){
      logdebug("ok")
      dir.create(paste(outdir,"step4","custer",sep="/"))
      for(i in names(csc$step4$custer)){
        
        #igene=eval(parse(text=paste0("csc$step4$custer$",i)))
        igene=csc$step4$custer[[i]] #20210926
        dir.create(paste(outdir,"step4","custer",i,sep="/"))
        noexistgene=setdiff(igene,row.names(pbmcall@assays$RNA))
        igene=intersect(igene,row.names(pbmcall@assays$RNA))
        write.table(igene,file = paste(outdir,"step4","custer",i,"use_genelist",sep="/"),sep='\n',row.names = F, col.names = F,quote = F)
        if(length(igene)==0){
          logwarn("no one gene in it please check it: %s" , noexistgene)
          
          next
        }
        logwarn("list of not exist gene: %s" , noexistgene)
        for (m in igene){
          logdebug(m)
          fig4.3=VlnPlot(pbmcall,features = m,pt.size=-1)
          savepdf(paste(outdir,"step4","custer",i,paste0(m,"_vlnplot"),sep="/"),fig4.3)
          fig4.4=FeaturePlot( pbmcall, 
                              features = m, 
                              cols = c("lightgrey", "blue"),
                              reduction = csc$step3$reduction)
          savepdf(paste(outdir,"step4","custer",i,paste0(m,"_reduction"),sep="/"),fig4.4)
        }
        if(length(names(pbmcall@assays))>1){
          pbmcall@active.assay <- names(pbmcall@assays)[2]}
        fig4.5=DotPlot(pbmcall, features = rev(igene), dot.scale = 8) +
          scale_color_gradient(low = "#132B43", high = "#56B1F7") +
          theme(axis.text.x  = element_text(angle = 90))
        savepdf(paste(outdir,"step4","custer",i,"use_genedotplot",sep="/"),fig4.5,widthle = 3+1.5*length(igene),heightle = 4+0.35*length(levels(pbmcall$seurat_clusters)))
      }
      
    }
    
  }
  if(!is.null(csc$step4$difcluster)){
    logdebug("difcluster为寻找不同分组的cluster之间的差异基因")
    dir.create(paste(outdir,"step4","difcluster",sep="/"))
    for (i in names(csc$step4$difcluster)){
      logdebug(i)
      # idif=eval(parse(text=paste0("csc$step4$difcluster$",i)))
      idif=csc$step4$difcluster[[i]] #20210926
      
      dir.create(paste(outdir,"step4","difcluster",i,sep="/"))
      difgenelist=FindMarkers(pbmcall,ident.1 = idif$a,ident.2 = idif$b,group.by = "seurat_clusters",test.use = idif$testuse)
      difgenelist$gene <- rownames(difgenelist)
      dif_filgenelist=subset(difgenelist, avg_log2FC > 0.5 & p_val_adj < 0.05)
      write.table(difgenelist,file=paste(outdir,"step4","difcluster",i,"dif_allgenelist.xls",sep="/"),sep='\t', row.names = F,col.names = TRUE,quote = F)
      write.table(dif_filgenelist,file=paste(outdir,"step4","difcluster",i,"dif_logFC0.5padj0.05genelist.xls",sep="/"),sep='\t',row.names = F, col.names = TRUE,quote = F)
      #展示差异的前十个基因
      arrange(dif_filgenelist,desc(avg_log2FC))
      fig4.2=FeaturePlot( pbmcall, 
                          features = head(arrange(dif_filgenelist,desc(avg_log2FC)),n = 10)$gene, 
                          cols = c("lightgrey", "blue"),
                          ncol = 4,reduction = csc$step3$reduction)
      savepdf(paste(outdir,"step4","difcluster",i,"dif_genereduction",sep="/"),fig4.2,widthle = 20,heightle = 15)
      if (!is.null(csc$step4$ClusterProfiler) & csc$step4$ClusterProfiler[1]=="true" ){
        dir.create(paste(outdir,"step4","difcluster",i,"clusterprofiler",sep="/"))
        setwd(paste(outdir,"step4","difcluster",i,"clusterprofiler",sep="/"))
        cat(csc$step4$ClusterProfiler[-1],
            "-f",paste(outdir,"step4","difcluster",i,"dif_logFC0.5padj0.05genelist.xls",sep="/"),
            "-n",paste0("\"",csc$title,"\""),
            "-o",paste(outdir,"step4","difcluster",i,"clusterprofiler",sep="/"),
            file = paste(outdir,"step4","difcluster",i,"clusterprofiler","run.sh",sep="/"),
            sep=" ", fill= F, labels=NULL, append=F)
        system(paste0("bash ",paste(outdir,"step4","difcluster",i,"clusterprofiler","run.sh",sep="/")))
      }
    }
    
  }
  if(!is.null(csc$step4$difident)){
    logdebug("difident为寻找样品间差异基因")
    dir.create(paste(outdir,"step4","difident",sep="/"))
    for (i in names(csc$step4$difident)){
      logdebug(i)
      # idif=eval(parse(text=paste0("csc$step4$difident$",i)))
      
      idif=csc$step4$difident[[i]] #20210926
      
      dir.create(paste(outdir,"step4","difident",i,sep="/"))
      
      
      difgenelist=FindMarkers(pbmcall,ident.1 = idif$a,ident.2 = idif$b,group.by = "orig.ident",test.use = idif$testuse)
      difgenelist$gene <- rownames(difgenelist)
      dif_filgenelist=subset(difgenelist, avg_log2FC > 0.5 & p_val_adj < 0.05)
      write.table(difgenelist,file=paste(outdir,"step4","difident",i,"dif_allgenelist.xls",sep="/"),sep='\t', row.names = F,col.names = TRUE,quote = F)
      write.table(dif_filgenelist,file=paste(outdir,"step4","difident",i,"dif_logFC0.5padj0.05genelist.xls",sep="/"),sep='\t',row.names = F, col.names = TRUE,quote = F)
      #展示差异的前十个基因
      arrange(dif_filgenelist,desc(avg_log2FC))
      fig4.2=FeaturePlot( pbmcall, 
                          features = head(arrange(dif_filgenelist,desc(avg_log2FC)),n = 10)$gene, 
                          cols = c("lightgrey", "blue"),
                          ncol = 4,reduction = csc$step3$reduction)
      savepdf(paste(outdir,"step4","difident",i,"dif_genereduction",sep="/"),fig4.2,widthle = 20,heightle = 15)
      if (!is.null(csc$step4$ClusterProfiler) & csc$step4$ClusterProfiler[1]=="true" ){
        dir.create(paste(outdir,"step4","difident",i,"clusterprofiler",sep="/"))
        setwd(paste(outdir,"step4","difident",i,"clusterprofiler",sep="/"))
        cat(csc$step4$ClusterProfiler[-1],
            "-f",paste(outdir,"step4","difident",i,"dif_logFC0.5padj0.05genelist.xls",sep="/"),
            "-n",paste0("\"",csc$title,"\""),
            "-o",paste(outdir,"step4","difident",i,"clusterprofiler",sep="/"),
            file = paste(outdir,"step4","difident",i,"clusterprofiler","run.sh",sep="/"),
            sep=" ", fill= F, labels=NULL, append=F)
        system(paste0("bash ",paste(outdir,"step4","difident",i,"clusterprofiler","run.sh",sep="/")))
      }
    }
    
    
  }
  loginfo("save %s",paste0(outtemp,"step4.RDS"))
  saveRDS(pbmcall,paste0(outtemp,"step4.RDS"))
  loginfo("step4 end")
  
}

#pbmcall=readRDS("/haplox/haprs/liyq/singlecell/test/mmu10/test/test/workout/temp/step4.RDS")
if (isstep(5)){
  loginfo("step5 start")
  dir.create(paste(outdir,"step5",sep="/"))
  if (!isstep(4)){
    loginfo("input step3 rds")
    pbmcall=getrds(4)
  }
  library(monocle)
  #从Seurat转换
  if(ncol(pbmcall@assays$RNA@data)>10000){
    mincells=sample(seq_len(ncol(pbmcall@assays$RNA@data)),size=10000)
    datamon=as(pbmcall@assays$RNA@data[,mincells,drop=FALSE],"sparseMatrix")
    
    pd<-new("AnnotatedDataFrame", data = pbmcall@meta.data[colnames(datamon),])
    
  }else{
    datamon=as(pbmcall@assays$RNA@data,'sparseMatrix')
    pd<-new("AnnotatedDataFrame", data = pbmcall@meta.data)
    
  }
  fd<-new("AnnotatedDataFrame", data = data.frame(gene_short_name = row.names(datamon), row.names = row.names(datamon)))
  cds <- newCellDataSet(datamon, phenoData = pd, featureData = fd)
  sce = cds
  
  sce <- estimateSizeFactors(sce)
  sce <- estimateDispersions(sce)
  
  #过滤细胞QC
  sce <- detectGenes(sce, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(sce),num_cells_expressed >= 10))
  pData(sce)$Total_mRNAs <- Matrix::colSums(Biobase::exprs(sce))
  sce <- sce[, pData(sce)$Total_mRNAs < 1e6 ]
  upper_bound <- 10^(mean(log10(pData(sce)$Total_mRNAs)) +
                       2*sd(log10(pData(sce)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(sce)$Total_mRNAs)) -
                       2*sd(log10(pData(sce)$Total_mRNAs)))
  fig5.1=qplot(Total_mRNAs, data = pData(sce), color = orig.ident, geom = "density") +
    geom_vline(xintercept = lower_bound) +
    geom_vline(xintercept = upper_bound) 
  savepdf(paste(outdir,"step5","QC",sep="/"),fig5.1)
  
  sce <- sce[,pData(sce)$Total_mRNAs > lower_bound &
               pData(sce)$Total_mRNAs < upper_bound]
  sce <- detectGenes(sce, min_expr = 0.1)
  
  
  #选择合适的基因来标记状态
  #diff_test_res <- differentialGeneTest(sce[expressed_genes,],
  #                                     fullModelFormulaStr = "~seurat_clusters")
  #sce_ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
  disp_table <- dispersionTable(sce)
  sce_ordering_genes <- subset(disp_table,mean_expression >= csc$step5$meanexpression &dispersion_empirical >= 1 * dispersion_fit)$gene_id
  sce <-
    setOrderingFilter(sce,
                      ordering_genes = sce_ordering_genes)
  fig5.2=plot_ordering_genes(sce)
  dir.create(paste(outdir,"step5","pseudotime",sep="/"))
  savepdf(paste(outdir,"step5","pseudotime","geneplot",sep="/"),fig5.2)
  
  
  sce <- reduceDimension(sce, method = 'DDRTree')
  sce <- orderCells(sce)
  fig5.3=plot_cell_trajectory(sce, color_by = "seurat_clusters")
  savepdf(paste(outdir,"step5","pseudotime","pseudotime_byclusters",sep="/"),fig5.3)
  fig5.4=plot_cell_trajectory(sce, color_by = "Pseudotime")
  savepdf(paste(outdir,"step5","pseudotime","pseudotime_byPseudotime",sep="/"),fig5.4)
  fig5.5=plot_cell_trajectory(sce, color_by = "seurat_clusters") +     facet_wrap(~seurat_clusters, nrow = ceiling(length(levels(sce$seurat_clusters))/4))
  savepdf(paste(outdir,"step5","pseudotime","pseudotime_splitbyclusters",sep="/"),fig5.5,widthle = 12,heightle = 4+3*ceiling(length(levels(sce$seurat_clusters))/4))
  fig5.6=plot_cell_trajectory(sce, color_by = "orig.ident")+     facet_wrap(~orig.ident, nrow = ceiling(length(union(sce$orig.ident,NULL))/4))
  savepdf(paste(outdir,"step5","pseudotime","pseudotime_splitbyorig.ident",sep="/"),fig5.6,widthle = 12,heightle = 4+3*ceiling(length(union(sce$orig.ident,NULL))/4))
  fig5.7=plot_cell_trajectory(sce, color_by = "State")
  savepdf(paste(outdir,"step5","pseudotime","pseudotime_bystate",sep="/"),fig5.7)
  # write.table(select(pData(sce),State,orig.ident,seurat_clusters),file=paste(outdir,"step5","pseudotime","cell_state_clusters.csv",sep="/"),sep='\t',row.names = T, col.names = NA,quote = F)
  write.table(pData(sce)[,c("State","orig.ident","seurat_clusters")],file=paste(outdir,"step5","pseudotime","cell_state_clusters.csv",sep="/"),sep='\t',row.names = T, col.names = NA,quote = F)
  
  #差异表达分析
  diff_test_res <- differentialGeneTest(sce[expressed_genes, ],
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  sig_genes <- subset(diff_test_res, qval < 0.01)
  sig_genes <- sig_genes[order(sig_genes$qval), ]
  dir.create(paste(outdir,"step5","diffgenes",sep="/"))
  #  write.table(select(sig_genes,pval,qval,num_cells_expressed),file=paste(outdir,"step5","diffgenes","cell_diffgene.csv",sep="/"),sep='\t',row.names = T, col.names = NA,quote = F)
  write.table(sig_genes[,c("pval","qval","num_cells_expressed")],file=paste(outdir,"step5","diffgenes","cell_diffgene.csv",sep="/"),sep='\t',row.names = T, col.names = NA,quote = F)
  
  pdf(paste(outdir,"step5","diffgenes","topgeneheatmap.pdf",sep="/"))
  fig5.8=plot_pseudotime_heatmap(sce[row.names(sig_genes)[1:csc$step5$genenum],],
                                 num_clusters = csc$step5$numclusters,
                                 cores = 1,
                                 show_rownames = T)
  
  dev.off()
  #BEAM分析
  BEAM_res <- BEAM(sce[expressed_genes, ], branch_point = csc$step5$pointid, cores = 16)
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  dir.create(paste(outdir,"step5","BEAM",sep="/"))
  write.table(BEAM_res,file=paste(outdir,"step5","BEAM","diffgenelist.xls",sep="/"),sep='\t',row.names = F, col.names = T,quote = F)
  pdf(paste(outdir,"step5","BEAM","topgeneheatmap.pdf",sep="/"))
  
  plot_genes_branched_heatmap(sce[row.names(subset(BEAM_res,
                                                   qval < 1e-4))[1:csc$step5$BEAMgn],],
                              branch_point = csc$step5$pointid,
                              num_clusters = csc$step5$BEAMnumclusters,
                              cores = 1,
                              use_gene_short_name = T,
                              show_rownames = T)
  
  
  dev.off()
  fig5.9=plot_genes_branched_pseudotime(sce[rownames(subset(BEAM_res, qval < 1e-8))[1:9], ], 
                                        branch_point = csc$step5$pointid, color_by = "State", 
                                        ncol = 3)
  savepdf(paste(outdir,"step5","BEAM","top10gene",sep="/"),fig5.9)
  
  length(intersect( rownames(sce),csc$step5$BEAMgenelist))
  
  if(!is.null(csc$step5$BEAMgenelist)&&  length(intersect( rownames(sce),csc$step5$BEAMgenelist))!=0){
    loginfo("get gene")
    fig5.10=plot_genes_branched_pseudotime(sce[intersect( rownames(sce),csc$step5$BEAMgenelist),],branch_point = csc$step5$pointid, color_by = "State", 
                                           ncol = 3)
    savepdf(paste(outdir,"step5","BEAM","custergene",sep="/"),fig5.10)}
  
  loginfo("save %s",paste0(outtemp,"step5.RDS"))
  saveRDS(sce,paste0(outtemp,"step5.RDS"))
  loginfo("step5 end")
  
}

if (isstep(6)){
  loginfo("step5 start")
  dir.create(paste(outdir,"step6",sep="/"))
  if (!isstep(5)){
    loginfo("input step3 rds")
    pbmcall=getrds(5)
  }
 library(cerebroApp) 
  pbmcall<- RunUMAP(pbmcall,  reduction.name = 'UMAP_3D',
                    reduction.key = 'UMAP3D_',
                    dims=1:csc$step3$dims,
                    n.components = 3,seed.use = 100
  )
  
  pbmcall <- BuildClusterTree(
    pbmcall,
    dims = 1:csc$step3$dims,
    reorder = FALSE,
    reorder.numeric = FALSE
  )
  pbmcall@misc$trees$seurat_clusters <- pbmcall@tools$BuildClusterTree
  
  pbmcall <- getMostExpressedGenes(
    pbmcall,
    assay = 'RNA',
    column_sample =  'orig.ident',
    column_cluster = 'seurat_clusters')
  


 getMarkerGenes_my<-function (object, assay = "RNA", organism = NULL, column_sample = "sample", 
                              column_cluster = "cluster", only_pos = TRUE, min_pct = 0.7, 
                              thresh_logFC = 0.25, thresh_p_val = 0.01, test = "wilcox", 
                              verbose = TRUE, ...) 
 {
   if (!requireNamespace("Seurat", quietly = TRUE)) {
     stop("Package 'Seurat' needed for this function to work. Please install it.", 
          call. = FALSE)
   }
   if (is.null(column_sample) | (column_sample %in% names(object@meta.data) == 
                                 FALSE)) {
     stop(paste0("Cannot find specified column (`object@meta.data$", 
                 column_sample, "`) that is supposed to contain sample information."), 
          call. = FALSE)
   }
   if (is.null(column_cluster) | (column_cluster %in% names(object@meta.data) == 
                                  FALSE)) {
     stop(paste0("Cannot find specified column (`object@meta.data$", 
                 column_cluster, "`) that is supposed to contain cluster information."), 
          call. = FALSE)
   }
   if (organism %in% c("hg", "mm") == FALSE) {
     message(paste0("[", format(Sys.time(), "%H:%M:%S"), 
                    "] No information about genes on cell surface because organism is ", 
                    "either not specified or not human/mouse."))
   }
   else {
     if (organism == "hg" || organism == "human") {
       temp_attributes <- "hgnc_symbol"
       temp_dataset <- "hsapiens_gene_ensembl"
     }
     else if (organism == "mm" || organism == "mouse") {
       temp_attributes <- "external_gene_name"
       temp_dataset <- "mmusculus_gene_ensembl"
     }
     attempt <- 1
     while (!exists("genes_on_cell_surface") && attempt <= 
            3) {
       try(genes_on_cell_surface <- biomaRt::getBM(attributes = temp_attributes, 
                                                   filters = "go", values = "GO:0009986", mart = biomaRt::useMart("ensembl", 
                                                                                                                  dataset = temp_dataset))[, 1])
     }
     if (!exists("genes_on_cell_surface")) {
       message(paste0("[", format(Sys.time(), "%H:%M:%S"), 
                      "] Genes in GO term \"cell surface\" (GO:0009986) could not be ", 
                      "retrieved, possibly due to the server not being reachable at the ", 
                      "moment."))
     }
   }
   if (is.factor(object@meta.data[[column_sample]])) {
     sample_names <- levels(object@meta.data[[column_sample]])
   }
   else {
     sample_names <- unique(object@meta.data[[column_sample]])
   }
   if (length(sample_names) > 1) {
     message(paste0("[", format(Sys.time(), "%H:%M:%S"), 
                    "] Get marker genes for ", length(sample_names), 
                    " samples..."))
     if (object@version < 3) {
       object <- Seurat::SetAllIdent(object, id = column_sample)
       object@ident <- factor(object@ident, levels = sample_names)
       if (utils::packageVersion("Seurat") < 3) {
         markers_by_sample <- Seurat::FindAllMarkers(object, 
                                                     only.pos = only_pos, min.pct = min_pct, thresh.use = thresh_logFC, 
                                                     return.thresh = thresh_p_val, test.use = test, 
                                                     print.bar = verbose, ...)
       }
       else {
         if ((assay %in% names(object@assays) == FALSE)) {
           stop(paste0("Assay slot `", assay, "` could not be found in provided Seurat ", 
                       "object."), call. = FALSE)
         }
         markers_by_sample <- Seurat::FindAllMarkers(object, 
                                                     assay = assay, only.pos = only_pos, min.pct = min_pct, 
                                                     logfc.threshold = thresh_logFC, return.thresh = thresh_p_val, 
                                                     test.use = test, verbose = verbose, ...)
       }
     }
     else {
       Seurat::Idents(object) <- factor(object@meta.data[[column_sample]], 
                                        levels = sample_names)
       if (utils::packageVersion("Seurat") < 3) {
         markers_by_sample <- Seurat::FindAllMarkers(object, 
                                                     only.pos = only_pos, min.pct = min_pct, thresh.use = thresh_logFC, 
                                                     return.thresh = thresh_p_val, test.use = test, 
                                                     print.bar = verbose, ...)
       }
       else {
         markers_by_sample <- Seurat::FindAllMarkers(object, 
                                                     assay = assay, only.pos = only_pos, min.pct = min_pct, 
                                                     logfc.threshold = thresh_logFC, return.thresh = thresh_p_val, 
                                                     test.use = test, verbose = verbose, ...)
       }
     }
     if (nrow(markers_by_sample) > 0) {
       markers_by_sample <- markers_by_sample %>% dplyr::select(c("cluster", 
                                                                  "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", 
                                                                  "p_val_adj")) %>% dplyr::rename(sample = .data$cluster)
       if (exists("genes_on_cell_surface")) {
         markers_by_sample <- markers_by_sample %>% dplyr::mutate(on_cell_surface = .data$gene %in% 
                                                                    genes_on_cell_surface)
       }
     }
     else {
       message(paste0("[", format(Sys.time(), "%H:%M:%S"), 
                      "] No marker genes found for any of the samples."))
       markers_by_sample <- "no_markers_found"
     }
   }
   else {
     message(paste0("[", format(Sys.time(), "%H:%M:%S"), 
                    "] Sample column provided but only 1 sample found."))
     markers_by_sample <- "no_markers_found"
   }
   if (is.factor(object@meta.data[[column_cluster]])) {
     cluster_names <- levels(object@meta.data[[column_cluster]])
   }
   else {
     cluster_names <- unique(object@meta.data[[column_cluster]]) %>% 
       sort()
   }
   if (length(cluster_names) > 1) {
     message(paste0("[", format(Sys.time(), "%H:%M:%S"), 
                    "] Get marker genes for ", length(cluster_names), 
                    " clusters..."))
     if (object@version < 3) {
       object <- Seurat::SetAllIdent(object, id = column_cluster)
       object@ident <- factor(object@ident, levels = cluster_names)
       if (utils::packageVersion("Seurat") < 3) {
         markers_by_cluster <- Seurat::FindAllMarkers(object, 
                                                      only.pos = only_pos, min.pct = min_pct, thresh.use = thresh_logFC, 
                                                      return.thresh = thresh_p_val, test.use = test, 
                                                      print.bar = verbose, ...)
       }
       else {
         if ((assay %in% names(object@assays) == FALSE)) {
           stop(paste0("Assay slot `", assay, "` could not be found in provided Seurat ", 
                       "object."), call. = FALSE)
         }
         markers_by_cluster <- Seurat::FindAllMarkers(object, 
                                                      assay = assay, only.pos = only_pos, min.pct = min_pct, 
                                                      logfc.threshold = thresh_logFC, return.thresh = thresh_p_val, 
                                                      test.use = test, verbose = verbose, ...)
       }
     }
     else {
       Seurat::Idents(object) <- factor(object@meta.data[[column_cluster]], 
                                        levels = cluster_names)
       if (utils::packageVersion("Seurat") < 3) {
         markers_by_cluster <- Seurat::FindAllMarkers(object, 
                                                      only.pos = only_pos, min.pct = min_pct, thresh.use = thresh_logFC, 
                                                      return.thresh = thresh_p_val, test.use = test, 
                                                      print.bar = verbose, ...)
       }
       else {
         markers_by_cluster <- Seurat::FindAllMarkers(object, 
                                                      assay = assay, only.pos = only_pos, min.pct = min_pct, 
                                                      logfc.threshold = thresh_logFC, return.thresh = thresh_p_val, 
                                                      test.use = test, verbose = verbose, ...)
       }
     }
     if (nrow(markers_by_cluster) > 0) {
       markers_by_cluster <- markers_by_cluster %>% dplyr::select(c("cluster", 
                                                                    "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", 
                                                                    "p_val_adj"))
       if (exists("genes_on_cell_surface")) {
         markers_by_cluster <- markers_by_cluster %>% 
           dplyr::mutate(on_cell_surface = .data$gene %in% 
                           genes_on_cell_surface)
       }
     }
     else {
       message(paste0("[", format(Sys.time(), "%H:%M:%S"), 
                      "] No marker genes found for any of the clusters."))
       markers_by_cluster <- "no_markers_found"
     }
   }
   else {
     message(paste0("[", format(Sys.time(), "%H:%M:%S"), 
                    "] Cluster column provided but only 1 cluster found."))
     markers_by_cluster <- "no_markers_found"
   }
   object@misc$marker_genes <- list(by_sample = markers_by_sample, 
                                    by_cluster = markers_by_cluster, parameters = list(only_positive = only_pos, 
                                                                                       minimum_percentage = min_pct, logFC_threshold = thresh_logFC, 
                                                                                       p_value_threshold = thresh_p_val, test = test))
   return(object)
 }
 
 pbmcall <- getMarkerGenes_my(
   pbmcall,
   assay = 'RNA',
   organism = ifelse(csc$step1$mttype=="mt","mm","hg"),
   column_sample = "orig.ident",
   column_cluster="seurat_clusters",
   name = 'cerebro_seurat',
   only_pos = F
 )
 exportFromSeurat(
   pbmcall,
   assay = "RNA",
   file=paste0(paste(outdir,"step6",sep="/"),"/cerebro_seurat",Sys.Date(),".crb"),
   experiment_name="sample",
   organism=csc$specis,
   column_sample = "orig.ident",
   column_cluster = "seurat_clusters",
   column_nUMI = "nCount_RNA",
   column_nGene = "nFeature_RNA",
   column_cell_cycle_seurat = NULL,
   column_cell_cycle_cyclone = NULL,
   add_all_meta_data = TRUE
 )
 loginfo("cerebroApp end")
 library(Matrix)
 dir.create(paste(outdir,"step6","circos",sep="/"))
 commanda=sprintf('Rscript %s -r %s -o %s -s gene -m %s',
                  csc$step6$circosbin,
 paste0(outdir,"/temp/step4.RDS"),
 paste(outdir,"step6","circos",".",sep="/"),
 paste0(outdir,"/step4/clusterall_genelist.xls")
 )
 syscomd(commanda,savef = paste(outdir,"shell","05.circos.sh",sep = "/"))
 commanda=sprintf( "perl %s  -t %s -c %s -m %s -d  %s -p result -a %s",
                   csc$step6$circos_perl_bin,
                   paste(outdir,"step6","circos","type2gene.txt",sep="/"),
                   paste(outdir,"step6","circos","gene.cellpercent.txt",sep="/"),
                   paste(outdir,"step6","circos","gene.ClusterMean.txt",sep="/"),
                   paste(outdir,"step6","circos","result",sep="/"),
                   #paste0(outdir,"/step4/clusterall_genelist.xls")
                   paste(outdir,"step6","circos","test",sep="/")
                   
                   )
 syscomd(commanda,savef = paste(outdir,"shell","05.circos.sh",sep = "/"))
 
  }
if(isstep(1)){
  commanda=sprintf("mv %s %s",paste(outdir,"00.cellranger",sep="/"),paste(outdir,"02.Cellranger",sep="/"))
  syscomd(commanda,savef = paste(outdir,"shell","06.all.sh",sep = "/"))
  commanda=sprintf("mv %s %s",paste(outdir,"00.fastp",sep="/"),paste(outdir,"01.Fastp",sep="/"))
  syscomd(commanda,savef = paste(outdir,"shell","06.all.sh",sep = "/"))
  commanda=sprintf("mv %s %s",paste(outdir,"step1",sep="/"),paste(outdir,"03.CellFilter",sep="/"))
  syscomd(commanda,savef = paste(outdir,"shell","06.all.sh",sep = "/"))
}
if(isstep(3)){
  commanda=sprintf("mv %s %s",paste(outdir,"step3",sep="/"),paste(outdir,"04.PCA_UMAP",sep="/"))
  syscomd(commanda,savef = paste(outdir,"shell","06.all.sh",sep = "/"))
}
if(isstep(4)){
  commanda=sprintf("mv %s %s",paste(outdir,"step4",sep="/"),paste(outdir,"05.MarkerGene",sep="/"))
  syscomd(commanda,savef = paste(outdir,"shell","06.all.sh",sep = "/"))
}
if(isstep(5)){
  commanda=sprintf("mv %s %s",paste(outdir,"step5",sep="/"),paste(outdir,"06.Pseudotime",sep="/"))
  syscomd(commanda,savef = paste(outdir,"shell","06.all.sh",sep = "/"))
}
if(isstep(6)){
  commanda=sprintf("mv %s %s",paste(outdir,"step6",sep="/"),paste(outdir,"07.Cerebro",sep="/"))
  syscomd(commanda,savef = paste(outdir,"shell","06.all.sh",sep = "/"))
}
if(isstep(7)){
  commanda=sprintf("mv %s %s",paste(outdir,"step7",sep="/"),paste(outdir,"08.Circos",sep="/"))
  syscomd(commanda,savef = paste(outdir,"shell","06.all.sh",sep = "/"))
}
if(isstep(10)){
  
}

if(isstep(10)){
  loginfo("step6 start")
  dir.create(paste(outdir,"step6",sep="/"))
  if (!isstep(4)){
    loginfo("input step4 rds")
    pbmcall=getrds(4)
  }
  
  library(CellChat)
  library(msigdbr)
  library(fgsea)
  library(ggalluvial)
  library(sys)
  #  outdir="/haplox/haprs/liyq/singlecell/allliver/getout/step5/test"
  #  pbmcall=readRDS("/haplox/haprs/liyq/singlecell/allliver/step4.3.RDS")
  data.input <- GetAssayData(pbmcall, assay = "SCT", slot = "data")
  labels <- Idents(pbmcall)
  meta=pbmcall@meta.data
  cell.use = rownames(meta)
  data.input = data.input[, cell.use]
  meta = meta[cell.use, ]
  #  meta$labels=cellchat@meta$labels
  
  meta$labels= paste0("cluster",meta$seurat_clusters)
  
  
  
  
  
  cellchat=createCellChat(object = data.input, meta = meta, group.by = "labels")
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  dplyr::glimpse(CellChatDB$interaction)
  
  # showDatabaseCategory(CellChatDB)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  CellChatDB.use <- CellChatDB 
  cellchat@DB <- CellChatDB.use
  
  #showDatabaseCategory(CellChatDB)
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 8)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # project gene expression data onto PPI network (optional)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dir.create(paste(outdir,"step5","net_circle",sep="/"),recursive = T)
  pdf(paste(outdir,"step5","net_circle","topgeneheatmap.pdf",sep="/"))
  
  dev.off()
  
  
  
  mat <- cellchat@net$weight
  allname="B cells
Dendritic cells
Endothelial cells
Fibroblasts
Hepatocytes
Macrophage 1
Macrophage 2
Mix 1
Mix 2
Mix 3
Mix 4
Myeloid cells
Naive T cells
Neutrophils
NK cells
NK2
Plasma B cells
Plasmacytoid DC
T cells
T cells 2"
  levels(cellchat@idents)=geneliss
  #cellchat@meta$labels=cellchat@idents
  #cellchat <- setIdent(cellchat, ident.use = "cell_types") 
  #geneliss=str_split(allname,"\n")[[1]]
  #rownames(mat)=geneliss
  #colnames(mat)=geneliss
  par(mfrow = c(3,4), xpd=TRUE)
  c(8:11,1:7,12:20)
  mat=mat[c(1:3,8:11,4:7,12:20),c(1:3,8:11,4:7,12:20)]
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    pdf(paste(outdir,"step5","net_circle",paste0(rownames(mat2)[i],"_circle.pdf"),sep ="/" ),width=5,height = 5)
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = "T cell",label.edge=T)
    dev.off()
    
  }
  
  mat2[c(13,19,20), ] <- mat[c(13,19,20), ]
  mat2=mat[,c(13,19,20) ]
  
  
  pathways.show <- c("CXCL") 
  vertex.receiver = seq(1,9) # a numeric vector. 
  netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout ="hierarchy")
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  par(mfrow=c(1,1))
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  LR.show <- pairLR.CXCL[1,]
  netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show,vertex.receiver = vertex.receiver, layout = "hierarchy")
  
  pathways.show.all <- cellchat@netP$pathways
  # check the order of cell identity to set suitable vertex.receiver
  levels(cellchat@idents)
  vertex.receiver = seq(1,8)
  dir.create(paste(outdir,"step5","net_circle","pathway",sep="/"),recursive = T)
  
  for (i in 1:length(pathways.show.all)) {
    # Visualize communication network associated with both signaling pathway and individual L-R pairs
    #  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
    # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
    gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
    ggsave(filename=paste(paste0(outdir,"step5","net_circle","pathway",pathways.show.all[i], "_L-R_contribution.pdf")), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
    pdf(paste(outdir,"step5","net_circle","pathway",paste0(pathways.show.all[i], "_hier.pdf"),sep="/"),width = 14,height = 4)
    netVisual_aggregate(cellchat, signaling = pathways.show.all[i], 
                        signaling.name="",
                        vertex.receiver = vertex.receiver,layout ="hierarchy",vertex.label.cex=0.8,label.dist=7)
    dev.off()
    pdf(paste(outdir,"step5","net_circle","pathway",paste0(pathways.show.all[i], "_circle.pdf"),sep="/"))
    
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
    dev.off()
  }
  
  dir.create(paste(outdir,"step5","net_circle","bubble",sep="/"),recursive = T)
  for(i in 1:nrow(mat)){
    #  pdf(paste(outdir,"step5","net_circle","bubble",paste0(rownames(mat2)[i], "_bubble.pdf"),sep="/"),width = 7,height = 7)
    gg=netVisual_bubble(cellchat,sources.use=i,remove.isolate = FALSE,font.size=15,font.size.title=15,angle.x=45)+
      theme(legend.text = element_text(size = 12,
                                       face = "bold"),
            legend.title  = element_text( size = 12,
                                          face = "bold"))
    # Sys.sleep(1)
    #  dev.off()
    # Sys.sleep(1)
    ggsave(filename=paste(outdir,"step5","net_circle",paste0(rownames(mat)[i], "_bubble.pdf"),sep="/"), plot=gg, width = 8, height = 8, units = 'in', dpi = 300)
  }
  
  plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  saveRDS(cellchat ,paste(outdir,"step5","cellchat.rds",sep = "/"))
  
  
}


