## Initialize packages
Sys.setenv("DISPLAY"=":0.0")
options(bitmapType="cairo-png")
suppressPackageStartupMessages(library(clusterProfiler)) 
suppressPackageStartupMessages(library(ReactomePA)) 
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(org.Mm.eg.db))
#suppressPackageStartupMessages(library(org.Rn.eg.db))
#suppressPackageStartupMessages(library(org.Ss.eg.db))
#suppressPackageStartupMessages(library(org.Dm.eg.db))
#suppressPackageStartupMessages(library(org.At.tair.db))
#suppressPackageStartupMessages(library(org.Gg.eg.db))
#suppressPackageStartupMessages(library(org.Mmu.eg.db))
#suppressPackageStartupMessages(library(org.Pt.eg.db))
#suppressPackageStartupMessages(library(org.Bt.eg.db))
#suppressPackageStartupMessages(library(org.Ce.eg.db))
#suppressPackageStartupMessages(library(org.Cf.eg.db))
#suppressPackageStartupMessages(library(org.Dr.eg.db))
#suppressPackageStartupMessages(library(org.Sc.sgd.db))
suppressPackageStartupMessages(library(pathview))
library(RColorBrewer)
suppressPackageStartupMessages(library(ggplot2))
library("getopt")
library("Cairo")
command=matrix(c(
"file","f",1,"character",
"name","n",1,"character",
"species","s",1,"character",
"databases","d",1,"character",
"outdir","o",1,"character",
"gene.col","g",1,"character",
"idtype","t",1,"character",
"log2fc.col","L",1,"character",
"log2fc.cut","l",1,"character",
"padj.col","Q",1,"character",
"gene.padj","q",1,"character",
"p.col","P",1,"character",
"gene.p","p",1,"character",
"padj.cut","C",1,"character",
"p.cut","c",1,"character",
"bgall","a",1,"character",
"bggene","b",1,"character",
"help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if (!is.null(args$help) || is.null(args$file) || is.null(args$name)) {
    cat(paste(getopt(command, usage = T), "\n"))
    q()
}

## Initialize data
if(!file.exists(args$outdir)){
    dir.create(args$outdir) 
}
setwd(args$outdir)
## Set organism
spe_db=unlist(strsplit(args$species, "[,]"))
db_type=unlist(strsplit(args$databases, "[,]"))
OrgDb = spe_db[1]
organism = spe_db[2]
rorganism = spe_db[3]
## Get Gene Table
mydata<-read.table(args$file,header=TRUE,sep="\t",fill=TRUE)
a=as.numeric(args$gene.col)
colnames(mydata)[a]<-"AccID"
if(!is.null(args$log2fc.col)){b=as.numeric(args$log2fc.col)}
if(!is.null(args$padj.col)){c=as.numeric(args$padj.col)}
if(!is.null(args$p.col)){d=as.numeric(args$p.col)}

## Use the data frame depend on the input args
if(is.null(args$log2fc.cut) && is.null(args$gene.padj) && is.null(args$gene.p)){ #just gene col
    nowdata = mydata[,a]
    nowdata<-as.data.frame(nowdata)
    colnames(nowdata)<-"AccID"
}else if(!is.null(args$log2fc.cut)){ # have log2fc.cut parameters
    if(!is.null(args$gene.padj)){ 
        nowdata =mydata[,c(a,b,c)]
        colnames(nowdata)<-c("AccID","log2FC","p.adjust")
        nowdata =nowdata[abs(nowdata$log2FC) >= as.numeric(args$log2fc.cut),]
        nowdata =nowdata[nowdata$p.adjust <= as.numeric(args$gene.padj),]
    }else if(is.null(args$gene.padj) && !is.null(args$gene.p)){
        nowdata =mydata[,c(a,b,d)]
        colnames(nowdata)<-c("AccID","log2FC","p.adjust")
        nowdata =nowdata[abs(nowdata$log2FC) >= as.numeric(args$log2fc.cut),]
        nowdata =nowdata[nowdata$p.adjust <= as.numeric(args$gene.p),]
    }else if(is.null(args$gene.padj) && is.null(args$gene.p)){
        nowdata =mydata[,c(a,b)]
        colnames(nowdata)<-c("AccID","log2FC")
        nowdata =nowdata[abs(nowdata$log2FC) >= as.numeric(args$log2fc.cut),]
    }
}else if(!is.null(args$gene.padj) && is.null(args$log2fc.cut)){
    nowdata =mydata[,c(a,c)]
    colnames(nowdata)<-c("AccID","p.adjust")
    nowdata =nowdata[nowdata$p.adjust <= as.numeric(args$gene.padj),]
}else if(is.null(args$gene.padj) && !is.null(args$gene.p) && is.null(args$log2fc.cut)){
    nowdata =mydata[,c(a,d)]
    colnames(nowdata)<-c("AccID","p.adjust")
    nowdata =nowdata[nowdata$p.adjust <= as.numeric(args$gene.p),]
}
## Gene ID Trans
if(args$idtype=="ENTREZID"){
    bi<-nowdata
    colnames(bi)[1]<-"ENTREZID"
}else{
    #print(nowdata$AccID)
    bi<-bitr(as.character(nowdata$AccID), fromType = args$idtype,toType = "ENTREZID",OrgDb = OrgDb, drop = TRUE)
    print("****check id")
    print(head(bi,3))
}

## Read background gene data
if(args$bgall=="false" && !is.null(args$bggene)){
    bggene<-read.table(args$bggene,header=TRUE,sep="\t",fill=TRUE)
    bggene = bggene[,a]
    bggene<-as.data.frame(bggene)
    colnames(bggene)<-"AccID"
    if(args$idtype=="ENTREZID"){
        bibggene<-bggene
        colnames(bibggene)<-"ENTREZID"
    }else{
        bibggene<-bitr(as.character(bggene$AccID), fromType = args$idtype,toType = "ENTREZID",OrgDb = OrgDb, drop = TRUE)
    }
}

## Function Enrich Analasys
if(!is.null(bi) && nrow(bi)>1){
    ## use all Gene as background
    if(args$bgall=="false" && (!is.null(bibggene) || !is.na(as.data.frame(bibggene)[1,1]))){
        ## GO Enrich & plot barplot & DAG
        universe=as.factor(bibggene$ENTREZID)
        egoall<-enrichGO(gene = bi$ENTREZID,universe=universe,OrgDb=OrgDb,ont="ALL",readable = TRUE,pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1)
        write.table(egoall,paste(args$name,".GO_Enrich.xls",sep=""),sep="\t",quote =FALSE,row.names=FALSE)
	print("good1")
	## GO CC
	## one of padj.cut, p.cut must be set.
        if(!is.null(args$padj.cut)){
            egoCC<-enrichGO(gene = bi$ENTREZID,universe=universe,OrgDb=OrgDb,ont ="CC",readable = TRUE,pAdjustMethod = "BH",pvalueCutoff = as.numeric(args$padj.cut),qvalueCutoff =as.numeric(args$padj.cut))
	    print("go cc good")
	    print(head(egoCC,3))
        }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
            egoCC<-enrichGO(gene = bi$ENTREZID,universe=universe,OrgDb=OrgDb,ont ="CC",readable = TRUE,pAdjustMethod = "none",pvalueCutoff = as.numeric(args$p.cut),qvalueCutoff =1)
        }
	print(egoCC)
        if(!is.null(egoCC) && !is.na(as.data.frame(egoCC)[1,1])){
	    print("check1")
            CairoSVG(paste(args$name,".CC.DAG.svg",sep=""),width=8, height=6)
            p<-plotGOgraph(egoCC)
            print(p)
	    print("cc_dag.svg.good!")
            dev.off()
        }else{
	    print("0 enriched go cc terms found")
	}
	## GO MF
        if(!is.null(args$padj.cut)){
            egoMF<-enrichGO(gene = bi$ENTREZID,universe=universe,OrgDb=OrgDb,ont ="MF",readable = TRUE,pAdjustMethod = "BH",pvalueCutoff = as.numeric(args$padj.cut),qvalueCutoff =as.numeric(args$padj.cut))
        }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
            egoMF<-enrichGO(gene = bi$ENTREZID,universe=universe,OrgDb=OrgDb,ont ="MF",readable = TRUE,pAdjustMethod = "none",pvalueCutoff = as.numeric(args$p.cut),qvalueCutoff =1)
        }
	print(egoMF)
        if(!is.null(egoMF) && !is.na(as.data.frame(egoMF)[1,1])){
            CairoSVG(paste(args$name,".MF.DAG.svg",sep=""),width=8, height=6)
            p<-plotGOgraph(egoMF)
            print(p)
            dev.off()
        }else{
	    print("0 enriched go MF terms found")
	}
        ## GO BP
	if(!is.null(args$padj.cut)){
            egoBP<-enrichGO(gene = bi$ENTREZID,universe=universe,OrgDb=OrgDb,ont ="BP",readable = TRUE,pAdjustMethod = "BH",pvalueCutoff = as.numeric(args$padj.cut),qvalueCutoff =as.numeric(args$padj.cut))
        }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
            egoBP<-enrichGO(gene = bi$ENTREZID,universe=universe,OrgDb=OrgDb,ont ="BP",readable = TRUE,pAdjustMethod = "none",pvalueCutoff = as.numeric(args$p.cut),qvalueCutoff =1)
        }
	print(egoBP)
        if(!is.null(egoBP) && !is.na(as.data.frame(egoBP)[1,1])){
            CairoSVG(paste(args$name,".BP.DAG.svg",sep=""),width=8, height=6)
            p<-plotGOgraph(egoBP)
            print(p)
            dev.off()
        }else{
	    print(" enriched go BF terms found")
	}
        gomax=50
        picture_data<-data.frame(Description=0,GeneRatio=0,pic=0,type=0)
        picture_data<-picture_data[-1,]
        print("check2")
	if(!is.null(egoBP) && !is.na(as.data.frame(egoBP)[1,1])){
	    print("here1")
            if(!is.null(args$padj.cut)){
                egoBP<-egoBP[egoBP$p.adjust <= as.numeric(args$padj.cut),]
            }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
                egoBP<-egoBP[egoBP$pvalue <= as.numeric(args$p.cut),]
            }
            t<-egoBP[order(-egoBP[,9]),]
            picBP<-t[1:10,c(2,3)]
            for(i in seq(dim(picBP)[1])){
                tmp = length(unlist(strsplit(as.character(picBP[i,1]),split="")))
                if(tmp>gomax){
                    gomax=tmp
                }
            }
            tmp<-unlist(lapply(picBP[,2],function(x){tmp<-as.numeric(strsplit(x,"/",fixed=TRUE)[[1]][1])/as.numeric(strsplit(x,"/",fixed=TRUE)[[1]][2])}))
            picBP$pic<-tmp
            picBP$type<-"BP"
            picture_data<-rbind(picture_data,picBP)
	    print(picture_data)
        }
        if(!is.null(egoCC) && !is.na(as.data.frame(egoCC)[1,1])){
	    print("here2")
            if(!is.null(args$padj.cut)){
                egoCC<-egoCC[egoCC$p.adjust <= as.numeric(args$padj.cut),]
            }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
                egoCC<-egoCC[egoCC$pvalue <= as.numeric(args$p.cut),]
            }
            t<-egoCC[order(-egoCC[,9]),]
            picCC<-t[1:10,c(2,3)]
            for(i in seq(dim(picCC)[1])){
                tmp = length(unlist(strsplit(as.character(picCC[i,1]),split="")))
                if(tmp>gomax){
                    gomax=tmp
                }
            }
            tmp<-unlist(lapply(picCC[,2],function(x){tmp<-as.numeric(strsplit(x,"/",fixed=TRUE)[[1]][1])/as.numeric(strsplit(x,"/",fixed=TRUE)[[1]][2])}))
            picCC$pic<-tmp
            picCC$type<-"CC"
            picture_data<-rbind(picture_data,picCC)
	    print(picture_data)
        }
        if(!is.null(egoMF) && !is.na(as.data.frame(egoMF)[1,1])){
	    print("here3")
            if(!is.null(args$padj.cut)){
                egoMF<-egoMF[egoMF$p.adjust <= as.numeric(args$padj.cut),]
            }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
                egoMF<-egoMF[egoMF$pvalue <= as.numeric(args$p.cut),]
            }
            t<-egoMF[order(-egoMF[,9]),]
            picMF<-t[1:10,c(2,3)]
            for(i in seq(dim(picMF)[1])){
                tmp = length(unlist(strsplit(as.character(picMF[i,1]),split="")))
                if(tmp>gomax){
                    gomax=tmp
                }
            }
            tmp<-unlist(lapply(picMF[,2],function(x){tmp<-as.numeric(strsplit(x,"/",fixed=TRUE)[[1]][1])/as.numeric(strsplit(x,"/",fixed=TRUE)[[1]][2])}))
            picMF$pic<-tmp
            picMF$type<-"MF"
            picture_data<-rbind(picture_data,picMF)
	    print(picture_data)
        }
	print("check3")
        picture_data=na.omit(picture_data)
        picture_data$Description<-as.factor(picture_data$Description)
        picture_data$type<-as.factor(picture_data$type)
        go_width=12+(gomax-50)/35
        go_height=10+(gomax-50)/35
        if(!is.null(picture_data) && nrow(picture_data) > 0){
            r<-brewer.pal(3,"Set2")
            g1<-ggplot(picture_data,aes(x=picture_data$Description,y=picture_data$pic))+geom_bar(stat="identity",aes(fill=picture_data$type))+theme_bw()+ggtitle(paste(args$name,".GO_Enrich",sep=""))+labs(fill="Type")+xlab("GO term")+ylab("GeneRatio")+scale_fill_manual(values = r)+theme(panel.grid.major=element_line(colour=NA))+facet_wrap(~picture_data$type,scales="free")+theme(panel.border=element_blank())+theme(strip.background = element_rect(colour="white", fill="white"))+theme(plot.margin = unit(c(1,1,1,4),"cm"))+theme(title=element_text(size=12,face="bold"),axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),axis.title.y=element_text(size=12,face="bold"),axis.text.y=element_text(size=10),legend.title = element_text(size = 12,face = "bold"),legend.text=element_text(size = 10),plot.title=element_text(hjust=0.5))
	   print("good2")
           ggsave(g1,file=paste(args$name,".go.png",sep=""),  device="png", width=go_width, height=go_height, dpi=700)
           ggsave(g1,file=paste(args$name,".go.pdf",sep=""),  device="pdf", width=go_width, height=go_height)
        }
    }else{
        ## GO Enrich & plot barplot & DAG
        egoall<-enrichGO(gene = bi$ENTREZID,OrgDb=OrgDb,ont="ALL",readable = TRUE,pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1)
        write.table(egoall,paste(args$name,".GO_Enrich.xls",sep=""),sep="\t",quote =FALSE,row.names=FALSE)
	print("good3")
        if(!is.null(args$padj.cut)){
            egoCC<-enrichGO(gene = bi$ENTREZID,OrgDb=OrgDb,ont ="CC",readable = TRUE,pAdjustMethod = "BH",pvalueCutoff = as.numeric(args$padj.cut),qvalueCutoff =as.numeric(args$padj.cut)) 
        }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
            egoCC<-enrichGO(gene = bi$ENTREZID,OrgDb=OrgDb,ont ="CC",readable = TRUE,pAdjustMethod = "BH",pvalueCutoff = as.numeric(args$p.cut),qvalueCutoff =1)
        }
        if(!is.null(egoCC) && !is.na(as.data.frame(egoCC)[1,1])){
            CairoSVG(paste(args$name,".CC.DAG.svg",sep=""),width=8, height=6)
            p<-plotGOgraph(egoCC)
            print(p)       
            dev.off()
        }
        if(!is.null(args$padj.cut)){
            egoMF<-enrichGO(gene = bi$ENTREZID,OrgDb=OrgDb,ont ="MF",readable = TRUE,pAdjustMethod = "BH",pvalueCutoff = as.numeric(args$padj.cut),qvalueCutoff =as.numeric(args$padj.cut))
        }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
            egoMF<-enrichGO(gene = bi$ENTREZID,OrgDb=OrgDb,ont ="MF",readable = TRUE,pAdjustMethod = "BH",pvalueCutoff = as.numeric(args$p.cut),qvalueCutoff =1)
        }
        if(!is.null(egoMF) && !is.na(as.data.frame(egoMF)[1,1])){
            CairoSVG(paste(args$name,".MF.DAG.svg",sep=""),width=8, height=6)
            p<-plotGOgraph(egoMF)
            print(p)
            dev.off()
        }
        if(!is.null(args$padj.cut)){
            egoBP<-enrichGO(gene = bi$ENTREZID,OrgDb=OrgDb,ont ="BP",readable = TRUE,pAdjustMethod = "BH",pvalueCutoff = as.numeric(args$padj.cut),qvalueCutoff =as.numeric(args$padj.cut))
        }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
            egoBP<-enrichGO(gene = bi$ENTREZID,OrgDb=OrgDb,ont ="BP",readable = TRUE,pAdjustMethod = "BH",pvalueCutoff = as.numeric(args$p.cut),qvalueCutoff =1)
        }
        if(!is.null(egoBP) && !is.na(as.data.frame(egoBP)[1,1])){
            CairoSVG(paste(args$name,".BP.DAG.svg",sep=""),width=8, height=6)
            p<-plotGOgraph(egoBP)
            print(p)
            dev.off()
        }
        gomax=50
        picture_data<-data.frame(Description=0,GeneRatio=0,pic=0,type=0)
        picture_data<-picture_data[-1,]
        if(!is.null(egoBP) && !is.na(as.data.frame(egoBP)[1,1])){    
            if(!is.null(args$padj.cut)){
                egoBP<-egoBP[egoBP$p.adjust <= as.numeric(args$padj.cut),]
            }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
                egoBP<-egoBP[egoBP$pvalue <= as.numeric(args$p.cut),]
            }
            t<-egoBP[order(-egoBP[,9]),]
            picBP<-t[1:10,c(2,3)]
            for(i in seq(dim(picBP)[1])){
                tmp = length(unlist(strsplit(as.character(picBP[i,1]),split="")))
                if(tmp>gomax){
                    gomax=tmp
                }
            }
            tmp<-unlist(lapply(picBP[,2],function(x){tmp<-as.numeric(strsplit(x,"/",fixed=TRUE)[[1]][1])/as.numeric(strsplit(x,"/",fixed=TRUE)[[1]][2])}))
            picBP$pic<-tmp
            picBP$type<-"BP"    
            picture_data<-rbind(picture_data,picBP)
        }    
        if(!is.null(egoCC) && !is.na(as.data.frame(egoCC)[1,1])){
            if(!is.null(args$padj.cut)){
                egoCC<-egoCC[egoCC$p.adjust <= as.numeric(args$padj.cut),]
            }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
                egoCC<-egoCC[egoCC$pvalue <= as.numeric(args$p.cut),]
            }
            t<-egoCC[order(-egoCC[,9]),]
            picCC<-t[1:10,c(2,3)]
            for(i in seq(dim(picCC)[1])){
                tmp = length(unlist(strsplit(as.character(picCC[i,1]),split="")))
                if(tmp>gomax){
                   gomax=tmp
                }
            }
            tmp<-unlist(lapply(picCC[,2],function(x){tmp<-as.numeric(strsplit(x,"/",fixed=TRUE)[[1]][1])/as.numeric(strsplit(x,"/",fixed=TRUE)[[1]][2])}))
            picCC$pic<-tmp
            picCC$type<-"CC"
            picture_data<-rbind(picture_data,picCC)
        }
        if(!is.null(egoMF) && !is.na(as.data.frame(egoMF)[1,1])){
            if(!is.null(args$padj.cut)){
                egoMF<-egoMF[egoMF$p.adjust <= as.numeric(args$padj.cut),]
            }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
                egoMF<-egoMF[egoMF$pvalue <= as.numeric(args$p.cut),]
            }
            t<-egoMF[order(-egoMF[,9]),]
            picMF<-t[1:10,c(2,3)]
            for(i in seq(dim(picMF)[1])){
                tmp = length(unlist(strsplit(as.character(picMF[i,1]),split="")))
                if(tmp>gomax){
                    gomax=tmp
                }
            }
            tmp<-unlist(lapply(picMF[,2],function(x){tmp<-as.numeric(strsplit(x,"/",fixed=TRUE)[[1]][1])/as.numeric(strsplit(x,"/",fixed=TRUE)[[1]][2])}))
            picMF$pic<-tmp
            picMF$type<-"MF"
            picture_data<-rbind(picture_data,picMF)
        }
        picture_data=na.omit(picture_data)
        picture_data$Description<-as.factor(picture_data$Description)
        picture_data$type<-as.factor(picture_data$type)
        go_width=12+(gomax-50)/35
        go_height=10+(gomax-50)/35
        if(!is.null(picture_data) && nrow(picture_data) > 0){
            r<-brewer.pal(3,"Set2")
            g1<-ggplot(picture_data,aes(x=picture_data$Description,y=picture_data$pic))+geom_bar(stat="identity",aes(fill=picture_data$type))+theme_bw()+ggtitle(paste(args$name,".GO_Enrich",sep=""))+labs(fill="Type")+xlab("GO term")+ylab("GeneRatio")+scale_fill_manual(values = r)+theme(panel.grid.major=element_line(colour=NA))+facet_wrap(~picture_data$type,scales="free")+theme(panel.border=element_blank())+theme(strip.background = element_rect(colour="white", fill="white"))+theme(plot.margin = unit(c(1,1,1,4),"cm"))+theme(title=element_text(size=12,face="bold"),axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),axis.title.y=element_text(size=12,face="bold"),axis.text.y=element_text(size=10),legend.title = element_text(size = 12,face = "bold"),legend.text=element_text(size = 10),plot.title=element_text(hjust=0.5))
           ggsave(g1,file=paste(args$name,".go.png",sep=""),  device="png", width=go_width, height=go_height, dpi=700)
           ggsave(g1,file=paste(args$name,".go.pdf",sep=""),  device="pdf", width=go_width, height=go_height)
        }
    }
    
## KEGG Enrich & plot
db_ava <- c("KEGG","BIOCARTA","BioCyc","PANTHER","PID","Reactome")
for (db in db_type){
if(db %in% db_ava){ #added by wwf
    kegg_db<-read.table(paste("/hapbin/users/liyq/pipline/singlecell/bin/singleRdata/clusterprofiler/",db,"_hsa_entrezID.xls",sep=""),header=T,sep="\t")
    print("*****database check")
    print(head(kegg_db,3))
    path2gene<-kegg_db[,c(1,2)]
    path2name<-kegg_db[,c(1,3)]
    n=0
    for(i in 1:nrow(bi)){
        if(bi$ENTREZID[i] %in% path2gene$GENE){
            n=n+1
        }
    }
    if(n>0){
        kkall<-enricher(gene=bi$ENTREZID,pvalueCutoff = 1,qvalueCutoff = 1,pAdjustMethod = "BH",TERM2GENE=path2gene,TERM2NAME=path2name)
        if (is.null(kkall)){
		print(paste("NOTE:No result of ",db))
		next
	}
	#print("***********check write")
	## renameID
	#kkall<- rename(kkall,c("ID"="#ID")
	write.table(kkall,paste(args$name,".",db,"_Enrich.xls",sep=""),sep="\t",quote =FALSE,row.names=FALSE)
        #print("***********write done")
	kegg<-read.table(paste(args$name,".",db,"_Enrich.xls",sep=""),header=T,sep="\t",fill=TRUE)
        #print("***********check read done")
	if(!is.null(args$padj.cut)){
	    #print("maybe here1")
            kegg<-kegg[kegg$p.adjust <= as.numeric(args$padj.cut),]
	    print("**************DataFrame after filter by padj.cut")
	    #print(dim(kegg))
	    #print("maybe here2")
        }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
            kegg<-kegg[kegg$pvalue <= as.numeric(args$p.cut),]
        }
        if(!is.null(nrow(kegg)) && nrow(kegg)>0 ){
	    #print("check good here1")
            kegg<-kegg[order(-kegg[,9]),]
	    #print("check good here2")
            kegg<-head(kegg,n=20L)
            #print("check good here3")
	    keggmax=50
            for(i in seq(dim(kegg)[1])){
	    	print("check good here3")
                tmp = length(unlist(strsplit(as.character(kegg[i,2]),split="")))
                if(tmp>keggmax){
                    keggmax<-tmp
                }
            }
	    print("check good here4")
            kegg_width=7+(keggmax-50)/20
	    print("check good here5")
            diff_gene<-as.numeric(unlist(strsplit(as.character(kegg$GeneRatio[1]),"/"))[2])
            print("check good here6")
	    kegg$GeneRatio<-kegg$Count/diff_gene
	    print("check good here7")
            p<-ggplot(kegg, aes(GeneRatio,Description))+theme_bw()+theme(panel.background = element_blank(),panel.border=element_rect(fill="transparent", color="black"))+geom_point(aes(colour=p.adjust,size=Count))+scale_colour_gradientn(colours=c("red","green")) +expand_limits(color=seq(0, 0.05, by=0.01))+ggtitle(paste(args$name,".",db,"_Enrich",sep="")) + xlab("GeneRatio") +ylab("")+theme(title=element_text(size=12,face="bold"),axis.text=element_text( color="black", size=10),axis.title.x=element_text(size=12),axis.text.y=element_text(size=12),legend.title = element_text(size=10,face="bold"),legend.text=element_text(size = 8),plot.title=element_text(hjust=0.5),legend.key=element_blank())
            ggsave(p,file=paste(args$name,".",db,".png",sep=""),  device="png", width=kegg_width, height=7, dpi=700)
            ggsave(p,file=paste(args$name,".",db,".pdf",sep=""),  device="pdf", width=kegg_width, height=7)
        }
    }
    }}# added by wwf
    ## Reactome Enrich & plot
    rtmp<-nowdata[,1]
    rtmp<-as.data.frame(rtmp)
    names(rtmp)<-args$idtype
    rmege<-merge(x=bi,y=rtmp,by.x=args$idtype,by.y=args$idtype,incomparables = NA)
    if(!is.na(rorganism)){
        rrall<-enrichPathway(as.character(rmege$ENTREZID),organism = rorganism,pvalueCutoff = 1,pAdjustMethod = "BH",readable = TRUE)
        if(!is.null(rrall)){
            write.table(rrall,paste(args$name,".Reactome_Enrich.xls",sep=""),sep="\t",quote =FALSE,row.names=FALSE)
            reac<-read.table(paste(args$name,".Reactome_Enrich.xls",sep=""),header=T,sep="\t")
            if(!is.null(args$padj.cut)){
                reac<-reac[reac$p.adjust <= as.numeric(args$padj.cut),]
            }else if(!is.null(args$p.cut) && is.null(args$padj.cut)){
                reac<-reac[reac$pvalue <= as.numeric(args$p.cut),]
            }
            if(!is.null(nrow(reac)) && nrow(reac)>0 ){
                reac<-reac[order(-reac[,9]),]
                reac<-head(reac,n=20L)
                reacmax<-50
                for(i in seq(dim(reac)[1])){
                    tmp = length(unlist(strsplit(as.character(reac[i,2]),split="")))
                    if(tmp>reacmax){
                        reacmax<-tmp
                    }
                }
                reac_width=8+(reacmax-50)/15
                diff_gene<-as.numeric(unlist(strsplit(as.character(reac$GeneRatio[1]),"/"))[2])
                reac$GeneRatio<-reac$Count/diff_gene
                p<-ggplot(reac, aes(GeneRatio,Description))+theme_bw()+theme(panel.background = element_blank(),panel.border=element_rect(fill="transparent", color="black"))+geom_point(aes(colour=p.adjust,size=Count))+scale_colour_gradientn(colours=c("red","green")) +expand_limits(color=seq(0, 0.05, by=0.01))+ggtitle(paste(args$name,".Reactome_Enrich",sep="")) + xlab("GeneRatio") +ylab("")+theme(title=element_text(size=12,face="bold"),axis.text=element_text( color="black", size=10),axis.title.x=element_text(size=12),axis.text.y=element_text(size=12),legend.title = element_text(size=10,face="bold"),legend.text=element_text(size = 8),plot.title=element_text(hjust=0.5),legend.key=element_blank())
                ggsave(p,file=paste(args$name,".reactome.png",sep=""),  device="png", width=reac_width, height=7, dpi=700)
                ggsave(p,file=paste(args$name,".reactome.pdf",sep=""),  device="pdf", width=reac_width, height=7)
            }
        }
    }
}
