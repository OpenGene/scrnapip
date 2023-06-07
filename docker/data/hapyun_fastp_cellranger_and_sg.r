suppressPackageStartupMessages(library(optparse))
library(RcppTOML)
library(logging)
library(stringr)
library(rlang)
library(Seurat)
#install.packages("logging")
#library(rmarkdown)
#library(htmltools)
#library("XML")
#library("bitops")
#library("RCurl")
#if (!require("pacman")) install.packages("pacman")
#pacman::p_load_gh("trinker/textreadr")
wdpath=getwd()
isft=function(x,errdo=T,csize=F){
  if(class(x)!="character"|| is_empty(x)){
    if(errdo==T){
      logerror("dont know x is what or x is empty")
      q()
    }else{
      logwarn("dont know x is what or x is empty")
      return(F)
    }
  }
  if(file.exists(x)){
    loginfo(sprintf("find %s ,ok" ,x))
    if(file.info(x)$size==0){
      if(csize){
        logerror(sprintf("find %s ,but is size is 0" ,x))
        q()
        
      }else{
        loginfo(sprintf("find %s ,but is size is 0" ,x))
      }
    }
    return(T)
  }else{
    if(errdo==T){
      logerror(sprintf("count not find %s ,check it" ,x))
      q()
    }else{
      if(errdo==2){
        return(F)
      }else{
        logwarn(sprintf("count not find %s ,but continue" ,x))
        return(F)
      }
    }
  }
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
isnullorf=function(x){
  if(is.null(x) || x == FALSE){
    x <- FALSE
  }
  return(x)
}

option_list=list(make_option(c("-i","--infile"),  action="store", help='The input file.'))
opt<-parse_args(OptionParser(usage="%prog [options] file\n",option_list=option_list))
#opt$infile="/thinker/nfs4/public/liyq/bin/singlecell/users/liyq/pipline/singlecell/test.ini"
#opt$infile="/thinker/nfs4/public/liyq/bin/singlecell/sgrdata/test.ini"   #shielded by zhangdx
csc<-parseTOML(opt$infile)
#setwd("/thinker/nfs4/public/liyq/bin/singlecell/users/liyq/pipline/singlecell/test/outfile")
wdpath=getwd()
if(csc$outpath$outpath %>% str_detect("^/")){
  loginfo("outfile in %s",csc$outpath$outpath)
  outdir=csc$outpath$outpath}else{
    outdir=paste(wdpath,outdir,sep = "/")
    loginfo("outfile in %s",outdir)
  }

dir.create(outdir,recursive = T)
dir.create(paste(outdir,"shell",sep="/"))
#判断是否同时为真
csc$run$fastp=isnullorf(csc$run$fastp)
csc$run$cellrangle=isnullorf(csc$run$cellrangle)
csc$run$celescope=isnullorf(csc$run$celescope)


if ((csc$run$fastp | csc$run$cellrangle) & csc$run$celescope) {
  logwarn("Either fastp or cellrangle is true, or both, and celescope is true. Exiting...")
  #q()
}else{
  loginfo(sprintf("start %s",ifelse(csc$run$celescope,"celescope","cellrangle")))
  
}
if(csc$run$fastp){
  
  if(!is.null(csc$fastp_cellrange)&&csc$run$fastp){
    dir.create(paste(outdir,"00.fastp",sep="/"))
    
    for(i in 1:length(csc$fastp_cellrange ) ){
      x=csc$fastp_cellrange[[i]]
      if(length(x$R1)>1){
        dir.create(paste(outdir,"00.fastp","temp",sep="/"))
        commanda=sprintf('cat %s  > %s',
                         paste(x$R1,collapse = " "),paste(outdir,"00.fastp","temp",basename(x$R1[1]),sep="/"))
        syscomd(commanda,savef = paste(outdir,"shell","01.QC.sh",sep = "/"))
        commanda=sprintf('cat %s  > %s',
                         paste(x$R2,collapse = " "),paste(outdir,"00.fastp","temp",basename(x$R2[1]),sep="/"))
        syscomd(commanda,savef = paste(outdir,"shell","01.QC.sh",sep = "/"))
        csc$fastp_cellrange[[i]]$R1=paste(outdir,"00.fastp","temp",basename(x$R1[1]),sep="/")
        csc$fastp_cellrange[[i]]$R2=paste(outdir,"00.fastp","temp",basename(x$R2[1]),sep="/")
      }
      dir.create(paste(outdir,"00.fastp",names(csc$fastp_cellrange[i]),sep = "/"))
      setwd(paste(outdir,"00.fastp",names(csc$fastp_cellrange[i]),sep = "/"))
      commanda=sprintf('%s  -h %s_fastp.html -j %s_fastp.json -i %s -I %s -o %s -O %s -w 8 -l %s -n %s -q 10 -b %s -B 0',
                       csc$fastp$fastppath,names(csc$fastp_cellrange[i]),names(csc$fastp_cellrange[i]),
                       csc$fastp_cellrange[[i]]$R1,csc$fastp_cellrange[[i]]$R2,
                       paste(outdir,"00.fastp",names(csc$fastp_cellrange[i]),paste0(names(csc$fastp_cellrange[i]),"_S1_L001_R1_001.fastq.gz"),sep = "/"),
                       paste(outdir,"00.fastp",names(csc$fastp_cellrange[i]),paste0(names(csc$fastp_cellrange[i]),"_S1_L001_R2_001.fastq.gz"),sep = "/"),
                       csc$fastp$longr,csc$fastp$ncode,csc$fastp$longr
      )
      syscomd(commanda,savef = paste(outdir,"shell","01.QC.sh",sep = "/"))
      
      
    }
    
    loginfo("")
  }
}
#cellrangle
if(csc$run$cellrangle){
  dir.create(paste(outdir,"00.cellranger",sep="/"))
  
  for(i in 1:length(csc$fastp_cellrange )){
    
    xsnam=names(csc$fastp_cellrange[i])
    dir.create(paste(outdir,"00.cellranger",xsnam,sep="/"))
    setwd(paste(outdir,"00.cellranger",xsnam,sep="/"))
    if(csc$cellrangle$usedocker){
      commanda=sprintf('docker run --rm --user %s -v %s:%s litd/docker-cellranger bash -c "cd %s && cellranger count --id=%s --transcriptome=%s --fastqs=%s --sample=%s --expect-cells=%s --localcores=%s --localmem=%s "',
                       csc$cellrangle$dockerusr,csc$cellrangle$dir,csc$cellrangle$dir,paste(outdir,"00.cellranger",xsnam,sep="/"),
                       xsnam,csc$cellrangle$ref,paste(outdir,"00.fastp",names(csc$fastp_cellrange[i]),sep = "/"),xsnam,
                       csc$cellrangle$expectcell,csc$cellrangle$localcores,csc$cellrangle$localmem
      )
      syscomd(commanda,savef = paste(outdir,"shell","01.QC.sh",sep = "/"))
      dir.create(paste(outdir,"00.cellranger",xsnam,"outs",sep="/"))
      commanda=sprintf('cd %s && ls %s| grep -v "possorted_genome_bam.bam" |xargs -i ln -s %s/{}',
                       paste(outdir,"00.cellranger",xsnam,"outs",sep="/"),
                       paste(outdir,"00.cellranger",names(csc$fastp_cellrange[i]),names(csc$fastp_cellrange[i]),"outs",sep = "/"),
                       paste(outdir,"00.cellranger",names(csc$fastp_cellrange[i]),names(csc$fastp_cellrange[i]),"outs",sep = "/")
      )
      syscomd(commanda,savef = paste(outdir,"shell","01.QC.sh",sep = "/"))
    }else{
      commanda=sprintf('cd %s && %s count --id=%s --transcriptome=%s --fastqs=%s --sample=%s --expect-cells=%s --localcores=%s --localmem=%s ',
                       paste(outdir,"00.cellranger",xsnam,sep="/"),csc$cellrangle$cellrangpath,
                       xsnam,csc$cellrangle$ref,paste(outdir,"00.fastp",names(csc$fastp_cellrange[i]),sep = "/"),xsnam,
                       csc$cellrangle$expectcell,csc$cellrangle$localcores,csc$cellrangle$localmem
      )
      syscomd(commanda,savef = paste(outdir,"shell","01.QC.sh",sep = "/"))
      dir.create(paste(outdir,"00.cellranger",xsnam,"outs",sep="/"))
      commanda=sprintf('cd %s && ls %s| grep -v "possorted_genome_bam.bam" |xargs -i ln -s %s/{}',
                       paste(outdir,"00.cellranger",xsnam,"outs",sep="/"),
                       paste(outdir,"00.cellranger",names(csc$fastp_cellrange[i]),names(csc$fastp_cellrange[i]),"outs",sep = "/"),
                       paste(outdir,"00.cellranger",names(csc$fastp_cellrange[i]),names(csc$fastp_cellrange[i]),"outs",sep = "/")
      )
      syscomd(commanda,savef = paste(outdir,"shell","01.QC.sh",sep = "/"))
      
      
    }
  }
  
}
if(csc$run$celescope){
  dir.create(paste(outdir,"00.celescope",sep="/"))
  for(i in 1:length(csc$fastp_cellrange ) ){
    x=csc$fastp_cellrange[[i]]
    if(length(x$R1)>1){
      dir.create(paste(outdir,"00.celescope","temp",sep="/"))
      commanda=sprintf('cat %s  > %s',
                       paste(x$R1,collapse = " "),paste(outdir,"00.celescope","temp",basename(x$R1[1]),sep="/"))
      syscomd(commanda,savef = paste(outdir,"shell","01.QC.sh",sep = "/"))
      commanda=sprintf('cat %s  > %s',
                       paste(x$R2,collapse = " "),paste(outdir,"00.celescope","temp",basename(x$R2[1]),sep="/"))
      syscomd(commanda,savef = paste(outdir,"shell","01.QC.sh",sep = "/"))
      csc$fastp_cellrange[[i]]$R1=paste(outdir,"00.celescope","temp",basename(x$R1[1]),sep="/")
      csc$fastp_cellrange[[i]]$R2=paste(outdir,"00.celescope","temp",basename(x$R2[1]),sep="/")
    }
    xsnam=names(csc$fastp_cellrange[i])
    dir.create(paste(outdir,"00.celescope",xsnam,sep="/"))
    setwd(paste(outdir,"00.celescope",xsnam,sep="/"))
    
       commanda=sprintf("bash %s %s %s %s %s %s %i ",
                      csc$celescope$sh,
                     paste(outdir,"00.celescope",xsnam,sep="/"),
                     csc$fastp_cellrange[[i]]$R1,
                     csc$fastp_cellrange[[i]]$R2,
                     names(csc$fastp_cellrange[i]),
                     csc$celescope$ref,
                     csc$celescope$cores)
    #system2
    syscomd(commanda,savef = paste(outdir,"shell","01.QC.sh",sep = "/"))
    
    
}

}
