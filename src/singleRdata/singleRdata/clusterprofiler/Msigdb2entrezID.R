suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
library("getopt")
command=matrix(c(
"file","f",1,"character",
"name","n",1,"character",
"gene.col","g",1,"character",
"idtype","t",1,"character",
"help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if (!is.null(args$help) || is.null(args$file) || is.null(args$name)) {
    cat(paste(getopt(command, usage = T), "\n"))
    q()
}
mydata<-read.table(args$file,header=TRUE,sep="\t")
a=as.numeric(args$gene.col)
print(a)
colnames(mydata)[a]<-"AccID"
bi<-bitr(as.character(mydata$AccID), fromType = args$idtype,toType = "ENTREZID",OrgDb = "org.Hs.eg.db", drop = TRUE)
print(head(bi,3))
write.table(bi,args$name,sep="\t",quote =FALSE,row.names=FALSE)
