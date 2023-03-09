#!/usr/bin/perl -w
#use strict;
use FindBin '$Bin';
use Cwd 'abs_path';

my($all,$top,$cell_percent,$cluster_mean,$odir,$prefix,$help);
use Getopt::Long;
GetOptions(
 "a|all=s"		=>	\$all,
 "t|top=s"		=>	\$top,
 "c|cellpercent=s"	=>	\$cell_percent,
 "m|clustermean=s"	=>	\$cluster_mean,
 "d|directory=s"	=>	\$odir,
 "p|prefix=s"		=>	\$prefix,
 "help"			=>	\$help,
);
sub usage{

print STDERR <<USAGE;

	Options:
		-a|all		<s> : all marker gene
		-t|top		<s> : selected marker gene
		-c|cellpercent	<s> : cell percent
		-m|clustermean	<s> : cluster mean
		-d|directory	<s> : odir
		-p|prefix	<s> : prefix
		-help

	Examples:
		perl $0 -a allmarker_gene.txt -t type2gene.txt -c new.cellpercent.txt -m new.ClusterMean.txt -d ./ -p test

USAGE
}
if(defined $help || !defined $all || !defined $cell_percent || !defined $cluster_mean || !defined $prefix || !defined $top){
	&usage;
	exit;
}
if(!defined $odir){
	$odir="./";
}
$odir=abs_path($odir);
`mkdir -p $odir` if (!-d $odir);
open (IN,"$top") or die ("Can not open $top !\n");
$h=<IN>;
$h=~s///g;
chomp $h;
@h=split/\t/,$h;
for($i=0;$i<=$#h;$i++){
	$key{$h[$i]}=$i;
}
while(<IN>){
	$_=~s///g;
	chomp;
	($gene,$type)=(split/\t/)[$key{'Gene'},$key{'Type'}];
	if(!defined $mark{$type}){
		push @type,$type;
		$mark{$type}=1;
	}
	if(!defined $mark{$gene}){
		push @gene,$gene;
		$mark{$gene}=1;
	}
	$g{$type}.="$gene;";
	$gene{$type}{$gene}=1;
}
close IN;
open OP,">$odir/$prefix.all.point";
print OP "type\tgene\tx\tlog2fc\tpadj\tcolour\tlabel\n";
open (ALL,"$all") or die ("Can not open $all !\n");
$a=<ALL>;
chomp $a;
@a=split/\t/,$a;
for($i=0;$i<=$#a;$i++){
	$col{$a[$i]}=$i;
}
while(<ALL>){
	chomp;
	($type,$gene,$fc,$padj)=(split/\t/)[$col{'cluster'},$col{'Gene'},$col{'avg_log2FC'},$col{'p_val_adj'}];
	#$type=~s/NKcell/Nkcell/g;
	if($type eq "LowQuality"){next;}
	$x=int(rand(500))/500;
	if($padj<0.05){
		$col="#FF6347EE";
	}else{
		$col="#222222EE";
	}
	if($gene{$type}{$gene}){
		$label=1;
		$x=0.5;
		$col="#FF0000FF";
		$line.="$type\t$gene\t$x\t$fc\t$padj\t$col\t$label\n";
	}else{
		$label=0;
		print OP "$type\t$gene\t$x\t$fc\t$padj\t$col\t$label\n";
	}
}
print OP $line;
close ALL;
close OP;
open (CELL,"$cell_percent") or die ("Can not open $cell_percent !\n");
$c=<CELL>;
chomp $c;
@p=split/\t/,$c;
#for($i=0;$i<=$#p;$i++){
#        $key{$p[$i]}=$i;
#}
while(<CELL>){
	chomp;
	@c=split/\t/;
	$sample=$c[1];
	for($i=2;$i<=$#c;$i++){
		$cell_percent{$c[0]}{$sample}{$p[$i]}=$c[$i];
	}
	if(!defined $mark{$sample}){
		push @sample,$sample;
		$mark{$sample}=1;
	}
}
close CELL;

open (MEAN,"$cluster_mean") or die ("Can not open $cluster_mean !\n");
$m=<MEAN>;
chomp $m;
@m=split/\t/,$m;
for($i=0;$i<=$#m;$i++){
	$key{$m[$i]}=$i;
}
$min=100;
$max=0;
while(<MEAN>){
	chomp;
	@c=split/\t/;
	$sample=$c[1];
	for($i=2;$i<=$#c;$i++){
		$cluster_mean{$c[0]}{$sample}{$m[$i]}=$c[$i];
		if($c[$i]>$max){
			$max=$c[$i];
		}
		if($c[$i]<$min){
			$min=$c[$i];
		}
	}
}
$max=int($max)+1;
open O,">$odir/$prefix.all.bubble";
print O "type\tgene\tsample\trank\ty\tcell_percent\tcluster_mean\talpha\theatmap\n";
$y=0;
$sample_count=@sample;
$gap=0.5/($sample_count-1);
for $g(@gene){
	for $type(@type){
		$count=0;
		for $s(@sample){
			$x=0.25+$count*$gap;
			if(!defined $cluster_mean{$type}{$s}{$g}){
#				print "[Warning] $type\t$s\t$g\tNo Value in $cluster_mean\n";
				$cluster_mean{$type}{$s}{$g}=0;
			}
			$alpha=int(int($cluster_mean{$type}{$s}{$g}+1)/$max*255);
			$al=uc(sprintf("%1x",$alpha));
			if(length($al)==1){
				$al="0".$al;
			}
			if(!defined $cell_percent{$type}{$s}{$g}){
#				print "[Warning] $type\t$s\t$g\tNo Value in $cell_percent\n";
				$cell_percent{$type}{$s}{$g}=0;
			}
			if($cell_percent{$type}{$s}{$g}<0.25){
				$cell_percent{$type}{$s}{$g}=0.5;
			}else{
				$cell_percent{$type}{$s}{$g}=1+$cell_percent{$type}{$s}{$g};
			}
			$heat=int($cluster_mean{$type}{$s}{$g}/0.1)+1;
			print O "$type\t$g\t$s\t$x\t$y\t$cell_percent{$type}{$s}{$g}\t$cluster_mean{$type}{$s}{$g}\t$al\t$heat\n";
			$count++;
		}
	}
	$y++;
}
close O;
open M,">$odir/$prefix.mat";
print M "\t".join("\t",@gene)."\n";
for $type(@type){
	for $s(@sample){
		print M "$type\_$s";
		for $g(@gene){
			print M "\t$cluster_mean{$type}{$s}{$g}";
		}
		print M "\n";
	}
}
close M;
open R,">$odir/$prefix.plot.R";
print R << "LINE"
library(circlize)
library(RColorBrewer)

circos.clear()
all<-read.csv("$odir/$prefix.all.bubble", header = T, sep='\\t', check.names=F);
pdf("$odir/$prefix.circos.pdf")
gaps <- rep(1,length(unique(all\$type))-1)
gaps <- append(gaps,20)

circos.par("start.degree"=80,"gap.degree"=gaps)

par(mar = c(1, 1, 1, 1))
all\$type <- factor(all\$type,levels=unique(all\$type))
circos.initialize(factors = all\$type, xlim = c(0, 1))
circos.trackPlotRegion(all\$type, ylim = c(0, 10), track.height = 0.05, bg.border = NA,
	panel.fun = function(x, y) {
		circos.text(0.5, 5, get.cell.meta.data("sector.index"))
	}
)
bg.col<-colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(all\$type)))
circos.track(all\$type, ylim = c(0,1), bg.col=bg.col, track.height = 0.05)
if(length(unique(all\$sample))>9){
	s.col<-colorRampPalette(brewer.pal(9,"Set1"))(length(unique(all\$sample)))
}else{
	s.col<-brewer.pal(9,"Set1")[1:length(unique(all\$sample))]
}
sample.col<-rep(s.col,length(unique(all\$gene))*length(unique(all\$type)))
b.col<-paste(sample.col,all\$alpha,sep="")
circos.track(all\$type, ylim = c(0,length(unique(all\$gene))), track.height = 0.3, bg.border = NA)
circos.trackPoints(all\$type, all\$rank, all\$y, col=b.col, pch = 16, cex=all\$cell_percent)
for(i in 1:length(unique(all\$gene))){
	circos.trackText(unique(all\$type)[i], 0, i-1, unique(all\$gene)[i], cex=0.5)
}
ht.col<-colorRampPalette(brewer.pal(9,"Reds"))(max(unique(all\$heatmap)))
circos.track(all\$type, ylim = c(0,length(unique(all\$gene))), track.height = 0.3, bg.border = NA)
sample.count <- length(unique(all\$sample))
gap <- 1/sample.count
for(j in 1:length(unique(all\$type))){
	if(j==length(unique(all\$type))){
		for(i in 1:length(unique(all\$gene))){
			circos.trackText(unique(all\$type)[j], 1.35, i-0.5, unique(all\$gene)[i], cex=0.5)
		}
	}
	subset<-subset(all,all\$type==unique(all\$type)[j])
	for(i in 1:length(unique(all\$gene))){
		sta <- 0
		for (k in 1:length(unique(all\$sample))){
			end <- sta+gap
			circos.rect(sta,i-1,end,i,col=ht.col[subset\$heatmap[sample.count*(i-1)+k]], border = T, sector.index=unique(all\$type)[j])
			sta <- end
		}
	}
}
point<-read.csv("$odir/$prefix.all.point", header = T, sep='\\t', check.names=F);
circos.track(point\$type, ylim = c(-max(point\$log2fc)/2,0), track.height = 0.18, bg.border = NA)
circos.trackPoints(point\$type, point\$x, -point\$log2fc/2, col=sprintf("%s",point\$colour), pch = 16, cex=0.2)
for(i in 1:length(point\$label)){
	if(point\$label[i]==1){
		circos.trackText(point\$type[i], point\$x[i], -point\$log2fc[i]/2, point\$gene[i], cex=0.25)
		circos.trackPoints(point\$type[i], point\$x[i], -point\$log2fc[i]/2, col="red", pch = 16, cex=0.25)
	}
}
#legend
text(-0.53,1.05,"p.adjust",cex=0.6)
points(-0.58,1,pch=16,col="#FF6347",cex=0.5)
text(-0.48,1,"p.adj<0.05 ",cex=0.5)
points(-0.58,0.95,pch=16,col="#222222",cex=0.5)
text(-0.48,0.95,"p.adj>=0.05",cex=0.5)

text(-0.7,1.05,"sample",cex=0.6)
for(i in 1:length(unique(all\$sample))){
	points(-0.75,1-0.05*i+0.05,pch=16,col=s.col[i],cex=0.5)
	text(-0.68,1-0.05*i+0.05,unique(all\$sample)[i],cex=0.5)
}
text(-0.85,1.05,"pct.exp",cex=0.6)
text(-1,1.05,"avg.exp",cex=0.6)
LINE
;
for($i=0;$i<$max;$i++){
	$j=$i+1;
	$alpha=int($i+1)/$max*255;
	$al=uc(sprintf("%1x",$alpha));
	if(length($al)==1){
		$al="0".$al;
	}
	if($i%2==0){
		$heat=($i/0.1)+1;
		print R "points(-1.05,1-0.05*$i/2,pch=16,col=\"#000000$al\",cex=0.5)\n";
		print R "text(-1,1-0.05*$i/2,\"$i\",cex=0.5)\n";
		print R "points(-1.05,1-0.05*$i/2-0.05*($max-1),pch=15,col=ht.col[$heat],cex=2)\n";
		print R "text(-1,1-0.05*$i/2-0.05*($max-1)+0.025,\"$i\",cex=0.5)\n";
	}
}
if($j%2==0){
	$heat=($j/0.1)+1;
	print R "points(-1.05,1-0.05*$j/2,pch=16,col=\"#000000FF\",cex=0.5)\n";
	print R "text(-1,1-0.05*$j/2,\"$j\",cex=0.5)\n";
	print R "text(-1,1-0.05*$j/2-0.05*($max-1)+0.025,\"$j\",cex=0.5)\n";
}
for($i=0;$i<=4;$i++){
        $size=$i*0.25;
        if($size==0){
                $n=0.5
        }else{
                $n=0.7+$size;
        }
	$size=sprintf("%.2f",$size);
        print R "points(-0.9,1-0.05*$i,pch=16,col=\"#000000\",cex=$n)\n";
        print R "text(-0.85,1-0.05*$i,\"$size\",cex=0.5)\n";
}
print "/usr/bin/Rscript $odir/$prefix.plot.R\n";
`/usr/bin/Rscript $odir/$prefix.plot.R`;
