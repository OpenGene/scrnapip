#awk -F "\t" '{print $7}' test/All_common_snp.anno.vcf.hg19_multianno.xls > test/gene.lst
#sort -u test/gene.lst > test/gene.xls
#head -100 test/gene.xls > test/gene_to_go.xls
/thinker/nfs3/public/wangwf/bin/anaconda3/envs/R3.4.3/bin/Rscript Msigdb2entrezID.R -f KEGG.MSigDB.txt -n KEGG_MsigDB_hsa_entrezID.xls -g 1 -t SYMBOL

## add symbol id 
python KEGG_id_convert.py
awk -F "\t" '{print $2"\t"$4"\t"$3}' KEGG_HAS.txt  > KEGG_HAS_PATHWAY.txt
##KEGG_HAS_PATHWAY.txt used for pathway annotations
