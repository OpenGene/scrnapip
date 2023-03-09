
d = {}
with open("KEGG_MsigDB_hsa_entrezID.xls",'r') as f:
	for line in f.readlines():
		line = line.strip()
		lst = line.split("\t")
		if not lst[0] in d:
			d.update({lst[0]:lst[1]})
		else:
			print("opps!")
			print line


w = open("tmp2",'w')

with open("KEGG.MSigDB.txt",'r') as f:
	c = 0
	for line in f.readlines():
		line = line.strip()
		lst = line.split("\t")
		if lst[0] in d:
			w.write(d[lst[0]]+"\t"+line+"\n")
		else:	
			c += 1
#			print("Oh no")
	print "No count :",c
	
w.close()	
