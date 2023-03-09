import mygene
from biothings_client import get_client

mg = get_client('gene')

import pandas as pd

kegg = pd.read_csv("KEGG_hsa_entrezID.xls",sep="\t")


query_df = mg.querymany(kegg.GENE.tolist(),scopes="_id",fields=["symbol"],species="human",as_dataframe=True)

query_df["GENE"] = query_df["_id"]
query_df = query_df.filter(["GENE",'symbol'])

query_df['GENE'] = query_df.GENE.astype('int64')
m = pd.merge(query_df,kegg,on = 'GENE',how = 'inner')

m = m.drop_duplicates(['GENE','symbol','PATHWAY','NAME'],keep='first')

m.to_csv("/haplox/haprs/wangweifeng/shell/clusterProfiler/database/KEGG_HAS.txt",sep = '\t',index=False)
