import pandas as pd
import numpy as np
import csv
import re


gofile = open("data/go.txt", "r")
# print gofile.read()
ans = pd.DataFrame(columns = ['id', 'name','index'])
list_goterm=[]
list_goterm_name=[]
for line in gofile:
    txt = "id:"
    list_=line.split()
    name="name:"
    if not list_:
        continue
    else:
        if(str(list_[0])==txt):
            list_goterm.append(list_[1])
        if(str(list_[0])==name):
            gene_name=""
            for i in range(1,len(list_)):
                gene_name=gene_name+" "+list_[i]
            list_goterm_name.append(gene_name)
                




df_goterm = pd.DataFrame(list_goterm)
df_name = pd.DataFrame(list_goterm_name)
for i in range(0,len(list_goterm)):
    ans = ans.append({'id': list_goterm[i], 'name': list_goterm_name[i],'index':list_goterm[i]},ignore_index=True)
print ans.set_index('index')
# ans = pd.DataFrame(columns = ['id', 'name'])
# res_df = pd.DataFrame({'id':df_goterm,'name':df_name})
# print res_df
ans.to_pickle('data/Goterm_Database.pkl')
raw_input()