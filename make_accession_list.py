import pickle
import numpy as np
import pandas as pd
import csv



root='data/'
accession_status_file='AccessionNumber_Structure_StatusFileWithAccessionIndex.pkl'
accession_status_file_path=root+accession_status_file

df1 = pd.read_pickle(accession_status_file_path)
print(df1.loc['Q8L607'])
df=df1[df1['status']==True]

outfile='List_Accession_No.txt'


with open(outfile,'w') as outputcsv_file:
    spamwriter = csv.writer(outputcsv_file,delimiter=',')
    for ele in df['accession']:
    	print ele
    	spamwriter.writerow([ele])
   