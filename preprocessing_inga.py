
import pandas as pd
import numpy as np
import csv


query = 'bp'
ont = 'P'
root='data/'
outfile='inga-predictions-'+query+'.csv'
# Get different terms from ontology pickle file
go_term_file=root+query+'.pkl'
terms_df = pd.read_pickle(go_term_file)
terms = terms_df['functions'].values.flatten()


# Get list of accessions passed in testing
test_df = pd.read_pickle(root+'multimodaltest-'+query+'.pkl')
accessions = test_df['accessions'].values.flatten()

inga_pred_file=root+'inga_pred_file.csv'

# Read prediction file from INGA
columns = ['accessions', 'annotations', 'go_ont', 'INGA_score', 'dist_to_ont', 'leaf', 'go_definition']
df = pd.read_csv(inga_pred_file, sep='\t', names=columns)
df = df[['accessions', 'annotations', 'go_ont']]
# Filter dataframe according to its ontology
df = df[df['go_ont'] == ont]
df=df[['accessions','annotations']]


inga_pred_dict=dict()
for i, row in enumerate(df.itertuples()):
    key=(row.accessions,row.annotations)
    if key in inga_pred_dict:
        pass
    else:
        inga_pred_dict[key]=1


# generate the csv file



# with open(outfile,'w') as outputcsv_file:
#     spamwriter = csv.writer(outputcsv_file,delimiter=',')
#     col1=['accession_no']
#     col1=col1+terms.tolist()
#     spamwriter.writerow(col1)
#     print len(col1)
#     raw_input()
#     for accession in accessions:
#         col1=[]
#         col1.append(accession)
#         for goterm in terms:
#             if (accession,goterm) in inga_pred_dict:
#                 col1.append('1')
#             else:
#                 col1.append('0')
#         print len(col1)
#         spamwriter.writerow(col1)

# ans = pd.DataFrame(columns = ['accessions', 'annotations'])
# with open(outfile,'w') as outputcsv_file:
#     spamwriter = csv.writer(outputcsv_file,delimiter=',')
#     col1=['accession_no']
#     col1=col1+terms.tolist()
#     spamwriter.writerow(col1)
#     print len(col1)
#     raw_input()
#     for accession in accessions:
#         col1=[]
#         col1.append(accession)
#         for goterm in terms:
#             if (accession,goterm) in inga_pred_dict:
#                 col1.append('1')
#             else:
#                 col1.append('0')
#         print len(col1)
#         spamwriter.writerow(col1)

terms_dict = {v: i for i, v in enumerate(terms)}
nb_classes = len(terms)

ans = pd.DataFrame(columns = ['accessions', 'annotations'])
# Get Matrix
res = dict()
# for i, row in enumerate(df.itertuples()):
#     if row.annotations in terms:
#         if row.accessions in res.keys():
#             res[row.accessions][terms_dict[row.annotations]] = 1
#         else:
#             res[row.accessions] = np.zeros((nb_classes), dtype=np.int32)
#             res[row.accessions][terms_dict[row.annotations]] = 1

for accession in accessions:
    for term in terms:
        if accession in res.keys():
            if (accession,term) in inga_pred_dict:
                res[accession][terms_dict[term]]=1.00
            else:
                res[accession][terms_dict[term]]=0.00


        else:
            res[accession] = np.zeros((nb_classes), dtype=np.int32)
            if (accession,term) in inga_pred_dict:
                res[accession][terms_dict[term]]=float(1.00)
            else:
                res[accession][terms_dict[term]]=float(0.00)


for key in res.keys():
    ans = ans.append({'accessions': key, 'annotations': res[key]}, ignore_index = True)

print(len(ans))
print ans
raw_input()

ans.to_pickle('pred_'+query+'.pkl')

