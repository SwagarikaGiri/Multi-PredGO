#Task is to make the csv file for individual default value

#imports
import numpy as np
import pandas as pd
import pickle
import csv
import sys,os

global RETURN_OBJECT
RETURN_OBJECT=dict()
accession1='P04222'
root='data/'
accession_status_file='AccessionNumber_Structure_StatusFileWithAccessionIndex.pkl'
deepgo_prefix_train='data/train-'
deepgo_prefix_test='data/test-'
pklsuffix='.pkl'
multipred_prefix_train='data/multimodaltrain-'
multipred_prefix_test='data/multimodaltest-'
multipred='data/combined-multimodal-'
accession_status_file_path=root+accession_status_file
df1 = pd.read_pickle(accession_status_file_path)
print df1
raw_input()
accession_object = df1.loc[accession1]
print accession_object
Accesion_No_IndexDict = dict()

def generate_dictionary(index_list):
    counter=0
    for ele in index_list:
        if ele not in Accesion_No_IndexDict:
            Accesion_No_IndexDict[ele]=counter
            counter=counter+1
        else:
            counter=counter+1


#level 2 function
def get_dataframe(baseCode, function,accession_number):
    testData=[]
    if(baseCode=='multipred'):
        PathDataset = multipred+str(function)+pklsuffix
        try:
            df = pd.read_pickle(PathDataset)
            index = df.index.values
            generate_dictionary(index)
            accession_index = Accesion_No_IndexDict[accession_number]
            next_acession_index=accession_index+1
            testData=df.loc[index[accession_index:next_acession_index]]
            return testData
        except:
            return "Sorry the data was not loaded 2"




#level 1 function 
def load_train_test_data(accession_object):
    accession_number=accession_object['accession']
    bp=accession_object['bp']
    cc=accession_object['cc']
    mf=accession_object['mf']
    if(accession_object['status'==True]):
        if(bp):
            testData=get_dataframe('multipred','bp',accession_number)
            if type(testData)!=str:
                RETURN_OBJECT['bp']=testData['gos'].values
               
        if(cc):
            testData=get_dataframe('multipred','cc',accession_number)
            if type(testData)!=str:
                RETURN_OBJECT['cc']=testData['gos'].values
               
        if(mf):
            testData=get_dataframe('multipred','mf',accession_number)
            if type(testData)!=str:
                RETURN_OBJECT['mf']=testData['gos'].values
                
        return RETURN_OBJECT

    elif(accession_object['status'==False]):
        return "Sorry Accession No Cannot be Accepted due to computational limitations2"




RETURN_OBJECT=load_train_test_data(accession_object)
print RETURN_OBJECT
raw_input()
# Now we need to make the csv file for the results that will have the the goterm and the name

filename='data/Goterm_Database.pkl'
df_goterm = pd.read_pickle(filename)



output_file=accession1+'Results.csv'
with open(output_file,'w') as outputcsv_file:
    spamwriter = csv.writer(outputcsv_file,delimiter=',')
    bp=accession_object['bp']
    cc=accession_object['cc']
    mf=accession_object['mf']
    if(bp):
        BP=RETURN_OBJECT['bp']
        for ele in BP[0]:
            col1=[]
            obj= df_goterm.loc[ele]
            col1.append('bp')
            col1.append(obj['id'])
            col1.append(obj['name'])
            spamwriter.writerow(col1)
    if(mf):
        MF=RETURN_OBJECT['mf']
        for ele in MF[0]:
            col1=[]
            obj= df_goterm.loc[ele]
            col1.append('mf')
            col1.append(obj['id'])
            col1.append(obj['name'])
            spamwriter.writerow(col1)
    if(cc):
        CC=RETURN_OBJECT['cc']
        for ele in CC[0]:
            col1=[]
            obj= df_goterm.loc[ele]
            col1.append('cc')
            col1.append(obj['id'])
            col1.append(obj['name'])
            spamwriter.writerow(col1)


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





