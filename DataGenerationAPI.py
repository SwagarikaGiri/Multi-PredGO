
#imports
import numpy as np
import pandas as pd
import pickle
import csv
import sys,os
import ForTestMultiPredModelAPI as PredictionModel


#files
root='data/'
deepgo_prefix_train='data/train-'
deepgo_prefix_test='data/test-'
pklsuffix='.pkl'
multipred_prefix_train='data/multimodaltrain-'
multipred_prefix_test='data/multimodaltest-'
multipred='data/combined-multimodal-'
accession_status_file='AccessionNumber_Structure_StatusFileWithAccessionIndex.pkl'
accession_status_file_path=root+accession_status_file
df1 = pd.read_pickle(accession_status_file_path)

Accesion_No_IndexDict = dict()
global RETURN_OBJECT
RETURN_OBJECT=dict()
PAYLOAD=dict()
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
            next_acession_index=accession_index+2
            testData=df.loc[index[accession_index:next_acession_index]]
            return testData
        except:
            return "Sorry the data was not loaded 2"


#level 1 function 
def load_train_test_data(accession_object,ontology):
    RETURN_OBJECT={}
    accession_number=accession_object['accession']
    ontology_flag=accession_object[ontology]
    if(accession_object['status'==True]):
        if (ontology == "bp" and ontology_flag):
            testData=get_dataframe('multipred','bp',accession_number)
            if type(testData)!=str:
                prediction_list = PredictionModel.main('bp',testData,'cpu:0')
                RETURN_OBJECT['bp']=prediction_list
        else:
            if(ontology == "bp"):
                return "Accession no does not have Biological Function"
        if (ontology=="cc" and ontology_flag):


            testData=get_dataframe('multipred','cc',accession_number)
            if type(testData)!=str:
                prediction_list=PredictionModel.main('cc',testData,'cpu:0')
                RETURN_OBJECT['cc']=prediction_list
        else:
            if(ontology=="cc"):
                return "Accession no does not have Cellular Component"
        
        if(ontology=="mf" and ontology_flag):
            testData=get_dataframe('multipred','mf',accession_number)
            if type(testData)!=str:
                prediction_list=PredictionModel.main('mf',testData,'cpu:0')
                RETURN_OBJECT['mf']=prediction_list
        else:
            if(ontology=="mf"):
                return "Accession no does not have Molecular Function"
        return RETURN_OBJECT

    elif(accession_object['status'==False]):
        return "Sorry Accession No Cannot be Accepted due to computational limitations2"


#root function
def analyze_accession_status(accession_number,ontology):
    accession_number=str(accession_number)
    ontology = str(ontology)
    try:
        df1 = pd.read_pickle(accession_status_file_path)
        accession_object = df1.loc[accession_number]
        PAYLOAD = load_train_test_data(accession_object,ontology)
    except:
        PAYLOAD={"Sorry errorenous data"}
    print PAYLOAD
    return PAYLOAD
    


if __name__=='__main__':
    message=analyze_accession_status('P31946')
    
    
