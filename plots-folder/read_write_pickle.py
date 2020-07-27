import pickle
import numpy as np
import pandas as pd
import random


with open('data-cafa/DeepGO-seq-probablitybppreds.pkl', 'rb') as f:
    x = pickle.load(f)
    print(x['gos'])
    gos=x['gos'].values
    accession=x['accession_no'].values
    df = x['predictions'].values
    # print(gos)
    # print(accession)
    # print(len(df[0]))
    # a=np.zeros((589,), dtype=int)
    # b=np.ones((589,), dtype=int)
    # for i in range(0,700):
    #     df[i]=a
    # for i in range(1001,1500):
    #     df[i]=b
    res_df = pd.DataFrame({'accession_no':accession,
		'gos':gos, 'predictions':df})
    print(res_df)
    raw_input()
    res_df.to_pickle('data-cafa/preds-deepgo-seq-bp.pkl')




# with open('data-cafa/preds_bp_inga_formated.pkl', 'rb') as f:
#     x = pickle.load(f)
#     print(x)





    

    
    
    


