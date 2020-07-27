import pandas as pd
import numpy as np

for ont in ['cc', 'bp', 'mf']:

    df = pd.read_csv(f'data-cafa/inga-preds-{ont}.csv')
    test_df = pd.read_pickle(f'data-cafa/multimodaltest-{ont}.pkl')

    print(df.shape)
    print(test_df.shape)
    print(test_df.labels.values[0].shape)
    print(1)
    var = input()

    df_2 = pd.read_pickle(f'data-cafa/{ont}.pkl')
    cols = ['accession_no']
    for row in df_2.itertuples():
        cols.append(row.functions)
        if (row.functions) not in df.columns:
            print(row.functions)

    print(2)
    var = input()
    df = df[cols]

    for i, col in enumerate(df.columns):
        if col != cols[i]:
            print(col)

    print(3)
    var = input()

    ele = []

    for row in df.itertuples():
        accession_no = row.accession_no
        annotation = []
        predictions = []
        for idx, item in enumerate(row[2:]):
            if item == 1:
                annotation.append(cols[idx+1])
            predictions.append(item)

        predictions = predictions
        ele.append((accession_no, annotation, predictions))

    new_df = pd.DataFrame(ele, columns=['accession_no', 'gos','predictions'])

    print(4)
    var = input()

    new_df.to_pickle(f'data-cafa/preds_{ont}_inga.pkl')