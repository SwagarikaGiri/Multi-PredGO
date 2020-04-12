


import numpy as np
import pandas as pd


from utils import (
    get_gene_ontology,
    get_go_set,
    get_anchestors,
    get_parents,
    DataGenerator,
    FUNC_DICT,
    MyCheckpoint,
    save_model_weights,
    load_model_weights,
    get_ipro)
from keras.preprocessing import sequence
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
import logging

function='mf'
DATA_ROOT='data/'
MAXLEN = 1000



global FUNCTION
FUNCTION = function
global GO_ID
GO_ID = FUNC_DICT[FUNCTION]
global go
go = get_gene_ontology('go.obo')
func_df = pd.read_pickle(DATA_ROOT + FUNCTION + '.pkl')
global functions
functions = func_df['functions'].values
global func_set
func_set = set(functions)
global all_functions
all_functions = get_go_set(go, GO_ID)


def compute_roc(preds, labels):
    # Compute ROC curve and ROC area for each class
    fpr, tpr, _ = roc_curve(labels.flatten(), preds.flatten())
    roc_auc = auc(fpr, tpr)
    return roc_auc

def compute_mcc(preds, labels):
    # Compute ROC curve and ROC area for each class
    mcc = matthews_corrcoef(labels.flatten(), preds.flatten())
    return mcc

def compute_performance(preds, labels, gos):
    preds = np.round(preds, 2)
    f_max = 0
    p_max = 0
    r_max = 0
    t_max = 0
    for t in xrange(1, 100):
        threshold = t / 100.0
        predictions = (preds > threshold).astype(np.int32)
        total = 0
        f = 0.0
        p = 0.0
        r = 0.0
        p_total = 0
        for i in range(labels.shape[0]):
            tp = np.sum(predictions[i, :] * labels[i, :])
            fp = np.sum(predictions[i, :]) - tp
            fn = np.sum(labels[i, :]) - tp
            all_gos = set()
            for go_id in gos[i]:
                if go_id in all_functions:
                    all_gos |= get_anchestors(go, go_id)
            all_gos.discard(GO_ID)
            all_gos -= func_set
            fn += len(all_gos)
            if tp == 0 and fp == 0 and fn == 0:
                continue
            total += 1
            if tp != 0:
                p_total += 1
                precision = tp / (1.0 * (tp + fp))
                recall = tp / (1.0 * (tp + fn))
                p += precision
                r += recall
        if p_total == 0:
            continue
        r /= total
        p /= p_total
        if p + r > 0:
            f = 2 * p * r / (p + r)
            if f_max < f:
                f_max = f
                p_max = p
                r_max = r
                t_max = threshold
                predictions_max = predictions
    return f_max, p_max, r_max, t_max, predictions_max




def load_data(test_df):
    test_df=test_df
    def reshape(values):
        values = np.hstack(values).reshape(
            len(values), len(values[0]))
        return values

    def normalize_minmax(values):
        mn = np.min(values)
        mx = np.max(values)
        if mx - mn != 0.0:
            return (values - mn) / (mx - mn)
        return values - mn



    def get_values(data_frame):
        labels = reshape(data_frame['labels'].values)
        ngrams = sequence.pad_sequences(
            data_frame['ngrams'].values, maxlen=MAXLEN)
        ngrams = reshape(ngrams)
        rep = reshape(data_frame['struct_feature'].values)
        emb = reshape(data_frame['embeddings'].values)
        data = (ngrams, rep ,emb)
        return data, labels
    test = get_values(test_df)
    return test,test_df


def reshape(values):
    values = np.hstack(values).reshape(
        len(values), len(values[0]))
    return values

# load test_df
test_df = pd.read_pickle(DATA_ROOT + 'multimodaltest' + '-' + FUNCTION + '.pkl')
test_gos = test_df['gos'].values
test, test_df = load_data(test_df)
test_data, test_labels = test
pred_df = pd.read_pickle('pred_'+FUNCTION+'.pkl')
preds= reshape(pred_df['annotations'].values)
f, p, r, t, preds_max = compute_performance(preds, test_labels, test_gos)
roc_auc = compute_roc(preds, test_labels)
mcc = compute_mcc(preds_max, test_labels)
logging.info('Fmax measure: \t %f %f %f %f' % (f, p, r, t))
logging.info('ROC AUC: \t %f ' % (roc_auc, ))
logging.info('MCC: \t %f ' % (mcc, ))
print('f :%.3f & p: %.3f & r: %.3f & roc_auc: %.3f & mcc: %.3f' % (
    f, p, r, roc_auc, mcc))



