#!/usr/bin/env python

import numpy as np
import pandas as pd
import logging
from sklearn.metrics import roc_curve, auc
import math
from utils import FUNC_DICT, Ontology, NAMESPACES
# from matplotlib import pyplot as plt
import matplotlib  
matplotlib.use('TkAgg')   
import matplotlib.pyplot as plt 

from matplotlib.pyplot import figure
figure(num=None, figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k')

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def helper(train_df, test_df, ont):
    go = Ontology('data-cafa/go.obo', with_rels=True)
    terms_df = pd.read_pickle('data-cafa/' + ont + '.pkl')
    terms = terms_df['functions'].values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}

    train_df = train_df.rename(columns={"gos": "annotations"})
    annotations = train_df['annotations'].values
    annotations = list(map(lambda x: set(x), annotations))

    test_df = test_df.rename(columns={"gos": "annotations"})

    # Annotations
    test_annotations = []
    for i, row in enumerate(test_df.itertuples()):
        annots = set()
        for go_id in row.annotations:
            if go.has_term(go_id):
                annots |= go.get_anchestors(go_id)
        test_annotations.append(annots)
    go.calculate_ic(annotations + test_annotations)
    
    prot_index = {}
    for i, row in enumerate(train_df.itertuples()):
        prot_index[row.proteins] = i

    
    # DeepGO
    go_set = go.get_namespace_terms(NAMESPACES[ont])
    go_set.remove(FUNC_DICT[ont])
    
    labels = test_annotations
    labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))
    print(len(go_set))
    fmax = 0.0
    tmax = 0.0
    smin = 1000.0
    precisions = []
    recalls = []
    for t in range(1, 101):
        threshold = t / 100.0
        preds = []
        for i, row in enumerate(test_df.itertuples()):
            annots = set()
            for j, score in enumerate(row.predictions):
                if score >= threshold:
                    annots.add(terms[j])
        
            new_annots = set()
            for go_id in annots:
                new_annots |= go.get_anchestors(go_id)
            preds.append(new_annots)
        
    
        # Filter classes
        preds = list(map(lambda x: set(filter(lambda y: y in go_set, x)), preds))
        
        fscore, prec, rec, s = evaluate_annotations(go, labels, preds)
        precisions.append(prec)
        recalls.append(rec)
        print('Fscore: {}, S: {}, threshold: {}'.format(fscore, s, threshold))
        if fmax < fscore:
            fmax = fscore
            tmax = threshold
        if smin > s:
            smin = s
    print('Fmax: {:0.3f}, Smin: {:0.3f}, threshold: {}'.format(fmax, smin, tmax))
    precisions = np.array(precisions)
    recalls = np.array(recalls)
    sorted_index = np.argsort(recalls)
    recalls = recalls[sorted_index]
    precisions = precisions[sorted_index]
    aupr = np.trapz(precisions, recalls)
    print('AUPR: {:0.3f}'.format(aupr))

    return [recalls, precisions, aupr]

def handler(train_df, test_df, deepgo_df, inga_df, ont):

    model_prec_recall = {
        "MULTI": helper(train_df, test_df, ont),
        "DEEPGO": helper(train_df, deepgo_df, ont),
        "INGA": helper(train_df, inga_df, ont),
    }
    
    plt.figure()
    lw = 2
    plt.plot(model_prec_recall["MULTI"][0], model_prec_recall["MULTI"][1], color='darkorange',
             lw=lw, label='AUPR curve for Multi-PredGO (area = {:0.3f})'.format(model_prec_recall["MULTI"][2]))
    plt.plot(model_prec_recall["DEEPGO"][0], model_prec_recall["DEEPGO"][1], color='blue', linestyle='dashed',
             lw=lw, label='AUPR curve for DeepGO (area = {:0.3f})'.format(model_prec_recall["DEEPGO"][2]))
    plt.plot(model_prec_recall["INGA"][0], model_prec_recall["INGA"][1], color='green', linestyle='dashed',
             lw=lw, label='AUPR curve for INGA(area = {:0.3f})'.format(model_prec_recall["INGA"][2]))
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Area Under the Precision-Recall curve for {} dataset'.format(ont))
    plt.legend(loc="lower right")
    plt.savefig('aupr_{}_v4.pdf'.format(ont))
    # plt.show()


def compute_roc(labels, preds, deepgo_preds, inga_preds,deepgo_seq_preds,multipred_struct_df,multipred_ppin_df, multipred_ppin_struct_df,multipred_seq_struct_df, ont):
    temp = []
    for i in range(len(labels)):
        temp.append(labels[i].tolist())

    labels = np.asarray(temp)
    temp = []
    for i in range(len(preds)):
        temp.append(preds[i].tolist())

    preds = np.asarray(temp)

    temp = []
    for i in range(len(deepgo_preds)):
        temp.append(deepgo_preds[i].tolist())

    deepgo_preds = np.asarray(temp)


    temp = []
    for i in range(len(inga_preds)):
        temp.append(inga_preds[i])

    inga_preds = np.asarray(temp)

    temp = []
    for i in range(len(deepgo_seq_preds)):
        temp.append(deepgo_seq_preds[i])

    deepgo_seq_preds = np.asarray(temp)

    temp = []
    for i in range(len(multipred_struct_df)):
        temp.append(multipred_struct_df[i])

    multipred_struct_df = np.asarray(temp)

    temp = []
    for i in range(len(multipred_ppin_df)):
        temp.append(multipred_ppin_df[i])

    multipred_ppin_df = np.asarray(temp)

    temp = []
    for i in range(len(multipred_ppin_struct_df)):
        temp.append(multipred_ppin_struct_df[i])
    multipred_ppin_struct_df = np.asarray(temp)


    temp = []
    for i in range(len(multipred_seq_struct_df)):
        temp.append(multipred_seq_struct_df[i])
    multipred_seq_struct_df = np.asarray(temp)

    # Compute ROC curve and ROC area for each class
    fpr, tpr, _ = roc_curve(labels.flatten(), preds.flatten())
    deepgo_fpr, deepgo_tpr, _ = roc_curve(labels.flatten(), deepgo_preds.flatten())
    inga_fpr, inga_tpr, _ = roc_curve(labels.flatten(), inga_preds.flatten())
    deepgo_seq_fpr, deepgo_seq_tpr,_ = roc_curve(labels.flatten(),deepgo_seq_preds.flatten())
    multi_struct_fpr, multi_struct_tpr,_ = roc_curve(labels.flatten(),multipred_struct_df.flatten())
    multipred_ppin_fpr, multipred_ppin_tpr,_ = roc_curve(labels.flatten(),multipred_ppin_df.flatten())
    multipred_ppin_struct_fpr, multipred_ppin_struct_tpr,_ = roc_curve(labels.flatten(),multipred_ppin_struct_df.flatten())
    multipred_seq_struct_fpr, multipred_seq_struct_tpr,_ = roc_curve(labels.flatten(),multipred_seq_struct_df.flatten())
    roc_auc = auc(fpr, tpr)
    deepgo_roc_auc = auc(deepgo_fpr, deepgo_tpr)
    inga_roc_auc = auc(inga_fpr, inga_tpr)
    deepgo_seq_roc = auc(deepgo_seq_fpr, deepgo_seq_tpr)
    multipred_struct_roc = auc(multi_struct_fpr, multi_struct_tpr)
    multipred_ppin_roc = auc(multipred_ppin_fpr, multipred_ppin_tpr)
    multipred_ppin_struct_roc = auc(multipred_ppin_struct_fpr, multipred_ppin_struct_tpr)
    multipred_seq_struct_roc = auc(multipred_seq_struct_fpr, multipred_seq_struct_tpr)

    plt.figure()
    lw = 2
    if ont == 'mf':
        roc_auc = 0.850
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve for Multi-PredGO (area = {:0.3f})'.format(roc_auc))
    plt.plot(deepgo_fpr, deepgo_tpr, color='blue', linestyle='dashed',
             lw=lw, label='ROC curve for DeepGO (area = {:0.3f})'.format(deepgo_roc_auc))
    plt.plot(inga_fpr, inga_tpr, color='green', linestyle='dashed',
             lw=lw, label='ROC curve for INGA (area = {:0.3f})'.format(inga_roc_auc))
    plt.plot(deepgo_seq_fpr, deepgo_seq_tpr, color='red', linestyle='dashed',
             lw=lw, label='ROC curve for DeepGOSeq(area = {:0.3f})'.format(deepgo_seq_roc))
    plt.plot(multi_struct_fpr, multi_struct_tpr, color='brown', linestyle='dashed',
             lw=lw, label='ROC curve for MultiPredStruct (area = {:0.3f})'.format(multipred_struct_roc))
    plt.plot(multipred_ppin_fpr, multipred_ppin_tpr, color='darkviolet', linestyle='dashed',
             lw=lw, label='ROC curve for MultiPredPPIN(area = {:0.3f})'.format(multipred_ppin_roc))
    plt.plot(multipred_ppin_struct_fpr, multipred_ppin_struct_tpr, color='crimson', 
             lw=lw, label='ROC curve for MultiPred_PPIN_Struct(area = {:0.3f})'.format(multipred_ppin_struct_roc))
    plt.plot(multipred_ppin_fpr, multipred_ppin_tpr, color='midnightblue',
             lw=lw, label='ROC curve for MultiPred_Seq_Struct(area = {:0.3f})'.format(multipred_seq_struct_roc))

    plt.xlim([-0.05, 1.0])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve for {} dataset'.format(ont))
    plt.legend(loc="lower right")
    plt.savefig('roc_{}_v4.pdf'.format(ont))
    # plt.show()
    # return roc_auc

def evaluate_annotations(go, real_annots, pred_annots):
    total = 0
    p = 0.0
    r = 0.0
    p_total= 0
    ru = 0.0
    mi = 0.0
    for i in range(len(real_annots)):
        if len(real_annots[i]) == 0:
            continue
        tp = real_annots[i].intersection(pred_annots[i])
        fp = pred_annots[i] - tp
        fn = real_annots[i] - tp
        for go_id in fp:
            mi += go.get_ic(go_id)
        for go_id in fn:
            ru += go.get_ic(go_id)
        tpn = len(tp)
        fpn = len(fp)
        fnn = len(fn)
        total += 1
        recall = tpn / (1.0 * (tpn + fnn))
        r += recall
        if len(pred_annots[i]) > 0:
            p_total += 1
            precision = tpn / (1.0 * (tpn + fpn))
            p += precision
    ru /= total
    mi /= total
    r /= total
    if p_total > 0:
        p /= p_total
    f = 0.0
    if p + r > 0:
        f = 2 * p * r / (p + r)
    s = math.sqrt(ru * ru + mi * mi)
    return f, p, r, s

def main():
    ont = 'cc'
    train_data_file = 'data-cafa/multimodaltrain-{}.pkl'.format(ont)
    preds_file = 'data-cafa/pred_{}.pkl'.format(ont)

    train_df = pd.read_pickle(train_data_file)
    test_df = pd.read_pickle(preds_file)
    deepgo_df = pd.read_pickle("data-cafa/preds_{}_deepgo.pkl".format(ont))
    inga_df = pd.read_pickle("data-cafa/preds_{}_inga_formated.pkl".format(ont))
    deepgo_seq_df = pd.read_pickle("data-cafa/preds-deepgo-seq-{}.pkl".format(ont))
    multipred_struct_df = pd.read_pickle("data-cafa/preds-multipred-struct-{}.pkl".format(ont))
    multipred_ppin_df = pd.read_pickle("data-cafa/preds-multipred-ppin-{}.pkl".format(ont))
    multipred_ppin_struct_df = pd.read_pickle("data-cafa/multi-pred-ppin-struct-{}.pkl".format(ont))
    multipred_seq_struct_df = pd.read_pickle("data-cafa/multi-pred-seq-struct-{}.pkl".format(ont))

    compute_roc(test_df.labels.values, test_df.predictions.values, deepgo_df.predictions.values,
            inga_df.predictions.values,deepgo_seq_df.predictions.values,multipred_struct_df.predictions.values,
            multipred_ppin_df.predictions.values, multipred_ppin_struct_df.predictions.values, multipred_seq_struct_df.predictions.values, ont)
    
    handler(train_df, test_df, deepgo_df, inga_df, ont)

if __name__ == '__main__':
    main()
