#!/usr/bin/env python

"""
python nn_hierarchical_network.py
"""
import numpy as np
import pandas as pd
import click as ck
import os
import csv
b_size=int(50)
n_epoch=int(10)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
from keras.models import Sequential, Model, load_model
from keras.layers import (
    Dense, Dropout, Activation, Input,
    Flatten, Highway, merge, BatchNormalization,Lambda)
from keras.layers.embeddings import Embedding
from keras.layers.convolutional import (
    Convolution1D, MaxPooling1D)
from keras.optimizers import Adam, RMSprop, Adadelta
from sklearn.metrics import classification_report
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
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.preprocessing import sequence
from keras import backend as K
import sys
from collections import deque
import time
import logging
import tensorflow as tf
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
from scipy.spatial import distance
from multiprocessing import Pool
from keras.layers import Reshape, TimeDistributed
from keras.engine.topology import Layer, InputSpec
from keras import initializations

# from keras_self_attention import SeqSelfAttention

# import os
# os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
# os.environ["CUDA_VISIBLE_DEVICES"] = "0"

import tensorflow as tf
from keras.backend.tensorflow_backend import set_session


""" alternative when GPU is present"""
# config = tf.ConfigProto()
# config.gpu_options.allow_growth = True
# sess = tf.Session(config=config)
# K.set_session(sess)
# print(tf.test.is_gpu_available())

config = tf.ConfigProto(intra_op_parallelism_threads=6, inter_op_parallelism_threads=6, allow_soft_placement=True, log_device_placement=True)
# sess = tf.Session(config=config)
# K.set_session(sess)

global graph
graph = tf.get_default_graph()
sess = tf.Session(graph=graph, config=config)
K.set_session(sess)


logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
sys.setrecursionlimit(100000)

DATA_ROOT = 'data/'
MAXLEN = 1000
REPLEN = 256
ind = 0

def main(function, test_df,device):
    print test_df
    print " check if it is here"
    org = None 
    param = 0 
    filename = 'ResultSequenceStructPPI.txt'
    train = False
    global FUNCTION
    FUNCTION = function
    global GO_ID
    GO_ID = FUNC_DICT[FUNCTION]
    global go
    go = get_gene_ontology('go.obo')
    global ORG
    ORG = org
    func_df = pd.read_pickle(DATA_ROOT + FUNCTION + '.pkl')
    global functions
    functions = func_df['functions'].values
    global func_set
    func_set = set(functions)
    global all_functions
    all_functions = get_go_set(go, GO_ID)
    logging.info('Functions: %s %d' % (FUNCTION, len(functions)))
    global go_indexes
    go_indexes = dict()
    #will be used for my prediction list
    indexes_for_prediction = dict()
    for ind, go_id in enumerate(functions):
        go_indexes[go_id] = ind
        indexes_for_prediction[ind]=go_id
    global node_names
    global FILENAME 
    FILENAME = filename
    global PARAMS
    node_names = set()
    global prediction_list
    with tf.device('/' + device):
        params = {
            'fc_output': 1024,
            'learning_rate': 0.001,
            'embedding_dims': 128,
            'embedding_dropout': 0.2,
            'nb_conv': 1,
            'nb_dense': 1,
            'filter_length': 128,
            'nb_filter': 32,
            'pool_length': 64,
            'stride': 32
        }
        PARAMS=params
        prediction_list=model(params,test_df, is_train=train)
    return prediction_list


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


def get_feature_model(params):
    embedding_dims = params['embedding_dims']
    max_features = 8001
    model = Sequential()
    model.add(Embedding(
        max_features,
        embedding_dims,
        input_length=MAXLEN,
        dropout=params['embedding_dropout']))
    for i in xrange(params['nb_conv']):
        model.add(Convolution1D(
            nb_filter=params['nb_filter'],
            filter_length=params['filter_length'],
            border_mode='valid',
            activation='relu',
            subsample_length=1))
    model.add(MaxPooling1D(
        pool_length=params['pool_length'], stride=params['stride']))
    model.add(Flatten())
    model.summary()
    return model


def merge_outputs(outputs, name):
    if len(outputs) == 1:
        return outputs[0]
    return merge(outputs, mode='concat', name=name, concat_axis=1)


def merge_nets(nets, name):
    if len(nets) == 1:
        return nets[0]
    return merge(nets, mode='sum', name=name)


def get_node_name(go_id, unique=False):
    name = go_id.split(':')[1]
    if not unique:
        return name
    if name not in node_names:
        node_names.add(name)
        return name
    i = 1
    while (name + '_' + str(i)) in node_names:
        i += 1
    name = name + '_' + str(i)
    node_names.add(name)
    return name


def get_function_node(name, inputs):
    output_name = name + '_out'
    net = Dense(256, name=name, activation='relu')(inputs)
    output = Dense(1, name=output_name, activation='sigmoid')(inputs)
    return net, output



def get_layers(inputs):
    q = deque()
    layers = {}
    name = get_node_name(GO_ID)
    layers[GO_ID] = {'net': inputs}
    for node_id in go[GO_ID]['children']:
        if node_id in func_set:
            q.append((node_id, inputs))
    while len(q) > 0:
        node_id, net = q.popleft()
        parent_nets = [inputs]
        # for p_id in get_parents(go, node_id):
        #     if p_id in func_set:
        #         parent_nets.append(layers[p_id]['net'])
        # if len(parent_nets) > 1:
        #     name = get_node_name(node_id) + '_parents'
        #     net = merge(
        #         parent_nets, mode='concat', concat_axis=1, name=name)
        name = get_node_name(node_id)
        net, output = get_function_node(name, inputs)
        if node_id not in layers:
            layers[node_id] = {'net': net, 'output': output}
            for n_id in go[node_id]['children']:
                if n_id in func_set and n_id not in layers:
                    ok = True
                    for p_id in get_parents(go, n_id):
                        if p_id in func_set and p_id not in layers:
                            ok = False
                    if ok:
                        q.append((n_id, net))

    for node_id in functions:
        childs = set(go[node_id]['children']).intersection(func_set)
        if len(childs) > 0:
            outputs = [layers[node_id]['output']]
            for ch_id in childs:
                outputs.append(layers[ch_id]['output'])
            name = get_node_name(node_id) + '_max'
            layers[node_id]['output'] = merge(
                outputs, mode='max', name=name)
    return layers



def multiply(x):
    image,mask = x
    return mask*image

def subtract(x):
    return 1-x

def attention(feature_model,inputs2):
    net = merge(
        [feature_model, inputs2], mode='concat',
        concat_axis=1, name='merged1')
    batch_size = K.shape(net)[0]
    feature_size = K.shape(net)[1]
    alpha =Dense(1,name='alpha',activation='sigmoid')(net)
    feature_model1 = Lambda(multiply)([feature_model,alpha])
    beta = Lambda(subtract)([alpha])
    inputs21=Lambda(multiply)([inputs2,beta])
    return feature_model1,inputs21




def get_model(params):
    logging.info("Building the model")
    inputs = Input(shape=(MAXLEN,), dtype='int32', name='inputs')
    inputs1 = Input(shape=(REPLEN,), dtype='float32', name='input1')
    inputs2 = Input(shape=(REPLEN,), dtype='float32', name='input2')
    feature_model = get_feature_model(params)(inputs)
    net = merge(
        [feature_model, inputs2], mode='concat',
        concat_axis=1, name='merged')
    net = merge(
        [net,inputs1], mode='concat',
        concat_axis=1, name='merged1')
    # net = Dense(1024, activation='relu')(feature_model)
    
    # net = merge(
    #     [feature_model, inputs2], mode='concat',
    #     concat_axis=1, name='merged')
    for i in xrange(params['nb_dense']):
        net = Dense(params['fc_output'], activation='relu')(net)
    layers = get_layers(net)
    output_models = []
    for i in range(len(functions)):
        output_models.append(layers[functions[i]]['output'])
    net =  merge(output_models, mode='concat', concat_axis=1)
    net = Dense(100, activation='relu')(net)
    net = Dense(len(functions), activation='sigmoid')(net)
    model = Model(input=[inputs,inputs1,inputs2], output=net)
    logging.info('Compiling the model')
    optimizer = RMSprop(lr=params['learning_rate'])

    model.compile(
        optimizer=optimizer,
        loss='binary_crossentropy')
    logging.info(
        'Compilation finished')
    return model

def write_file(filename, numpy_array):
    output_file=filename
    with open(output_file,'w') as outputcsv_file:
        spamwriter = csv.writer(outputcsv_file,delimiter=' ')
        for i in range(0,len(numpy_array)):
            col1=[]
            col1=numpy_array[i]
            spamwriter.writerow(col1)


# it is the important function as it creates the output function
def find_the_predicted_go_term(prediction,functions):
    predicted_goterms=[]
    index=0
    for ele in prediction[0]:
        if ele ==1:
            predicted_goterms.append(functions[index])
        index=index+1
    return predicted_goterms
    

def model(params,test_df, batch_size=b_size, nb_epoch=n_epoch, is_train=True):
    # set parameters:
    nb_classes = len(functions)
    start_time = time.time()
    logging.info("Loading Data")
    test, test_df = load_data(test_df)
    test_gos = test_df['gos'].values
    test_data, test_labels = test
    logging.info("Data loaded in %d sec" % (time.time() - start_time))
    logging.info("Test data size: %d" % len(test_data[0]))
    model_path = (DATA_ROOT + 'models/model_' + FUNCTION + '.h5') 
    checkpointer = ModelCheckpoint(
        filepath=model_path,
        verbose=1, save_best_only=True)
    earlystopper = EarlyStopping(monitor='val_loss', patience=10, verbose=1)
    logging.info('Starting training the model')
    batch_size_new=2
    test_generator = DataGenerator(batch_size_new, nb_classes)
    test_generator.fit(test_data, test_labels)
    logging.info('Loading best model')
    pred={}
    start_time = time.time()
    with graph.as_default():
        model = load_model(model_path)
        logging.info('Loading time: %d' % (time.time() - start_time))
        start_time = time.time()
        preds = model.predict_generator(
            test_generator, val_samples=len(test_data[0]))
    running_time = time.time() - start_time
    logging.info('Running time: %d %d' % (running_time, len(test_data[0])))
    logging.info('Computing performance')
    f, p, r, t, preds_max = compute_performance(preds, test_labels, test_gos)
    roc_auc = compute_roc(preds, test_labels)
    mcc = compute_mcc(preds_max, test_labels)
    logging.info('Fmax measure: \t %f %f %f %f' % (f, p, r, t))
    logging.info('ROC AUC: \t %f ' % (roc_auc, ))
    logging.info('MCC: \t %f ' % (mcc, ))
    print('f :%.3f & p: %.3f & r: %.3f & roc_auc: %.3f & mcc: %.3f' % (
        f, p, r, roc_auc, mcc))
    proteins = test_df['proteins']
    predictions = list()
    for i in xrange(preds_max.shape[0]):
        predictions.append(preds_max[i])
    counter2=0
    for ele in test_labels[0]:
        if ele==1:
            counter2=counter2+1
    print counter2
    counter=0
    for ele in predictions[0]:
        if ele ==1:
            counter=counter+1
    print counter
    prediction_list=find_the_predicted_go_term(predictions,functions)
    print prediction_list
    return prediction_list


def load_prot_ipro():
    proteins = list()
    ipros = list()
    with open(DATA_ROOT + 'swissprot_ipro.tab') as f:
        for line in f:
            it = line.strip().split('\t')
            if len(it) != 3:
                continue
            prot = it[1]
            iprs = set(it[2].split(';'))
            proteins.append(prot)
            ipros.append(iprs)
    return pd.DataFrame({'proteins': proteins, 'ipros': ipros})


def performanc_by_interpro():
    pred_df = pd.read_pickle('test'+FUNCTION+'preds.txt')
    ipro_df = load_prot_ipro()
    df = pred_df.merge(ipro_df, on='proteins', how='left')
    ipro = get_ipro()

    def reshape(values):
        values = np.hstack(values).reshape(
            len(values), len(values[0]))
        return values

    for ipro_id in ipro:
        if len(ipro[ipro_id]['parents']) > 0:
            continue
        labels = list()
        predictions = list()
        gos = list()
        for i, row in df.iterrows():
            if not isinstance(row['ipros'], set):
                continue
            if ipro_id in row['ipros']:
                labels.append(row['labels'])
                predictions.append(row['predictions'])
                gos.append(row['gos'])
        pr = 0
        rc = 0
        total = 0
        p_total = 0
        for i in xrange(len(labels)):
            tp = np.sum(labels[i] * predictions[i])
            fp = np.sum(predictions[i]) - tp
            fn = np.sum(labels[i]) - tp
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
                pr += precision
                rc += recall
        if total > 0 and p_total > 0:
            rc /= total
            pr /= p_total
            if pr + rc > 0:
                f = 2 * pr * rc / (pr + rc)
                print('%s\t%d\t%f\t%f\t%f' % (
                    ipro_id, len(labels), f, pr, rc))


def function_centric_performance(functions, preds, labels):
    preds = np.round(preds, 2)
    for i in xrange(len(functions)):
        f_max = 0
        p_max = 0
        r_max = 0
        x = list()
        y = list()
        for t in xrange(1, 100):
            threshold = t / 100.0
            predictions = (preds[i, :] > threshold).astype(np.int32)
            tp = np.sum(predictions * labels[i, :])
            fp = np.sum(predictions) - tp
            fn = np.sum(labels[i, :]) - tp
            sn = tp / (1.0 * np.sum(labels[i, :]))
            sp = np.sum((predictions ^ 1) * (labels[i, :] ^ 1))
            sp /= 1.0 * np.sum(labels[i, :] ^ 1)
            fpr = 1 - sp
            x.append(fpr)
            y.append(sn)
            precision = tp / (1.0 * (tp + fp))
            recall = tp / (1.0 * (tp + fn))
            f = 2 * precision * recall / (precision + recall)
            if f_max < f:
                f_max = f
                p_max = precision
                r_max = recall
        num_prots = np.sum(labels[i, :])
        roc_auc = auc(x, y)
        print('%s %f %f %f %d %f' % (
            functions[i], f_max, p_max, r_max, num_prots, roc_auc))


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


def get_gos(pred):
    mdist = 1.0
    mgos = None
    for i in xrange(len(labels_gos)):
        labels, gos = labels_gos[i]
        dist = distance.cosine(pred, labels)
        if mdist > dist:
            mdist = dist
            mgos = gos
    return mgos


def compute_similarity_performance(train_df, test_df, preds):
    logging.info("Computing similarity performance")
    logging.info("Training data size %d" % len(train_df))
    train_labels = train_df['labels'].values
    train_gos = train_df['gos'].values
    global labels_gos
    labels_gos = zip(train_labels, train_gos)
    p = Pool(64)
    pred_gos = p.map(get_gos, preds)
    total = 0
    p = 0.0
    r = 0.0
    f = 0.0
    test_gos = test_df['gos'].values
    for gos, tgos in zip(pred_gos, test_gos):
        preds = set()
        test = set()
        for go_id in gos:
            if go_id in all_functions:
                preds |= get_anchestors(go, go_id)
        for go_id in tgos:
            if go_id in all_functions:
                test |= get_anchestors(go, go_id)
        tp = len(preds.intersection(test))
        fp = len(preds - test)
        fn = len(test - preds)
        if tp == 0 and fp == 0 and fn == 0:
            continue
        total += 1
        if tp != 0:
            precision = tp / (1.0 * (tp + fp))
            recall = tp / (1.0 * (tp + fn))
            p += precision
            r += recall
            f += 2 * precision * recall / (precision + recall)
    return f / total, p / total, r / total







if __name__ == '__main__':
    main()
