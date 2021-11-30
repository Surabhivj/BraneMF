import os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, RawTextHelpFormatter
import glob
import pickle
from sklearn.metrics import roc_auc_score, roc_curve,precision_recall_curve,auc
import sys
import os
from sklearn.linear_model import LogisticRegression
from sklearn import metrics, model_selection, pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils import shuffle
from sklearn import preprocessing
import numpy as np
from scipy.spatial import distance

path = "./baselines"
filenames = glob.glob(path + "/*.fmt")
samples_file_path = "./ppi_pred_samples.pkl"

with open(samples_file_path, 'rb') as f:
    temp = pickle.load(f)
    train_samples, train_labels = temp['train']['edges'], temp['train']['labels']
    test_samples, test_labels = temp['test']['edges'], temp['test']['labels']
            
_operators = ["cosine"]#["hadamard", "average", "l1", "l2"]#,"hamming","cosine"]

def extract_feature_vectors_from_embeddings(edges,lab, embeddings, binary_operator):
    features = []
    labels = []
    lab = lab.copy()
    for i in range(len(edges)):
        edge = edges[i]
        labs = lab[i]
        #embeddings = embeddings[~np.all(embeddings == 0, axis=1)]
        if all (k in embeddings for k in (str(edge[0]),str(edge[1]))):
            vec1 = np.asarray(embeddings[str(edge[0])],dtype=np.float64)
            vec2 = np.asarray(embeddings[str(edge[1])],dtype=np.float64)
            value = 0
            if binary_operator == "hadamard":
                value = [vec1[i]*vec2[i] for i in range(len(vec1))]
            if binary_operator == "average":
                value = 0.5 * (vec1 + vec2)
            if binary_operator == "l1":
                value = abs(vec1 - vec2)
            if binary_operator == "l2":
                value = abs(vec1 - vec2)**2
            if binary_operator == "hamming":
                value = 1.0 - distance.hamming(vec1, vec2)
            if binary_operator == "cosine":
                value = distance.cosine(vec1, vec2)

            features.append(value)
            labels.append(labs)

    features = np.asarray(features)
    labels = np.asarray(labels)
        # Reshape the feature vector if it is 1d vector
    if binary_operator in ["hamming", "cosine"]:
        features = features.reshape(-1, 1)

    return features,labels


for filename in filenames:
    embedding_file = filename
    
    print(filename)
    
    embeddings = {}
    with open(embedding_file, 'r') as fin:
        num_of_nodes, dim = fin.readline().strip().split()
        for line in fin.readlines():
            tokens = line.strip().split()
            embeddings[tokens[0]] = [float(v) for v in tokens[1:]]
    print(len(embeddings))
    for i in embeddings.copy() :
        if np.sum(embeddings[i]) == 0:
            embeddings.pop(i)
    print(len(embeddings))

    scores = {op: {'AUPR': [], 'test': []} for op in _operators}

    for op in _operators:
        print("....................training model ................")
        train_features,train_labels = extract_feature_vectors_from_embeddings(edges=train_samples.copy(),lab=train_labels.copy(),
                                                                          embeddings=embeddings,
                                                                          binary_operator=op)

        test_features,test_labels = extract_feature_vectors_from_embeddings(edges=test_samples.copy(),lab=test_labels.copy(),
                                                                         embeddings=embeddings,
                                                                         binary_operator=op)

        clf = LogisticRegression(max_iter=1000)
        clf.fit(train_features, train_labels)
        print("....................testing performance ................")
        train_preds = clf.predict_proba(train_features)[:, 1]
        test_preds = clf.predict_proba(test_features)[:, 1]
        zipped = list(zip(test_labels, test_preds))
        output_file1 =  filename  + "_"+ str(op)+ "_"+ str('_lp_preds')
        np.savetxt(output_file1, zipped, fmt='%i,%i')
        train_roc = roc_auc_score(y_true=train_labels, y_score=train_preds)
        test_roc = roc_auc_score(y_true=test_labels, y_score=test_preds)
        pr, re, thresholds = precision_recall_curve(test_labels, test_preds)
        aupr = auc(re,pr)
        scores[op]['AUPR'].append(aupr)
        #scores[op]['train'].append(train_roc)
        scores[op]['test'].append(test_roc)
        #del train_features
        #del train_labels
        #del test_features
        #del test_labels

    print(scores)
    
    output_file =  filename  + str('_lp_scores')
    with open(output_file, "w") as fp:
        for op in ["cosine"]:#["hadamard", "average", "l1", "l2"]:#, "hamming","cosine"]:    
            fp.write("{} {} {}\n".format(op, scores[op]['AUPR'][0],scores[op]['test'][0]))
    
    
