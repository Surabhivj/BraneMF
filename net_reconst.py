# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 15:40:56 2021
@author: jagtaps
project: BraneMF
"""
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances,paired_distances
import numpy as np
from sklearn.metrics import precision_recall_curve,recall_score,matthews_corrcoef
from sklearn.metrics import roc_curve,average_precision_score,auc
import scipy.io as sio
from sklearn.metrics import pairwise_distances
import ast
import networkx as nx
import statistics
from sklearn import preprocessing
from scipy.spatial import distance_matrix
import argparse

import glob
import os


def parse_args():
    parser = argparse.ArgumentParser(description="performs clustering and cluster enrichment/evaluation")

    parser.add_argument('--emb',type=str, nargs='?',
                        help='Embedding directory')
    parser.add_argument('--refnet', type=str, nargs='?',default='./data/yeast_string_refnet_2021.txt',
                        help='reference network of yeast')
    parser.add_argument('--genes',nargs='?', type=str, default='data\\yeast_string_genes.txt',
                        help='gene list')  
    return parser.parse_args()


def net_recons(emb,refnet):
    
    nodes = list(refnet.nodes)
    nodes.sort()
    a1 = nx.adjacency_matrix(refnet,nodelist = nodes)
    a1 = a1.todense()
    np.fill_diagonal(a1, 1)
    refnet_adjmat_dat = pd.DataFrame(a1)

    Mat_cosine = pairwise_distances(emb,metric='cosine')       
    min_max_scaler = preprocessing.MinMaxScaler()
    MatNormDF_cosine = min_max_scaler.fit_transform(Mat_cosine)
    MatNormDF_cosine = 1-MatNormDF_cosine
    MatNormDF_cosine = pd.DataFrame(MatNormDF_cosine)
    MatNormDF_melt_cosine = MatNormDF_cosine.stack().reset_index()
    probs_cosine  = MatNormDF_melt_cosine[0]

    refnet_df_melt = refnet_adjmat_dat.stack().reset_index()
    y = refnet_df_melt[0]
    y = np.where(y>0.0, 1, y)

    d = {'y': y,'prob': probs_cosine}
    dat = pd.DataFrame(d)
    
    return dat

def mcc_t(dat):
    y = dat['y']
    prob_s = dat['prob']
    thres = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    mcc = []
    for i in thres:
        prob = np.array([1 if item > i else 0 for item in prob_s])
        m = matthews_corrcoef(y,prob)
        mcc.append(m)
    p_d = {'thres': thres, 'dat_mcc': mcc}
    pk_d = pd.DataFrame(p_d)
    return(pk_d)


def pr_at_k(dat):
    k = [1,10,100,1000,10000,100000,1000000]
    dat_dist = dat.sort_values(ascending=False, by='prob')
    s_dist = []
    for i in k:
        ss_dist = ((sum(dat_dist.y[:i].values))/i)
        s_dist.append(ss_dist)
 
    p = {'k': k, 'precision_dist': s_dist}
    pk = pd.DataFrame(p)
    return pk

def main(args):
    
    #reading arguments
    emb = args.emb
    refnet = args.refnet
    genes = args.genes
    refnet = nx.read_weighted_edgelist(refnet)
    emb = pd.read_csv(emb,delimiter = "\s",index_col= 0, header = None, skiprows= 1)
    genes = pd.read_csv(args.genes,delimiter = ' ',header= None)
    
    #preprocessing
    
    genes = genes[0]
    gene_dic = {genes[i]:i for i in range(0, len(genes))}
    refnet = nx.relabel_nodes(refnet, gene_dic)
    refnet = nx.Graph(refnet.subgraph(list(range(0, len(genes)))))
    dist_nodes = set(list(range(0, len(genes)))) - set(refnet.nodes)
    refnet.add_nodes_from(dist_nodes)
    print("\n")    
    print("\n")
    #performning reconstruction
    branemf_dat = net_recons(emb.values,refnet)
    print("...........Computing Precision@k...........\n")

    #computing Precision@k
    branemf_dat_pk = pr_at_k(branemf_dat)
    print("\n")
    print(branemf_dat_pk.to_markdown())
    print("\n")
    print("...........Computing Mathews correlation coefficient...........\n")

    #computing Mathews correlation coefficient
    branemf_dat_mcc = mcc_t(branemf_dat)
    print(branemf_dat_mcc.to_markdown())
    print("|........................|")

if __name__ == '__main__':
    args = parse_args()
    main(args)
