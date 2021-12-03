# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 10:41:56 2021

@author: jagtaps
"""

import pandas as pd
import numpy as np
import networkx as nx
import statistics
import pickle
import argparse
import random
import ast


def parse_args():
    parser = argparse.ArgumentParser(description="sample the test and train data with negative datasets")

    parser.add_argument('--new', nargs='?',
                        help='updated PPIs of 2021')
    parser.add_argument('--old', nargs='?',
                        help='old PPIs of 2015 dataset')
    parser.add_argument('--genes',nargs='?', default='data\\yeast_string_genes.txt',
                        help='gene_list')
  
    return parser.parse_args()


def main(args):
    #'data\\yeast_string_refnet_2021.txt'
    #data\\yeast_string_genes.txt
    #data\\joinedFile.txt
    refnet = nx.read_weighted_edgelist(args.new)
    genes = pd.read_csv(args.genes,delimiter = ' ',header= None)
    genes = genes[0]

    gene_dic = {genes[i]:i for i in range(0, len(genes))}
    refnet = nx.relabel_nodes(refnet, gene_dic)
    refnet = nx.Graph(refnet.subgraph(list(range(0, len(genes)))))
    dist_nodes = set(list(range(0, len(genes)))) - set(refnet.nodes)
    refnet.add_nodes_from(dist_nodes)

    input_nets = nx.read_weighted_edgelist(args.old)
    genes2 = list(input_nets.nodes)
    genes2_fmt = [int(i) for i in genes2]
    genes2_fmt = [x - 1 for x in genes2_fmt]
    gene_dic2 = dict(zip(genes2, genes2_fmt))
    input_nets = nx.relabel_nodes(input_nets, gene_dic2)

    refnet_edges =[(a,b) for a, b, attrs in refnet.edges(data=True) if attrs["weight"] >= 500]
    G = nx.Graph()
    G.add_nodes_from(refnet)
    G.add_edges_from(refnet_edges)

    test_graph = nx.difference(G,input_nets)
    train_graph = nx.intersection(G,input_nets)


    pos_test_edges =list(test_graph.edges)
    pos_train_edges =list(train_graph.edges)

    non_edges_refnet = list(nx.non_edges(refnet))
    non_edges_input = list(nx.non_edges(input_nets))

    non_edges = list(set(non_edges_refnet) & set(non_edges_input))

    np.random.shuffle(non_edges)

    neg_train_edges = non_edges[:len(pos_train_edges)]
    neg_test_edges = non_edges[len(pos_train_edges):(len(pos_train_edges)+ len(pos_test_edges))]

    # Prepare the positive and negative samples for training set
    train_samples = pos_train_edges + neg_train_edges
    train_labels = [1 for _ in pos_train_edges] + [0 for _ in neg_train_edges]
    # Prepare the positive and negative samples for test set
    test_samples = pos_test_edges + neg_test_edges
    test_labels = [1 for _ in pos_test_edges] + [0 for _ in neg_test_edges]

    #save_file_path = "train.pkl"
    
    train = {'edges':train_samples, 'labels': train_labels}
    test = {'edges':test_samples, 'labels': test_labels}

    pickle.dump(train, open("./data/train.pkl", "wb"))
    pickle.dump(test, open("./data/test.pkl", "wb"))
    '''
    with open(save_file_path, 'wb') as f:
        pickle.dump({'train': {'edges':train_samples, 'labels': train_labels }},f)
        pickle.dump({'test': {'edges':test_samples, 'labels': test_labels}}, f)
        f.close()
    '''

if __name__ == '__main__':
    args = parse_args()
    main(args)