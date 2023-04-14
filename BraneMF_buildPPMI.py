#import packages
import pandas as pd
import wget
from datetime import date
import gzip
import networkx as nx
import scipy

import numpy as np
import os
import argparse
import time

from scipy.io import savemat

#import tensorflow as tf
#import tensorflow_probability as tfp



def parse_args():
    parser = argparse.ArgumentParser(description="performs PPI prediction and evaluation")

    parser.add_argument('--O', nargs='?',
                        help='organism tax ID, e.g. 9606 : homosapeins, 4932 : yeast',type=int)
    
    return parser.parse_args()

def main(args):
    
    import functions as fun

    start = time.time()
    
    isExist = os.path.exists("Network_files")

    if not isExist:
          os.makedirs("Network_files")


    #To keep the information of date for the downloaded files
    today = date.today()
    dat = today.strftime("%d%B%Y")

    #download organisms in string database

    #org_url = "https://stringdb-static.org/download/species.v11.5.txt"
    #org_url_file = "organisms_" + date + ".txt.gz"
    #response = wget.download(org_url, org_url_file)

    #select the taxon ID of organism of your interest
    #for example, 9606 is the tax ID for human 

    org = str(args.O)

    #download protein annotation file

    protein_anno_url = "https://stringdb-static.org/download/protein.enrichment.terms.v11.5/" + org + ".protein.enrichment.terms.v11.5.txt.gz"
    protein_anno_url_file = os.path.join("Network_files", "org_" + org + "_prot_annot_" + dat + ".txt.gz")
    response = wget.download(protein_anno_url, protein_anno_url_file)
    prot_anno_dat = pd.read_csv(protein_anno_url_file, sep = "\t",low_memory=False)

    # To download  metadata
    String_info_url = "https://stringdb-static.org/download/protein.info.v11.5/" + org + ".protein.info.v11.5.txt.gz"
    String_info_url_file = os.path.join("Network_files", "Protein" + org + "_info_" + dat + ".txt.gz")
    response = wget.download(String_info_url, String_info_url_file)
    string_info_dat = pd.read_csv(String_info_url_file, sep= '\t')
    

    #link of files from databases
    # To download networks
    String_network_url = "https://stringdb-static.org/download/protein.links.detailed.v11.5/" + org + ".protein.links.detailed.v11.5.txt.gz"
    String_url_file = os.path.join("Network_files","String_PPI_" + org + "_" + dat + ".txt.gz")
    response = wget.download(String_network_url, String_url_file)
    string_ppi_dat = pd.read_csv(String_url_file, sep= ' ')
    string_ppi_dat = string_ppi_dat[string_ppi_dat['combined_score'] > 899]
    
    string_ppi_dat_net = string_ppi_dat[['protein1','protein2']]
    G_C = nx.from_pandas_edgelist(string_ppi_dat_net, 'protein1', 'protein2')
    node_list = list(G_C.nodes)
    
    #print(node_list)
    #string_info_dat

    node_dict_ = dict(zip(list(string_info_dat['#string_protein_id']),list(range(len(string_info_dat)))))
    
    node_dict = {key: node_dict_[key] for key in node_list}
    all_nodes = list(node_dict.values())


    net_type = ['experimental', 'neighborhood', 'fusion','cooccurence', 'coexpression', 'database']
    
    
    isExist = os.path.exists("PPMI_Matrices")

    if not isExist:
          os.makedirs("PPMI_Matrices")

    #computing RW matrices
    window_size = [1,2,4,6,8,10]
    

    for w in window_size:
        print("\n")
        print("window_size :" + str(w))
        
        M_vec_list = []
        M_list = []

        for typ in range(len(net_type)):

            ppi = string_ppi_dat[['protein1','protein2',net_type[typ]]].copy()
            ppi[net_type[typ]] = ppi[net_type[typ]].div(100).round(2)
            pp = ppi.copy()
            ppi_net = pp.rename(columns = {net_type[typ]: 'weight'}).copy()
            G = nx.Graph()
            G = nx.from_pandas_edgelist(ppi_net, 'protein1', 'protein2',edge_attr=True)
            G = nx.relabel_nodes(G, node_dict)
            A = nx.adjacency_matrix(G, nodelist= all_nodes, weight='weight').copy()
            ppi_file = os.path.join("Network_files", net_type[typ] +  "_" + org  + '_' + dat + ".txt")
            nx.write_edgelist(G,ppi_file, data=True)

            print("------------Computing RW matrix for " + net_type[typ] + "------------------!!!")
            #print("window_size :" + str(w))

            file_name = os.path.join("Network_files", net_type[typ] +  "_" + org  + '_' + dat + ".txt")
     
            M = fun.PPMI_matrix(G,w,1,all_nodes)
            #M_vec = np.reshape(M, -1)
            M_vec_list.append(M.todense())
            M_list.append(M)
            
        net_type_dic = {'experimental' : M_vec_list[0], 'neighborhood': M_vec_list[1], 'fusion' :M_vec_list[2],'cooccurence': M_vec_list[3], 'coexpression': M_vec_list[4], 'database': M_vec_list[5]}
        
        mat_name = os.path.join("PPMI_Matrices",'Org' + org + 'BraneMF_w' + str(w) + '_' + dat + '.mat')
        
        savemat(mat_name, net_type_dic)
            
        
if __name__ == '__main__':
    args = parse_args()
    main(args)
    
