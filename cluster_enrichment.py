from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
import scipy.io as sio
from sklearn import metrics
from scipy.stats import ttest_ind
import argparse
import scipy.stats as stats
import glob
import os
import gseapy
import statistics as st




def parse_args():
    parser = argparse.ArgumentParser(description="performs clustering and cluster enrichment/evaluation")

    parser.add_argument('--emb', nargs='?',
                        help='Embedding directory')
    parser.add_argument('--k', nargs='?',
                        help='number of clusters')
    parser.add_argument('--sim', nargs='?',default=20,
                        help='number of simulations')
    parser.add_argument('--genes',nargs='?', default='data\\yeast_string_genes.txt',
                        help='functional annotation file')
  
    return parser.parse_args()

def main(args):
    

    file = args.emb
    #filenames = glob.glob(path + "/*.fmt") 
    genes = pd.read_csv(args.genes,delimiter = ' ',header= None)
    genes = genes[0]
    p_val_all_file = []
    es_all_file = []
    zscore_all_file = []
    results = []
    enriched_file = []
    #for file in range(len(filenames)):
        #print(str(filenames[file]))
    emb = pd.read_csv(file,delimiter = "\s",index_col= 0, header = None, skiprows= 1)
    enriched_file_per = []
    sim = int(args.sim)
    for s in range(sim):
        p_val_all = []
        z_score_all = []
        es_all = []
        k = int(args.k)
        enriched = []
        kmeans = KMeans(n_clusters=k).fit(emb.values)
        silh = metrics.silhouette_score(emb.values,kmeans.labels_)
        d = pd.DataFrame({'proteins' :genes,'cluster': kmeans.labels_})
        f = str(file)+"_k_" + str(k) +"_sim_" + str(s)+"_clustering.res"
        d = d.dropna()
        d.to_csv(f,index = False,line_terminator='\n')
        enrich_scores = []
        #cluster_size = []
        for clus in range(k):
            protein_list = list(d[d.cluster == clus].proteins.values)
            #print(protein_list)
            enrich_score = gseapy.enrichr(gene_list=protein_list, description='pathway',organism='Yeast', gene_sets='GO_Biological_Process_2018', outdir='test')
            enrich_score = enrich_score.res2d
            enrich_score = enrich_score.rename(columns={'Gene_set': 'Gene_set','Term':'Term',
            'Old P-value':'Old_P_value','Old Adjusted P-value':'Old_Adjusted_P_value',
            'Overlap':'Overlap','P-value':'P_value','Adjusted P-value':'Adjusted_P_value',
            'Z-score':'Z_score','Combined Score':'Combined_Score','Genes':'Genes'})
            enrich_score = enrich_score[enrich_score.Adjusted_P_value < 0.05]
            if len(enrich_score.Gene_set) != 0:
                enriched.append(1)
                z_score_all.append(enrich_score.Z_score)
                p_val_all.append(enrich_score.Adjusted_P_value)
                es_all.append(enrich_score.Combined_Score)
                enrich_score["cluster_size"] = [len(protein_list)]*len(enrich_score.Gene_set)
                enrich_scores.append(enrich_score)
                print("Cluster : " +str(clus)+".............Done!")

        flat_pval_list = [item for sublist in p_val_all for item in sublist]
        pval_per_sim = np.sum(np.asarray(flat_pval_list))/len(flat_pval_list)
        flat_zscore_list = [item for sublist in z_score_all for item in sublist]
        zscore_per_sim = np.sum(np.asarray(flat_zscore_list))/len(flat_pval_list)
        flat_es_list = [item for sublist in es_all for item in sublist]
        es_per_sim = np.sum(np.asarray(flat_es_list))/len(flat_pval_list)
        terms_per_sim = np.sum(np.asarray(flat_es_list))
        enriched_file_per.append(np.sum(np.asarray(enriched)))
        p_val_all_file.append(pval_per_sim)
        es_all_file.append(es_per_sim)
        zscore_all_file.append(zscore_per_sim)
        f = open(str(file)+"_k_" + str(k) +"_sim_" + str(s)+"_enrichment.res", 'a')
        for df in enrich_scores:
            df.to_csv(f,index = False,line_terminator='\n')
        f.close()
    dat_res = pd.DataFrame({'method':file,'Total_enriched_clusters':enriched_file_per,
    'simmulation':list(range(20)),'mean_adj_pval':p_val_all_file,'mean_es':es_all_file,
    'Z_score':zscore_all_file})
    fname = str(file)+"clustering_res_"+str(k)+".result"
    dat_res.to_csv(fname,index=False)

if __name__ == '__main__':
    args = parse_args()
    main(args)