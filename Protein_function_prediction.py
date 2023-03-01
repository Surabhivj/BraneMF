import numpy as np
#import scipy.io as sio
import logging
import pandas as pd
import argparse
from cv import cross_validation
import os
import ast


#import pickle5 as pickle

import theano

logger = logging.getLogger(__name__)
theano.config.exception_verbosity='high'

def parse_args():
    parser = argparse.ArgumentParser(description="performs protein function prediction using cross validation ")

    parser.add_argument('--emb_dir', type=str,nargs='?',
                        help='Embedding directory')
    parser.add_argument('--anno', type=str, nargs='?',
                        help='annotation file')
    parser.add_argument('--trials', type=int, nargs='?',
                        help='number of CV trials')
    

    return parser.parse_args()


def write_res(perf, fout):
  fout = open(str(fout)+'.txt', 'w')
  fout.write('aupr[micro] aupr[macro] F_max accuracy\n')         
  for ii in range(0, len(perf['fmax'])):
    fout.write('%0.5f %0.5f %0.5f %0.5f\n' % (perf['pr_micro'][ii], perf['pr_macro'][ii], perf['fmax'][ii], perf['acc'][ii]))
    #fout.write('\n')
  fout.close()
  

def main(args):
    anno_file =  pd.read_csv(args.anno,delimiter = "\t",index_col= 0)
    anno_file = anno_file.sort_index()
    emb_dir = args.emb_dir
    n_trials = args.trials
    
    
    emb_files = []
    for file in os.listdir(emb_dir):
        if file.endswith(".emb"):
            emb_files.append(file)
        
    
    for emb_file in emb_files:

        file = os.path.join(emb_dir,emb_file)
        print(file)
        emb = pd.read_csv(file, delimiter = "\s", index_col= 0, header = None, skiprows= 1)
        idx = list(emb.index)
        emb= pd.DataFrame(emb.values,index = idx)

        emb = emb.sort_index()
        emb = emb[emb.index.isin(anno_file.index)]
        anno = anno_file[anno_file.index.isin(emb.index)]
        print(emb.shape)

        out_file = os.path.join(emb_dir,emb_file.split('.')[0] + '_' +  args.anno.split('_')[2] + '_ProtFunc.res')
        pref_level= cross_validation(emb.values,anno.values,n_trials = n_trials)
        write_res(pref_level,out_file)

if __name__ == '__main__':
    args = parse_args()
    main(args)