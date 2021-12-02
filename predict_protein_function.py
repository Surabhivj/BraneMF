# -*- coding: utf-8 -*-
"""
Created on Wed May 12 00:43:42 2021

@author: jagtaps
"""
import numpy as np
import scipy.io as sio
import logging
import pandas as pd
import argparse
from code.cv import cross_validation

logger = logging.getLogger(__name__)
theano.config.exception_verbosity='high'

def parse_args():
    parser = argparse.ArgumentParser(description="performs protein function prediction using cross validation ")

    parser.add_argument('--emb', type=str,nargs='?',
                        help='Embedding directory')
    parser.add_argument('--anno', type=str, nargs='?',
                        help='annotation file')
    parser.add_argument('--trials', type=int, nargs='?',
                        help='number of CV trials')
    parser.add_argument('--level', type=str, nargs='?',default = '1',
                        help='levels in annotation file')
    return parser.parse_args()


def write_res(perf, fout):
  fout = open(str(fout)+'.txt', 'w')
  fout.write('aupr[micro] aupr[macro] F_max accuracy\n')         
  for ii in range(0, len(perf['fmax'])):
    fout.write('%0.5f %0.5f %0.5f %0.5f\n' % (perf['pr_micro'][ii], perf['pr_macro'][ii], perf['fmax'][ii], perf['acc'][ii]))
    #fout.write('\n')
  fout.close()
  
anno_file = args.anno
emb_file = args.emb
n_trials = args.trials
level = "level " + args.level

anno = sio.loadmat(anno_file)
anno_level = anno[level]
emb = pd.read_csv(emb_file,delimiter = "\s",index_col= 0, header = None, skiprows= 1)
out_file = emb_file + "." + level
pref_level= cross_validation(emb.values,anno_level,n_trials = n_trials)
write_res(pref_level,out_file)


