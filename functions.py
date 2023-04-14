import numpy as np
import scipy.sparse as sparse
from sklearn.preprocessing import normalize
from scipy.sparse import csgraph

import logging
import scipy.io as sio
import networkx as nx
#import theano
#from theano import tensor as T

logger = logging.getLogger(__name__)



def PPMI_matrix(G, window, b, all_nodes):

    A = nx.adjacency_matrix(G, nodelist= all_nodes, weight='weight').copy()
    degree = np.array([val for (node, val) in G.degree()])
    d_rt = np.diag(np.diag(degree)).astype(float) 
    L = nx.normalized_laplacian_matrix(G)
    n = A.shape[0]
    vol = float(A.sum())
    X = sparse.identity(n) - L
    S = np.zeros_like(X)
    X_power = sparse.identity(n)
    for i in range(window):
        X_power = X_power.dot(X)
        S += X_power
    S *= vol / window / b
    D_rt_inv = sparse.diags(d_rt ** -1)
    M = D_rt_inv.dot(D_rt_inv.dot(S).T)
    Y = normalize(M, axis =1, norm ='l1')
    return sparse.csr_matrix(Y)
   

    
    
