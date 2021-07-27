
import scipy.io
import scipy.sparse as sparse
from scipy.sparse import csgraph
import numpy as np
import argparse
import logging
import theano
import networkx as nx
from theano import tensor as T

logger = logging.getLogger(__name__)
theano.config.exception_verbosity='high'

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",default=['data\yeast_string_fusion_adjacency.txt',
    'data\yeast_string_neighborhood_adjacency.txt','data\yeast_coocurence_fusion_adjacency.txt',
    'data\yeast_string_coexpression_adjacency.txt','data\yeast_experimental_fusion_adjacency.txt',
    'data\yeast_string_database_adjacency.txt'], nargs='+',help=".txt edgelist input file path")
    parser.add_argument('--matfile-variable-name', default='network',
            help='variable name of adjacency matrix inside a .mat file.')
    parser.add_argument("--output", type=str,default='test.emb',help="embedding output file path")
    parser.add_argument("--dim", default=500, type=int,
            help="dimension of embedding")
    parser.add_argument("--window", default=3,
            type=int, help="context window size")
    parser.set_defaults(large=False)
    return parser.parse_args()

args = parse_arguments()
logging.basicConfig(level=logging.INFO,
format='%(asctime)s %(message)s') # include timestamp

NetList = args.input
dim = args.dim

def RandomWalkMatrix(NetList):
    window = 3
    b = 1
    Nets =[]

    for net in NetList:
        EdgeList = nx.read_weighted_edgelist(net)
        AdjMat = nx.adjacency_matrix(EdgeList)
        Nets.append(AdjMat)
    
    n = Nets[0].shape[0]
    s = (n,n)
    M_sum = np.zeros(s)
    
    for A in Nets:
        n = A.shape[0]
        vol = float(A.sum())
        L, d_rt = csgraph.laplacian(A, normed=True, return_diag=True)
        X = sparse.identity(n) - L
        S = np.zeros_like(X)
        X_power = sparse.identity(n)
        
        for i in range(window):
            logger.info("Compute matrix %d-th power", i+1)
            X_power = X_power.dot(X)
            S += X_power
        
        S *= vol 
        D_rt_inv = sparse.diags(d_rt ** -1)
        M = D_rt_inv.dot(D_rt_inv.dot(S).T)
        M_sum += M.todense()
    
    M_sum *= M_sum / window / b
    m = T.matrix()
    f = theano.function([m], T.log(T.maximum(m, 1)))
    Y = f(M_sum.astype(theano.config.floatX))
    return sparse.csr_matrix(Y)

def LearnEmbeddings(Y,dim):
    u, s, v = sparse.linalg.svds(Y, dim, return_singular_vectors="u")
    Embedding = sparse.diags(np.sqrt(s)).dot(u.T).T
    np.save(args.ouput, Embedding, allow_pickle=False)


