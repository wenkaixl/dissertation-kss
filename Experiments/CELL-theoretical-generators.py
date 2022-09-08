import numpy as np
import matplotlib.pyplot as plt
import time
import pickle
import os

os.chdir("..")
os.chdir("./Data")

import networkx as nx
import igraph
import graphkernels as gk
import pandas as pd

from cell.utils import link_prediction_performance
from cell.cell import Cell, EdgeOverlapCriterion, LinkPredictionCriterion
from cell.graph_statistics import compute_graph_statistics

import torch
import scipy
import scipy.sparse as sp
from scipy.sparse import load_npz

import warnings
warnings.filterwarnings("ignore")

### load training graphs from explicit null model
train_graph_list = [None] * 500

## for E2S-model
#df_e2s = pd.read_csv("CELL_originals_E2S.csv", header = None)
#df_e2s = df_e2s.iloc[1:100,]

for i in range(500):
    ## for Barabasi-Albert
    #g = igraph.Graph.Barabasi(n=20, m=1, outpref=False, directed=False, power=1, zero_appeal=1, implementation='psumtree', start_from=None)             
    #g = g.to_networkx()
    #train_graph = sp.csc_matrix(nx.to_numpy_matrix(g))
    #train_graph_list[i] = train_graph

    ## for Geometric Random Graph
    #g = nx.random_geometric_graph(20, 0.3, dim=2, pos=None, p=2, seed=None)
    #train_graph = sp.csc_matrix(nx.to_numpy_matrix(g))
    #train_graph_list[i] = train_graph

    ## for E2S-model
    # create ERGM with R package, save into csv, then load to train CELL
    vec = df_e2s.iloc[i,]
    vec = np.array(vec)
    vec = vec.astype(int)
    mat = vec.reshape((20,20))
    train_graph = sp.csc_matrix(mat)
    train_graph_list[i] = train_graph

### generate samples from CELL
adjacency_overview = pd.DataFrame(index=range(500),columns=range(20*20))
original_overview = pd.DataFrame(index=range(500),columns=range(20*20))

for i in range(500):
    ### create model
    model_cell = Cell(A=train_graph_list[i],
                 H = 10,
                 callbacks=[EdgeOverlapCriterion(invoke_every=10, edge_overlap_limit=.5)])

    ### train model 
    model_cell.train(steps=200,
                optimizer_fn=torch.optim.Adam,
                optimizer_args={'lr': 0.1,
                                'weight_decay': 1e-7})

    generated_graph = model_cell.sample_graph()
    adjacency_matrix = generated_graph.toarray()
    adjacency_vector = adjacency_matrix.flatten()
    adjacency_overview.iloc[i,range(20*20)] = adjacency_vector
        
    original_adjacency_matrix = train_graph_list[i].toarray()
    original_adjacency_vector = original_adjacency_matrix.flatten()
    original_overview.iloc[i,range(20*20)] = original_adjacency_vector
    
    if i%10 == 0:
        print(i)

#adjacency_overview.to_csv("CELL_graphs_from_geometric_no_torus.csv", index = False, header = False)
#original_overview.to_csv("CELL_originals_geometric_no_torus.csv", index = False, header = False) 