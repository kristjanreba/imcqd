import numpy as np
import networkx as nx
import pandas as pd
import re
from pathlib import Path
import tensorflow as tf
from scipy.sparse import csr_matrix
from sklearn.model_selection import train_test_split
from spektral.utils.convolution import normalized_adjacency

epsilon = 1e-6

def load_data(path):
    '''
    Read file that contains graphs and
    best Tlimit value for each graph.

    return: list of adjacency matrices and a list of Tlimit values
    '''
    graphs = []
    Tlimits = []
    f = open(path, 'r')
    lines = f.readlines()
    counter = 0
    n = int(lines[counter])
    counter += 1
    for i in range(n):
        n_vertices, n_edges = lines[counter].split(" ")
        n_vertices, n_edges = int(n_vertices), int(n_edges)
        counter += 1
        g = np.zeros((n_vertices, n_vertices))
        for j in range(n_edges):
            u, v = lines[counter].split(" ")
            u, v = int(u), int(v)
            g[u,v] = 1
            counter += 1
        graphs.append(g)
        Tlimits.append(float(lines[counter]))
        counter += 1
    Tlimits = np.array(Tlimits).reshape(-1,1)
    return graphs, Tlimits

def save_data(path, graphs, Tlimits):
    '''
    graphs -> list of graps in adjacency matrix format
    Tlimits -> list of float values for each graph

    function creates a save file of the format that is readable by c++ part
    of this project.
    '''
    f = open(path, 'w')
    n = len(graphs)
    f.write(str(n) + "\n")
    for i in range(n):
        g = graphs[i]
        n_vertices = g.shape[0]
        n_edges = np.sum(g)
        f.write("{} {}\n".format(n_vertices, int(n_edges)))
        for j in range(n_vertices):
            for k in range(j):
                if g[j,k]: f.write("{} {}\n".format(j, k))
        f.write(str(Tlimits[i]) + "\n")

def split_data(X, y, validation_size=0.2, test_size=0.1):
    train_size = 1. - (validation_size + test_size)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=1-train_size, shuffle=True)
    X_val, X_test, y_val, y_test = train_test_split(X_test, y_test, test_size=test_size/(test_size+validation_size), shuffle=False)
    return X_train, X_val, X_test, np.array(y_train), np.array(y_val), np.array(y_test)

def get_ones_as_feature_vectors(A_list):
    return[np.ones((a.shape[0],1)) for a in A_list]

def cast_list_to_float32(X):
    return [tf.dtypes.cast(x, tf.float32) for x in X]

def load_and_preprocess_train_data(path, val_size=0.1, test_size=0.1):
    A, y = load_data(path) # A is a list of graphs, y is a list of float values
    y += epsilon
    A = cast_list_to_float32(A)
    y = cast_list_to_float32(y)
    A_train, A_val, A_test, y_train, y_val, y_test = split_data(A, y, validation_size=val_size, test_size=test_size)    
    X_train = get_ones_as_feature_vectors(A_train)
    X_val = get_ones_as_feature_vectors(A_val)
    X_test = get_ones_as_feature_vectors(A_test)
    print("Train size:", len(A_train))
    print("Val_size:", len(A_val))
    print("Test_size:", len(A_test))

    X_train = cast_list_to_float32(X_train)
    X_val = cast_list_to_float32(X_val)
    X_test = cast_list_to_float32(X_test)

    # Preprocessing
    A_train = [normalized_adjacency(csr_matrix(a)) for a in A_train]
    A_val = [normalized_adjacency(csr_matrix(a)) for a in A_val]
    A_test = [normalized_adjacency(csr_matrix(a)) for a in A_test]

    F = X_train[0].shape[-1] # number of feauteres for each node (we dont have any features so we set a single feature to 1)
    n_out = 1 # regression requires 1 output value (Tlimit)

    return A_train, A_val, A_test, X_train, X_val, X_test, y_train, y_val, y_test, F, n_out

def load_and_preprocess_test_data(path):
    A_list, _ = load_data(path) # A is a list of graphs, y is a list of float values
    A_list = cast_list_to_float32(A_list)
    X_list = get_ones_as_feature_vectors(A_list)
    X_list = cast_list_to_float32(X_list)
    A_list = [normalized_adjacency(csr_matrix(a)) for a in A_list]
    return A_list, X_list

def load_dimacs(path):
    f = open(path, 'r')
    lines = f.readlines()
    g = None
    for l in lines:
        if l[0] == 'c': continue
        elif l[0] == 'p':
            _, _, n_vertices, _ = l.split()
            n_vertices = int(n_vertices)
            g = np.zeros((n_vertices, n_vertices))
        elif l[0] == 'e':
            _, u, v = l.split(" ")
            u, v = int(u)-1, int(v)-1 # substract one to make name of vertices start at 0
            g[u,v] = 1
        else:
            print('Unknown line in dimacs graph.')
            exit()
    return g

def save_graph(path, g, Tlimit):
    save_data(path, [g], [Tlimit])

def save_clq_to_csv(path_in, path_out):
    directory_in_str = path_in
    save_file = path_out
    pathlist = Path(directory_in_str).rglob('*.clq')
    pathlist = sorted(pathlist)
    graphs = []
    Tlimits = []
    i = 0
    for path in pathlist:
        path_in_str = str(path)
        print(i, path_in_str)
        g = load_dimacs(path_in_str)
        graphs.append(g)
        Tlimits.append(0)
        i += 1
    save_data(save_file, graphs, Tlimits)

def save_predictions(path, pred):
    f = open(path, 'w')
    n = len(pred)
    f.write(str(n) + "\n")
    for p in pred:
        f.write(str(p) + "\n")

def load_results(path):
    steps = []
    times = []
    f = open(path, 'r')
    lines = f.readlines()
    n = int(lines[0])
    for l in lines[1:]:
        s, t = l.split()
        s, t = int(s), float(t)
        steps.append(s)
        times.append(t)
    return steps, times

def load_dist(path):
    Tlimits = []
    times = []
    f = open(path, 'r')
    lines = f.readlines()
    n = int(lines[0])
    for l in lines[1:]:
        s, t = l.split()
        s, t = float(s), float(t)
        Tlimits.append(s)
        times.append(t)
    return Tlimits, times


# Classical ML

def get_features(A):
    g = nx.from_numpy_matrix(A)
    d = {}
    d['num_nodes'] = nx.number_of_nodes(g)
    d['num_edges'] = nx.number_of_edges(g)
    d['density'] = nx.density(g)
    d['average_degree'] = np.sum(d for n, d in g.degree()) / float(g.number_of_nodes())
    d['max_degree'] = np.max([d for n, d in g.degree()])
    d['min_degree'] = np.min([d for n, d in g.degree()])
    d['avg_clustering'] = nx.average_clustering(g)
    return d

def create_basic_data(path_in, path_out):
    print('Creating features ... ', end='')
    graphs, Tlimits = load_data(path_in)
    list_of_dict = []
    c = 0
    for g in graphs:
        d = get_features(g)
        d['Tlimit'] = Tlimits[c][0]
        list_of_dict.append(d)
        c += 1
    df = pd.DataFrame(list_of_dict)
    df.to_csv(path_out)
    print('Done')

def load_basic_data(path):
    df = pd.read_csv(path, index_col=0)
    return df


if __name__ == "__main__":
    pass
