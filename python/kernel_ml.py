import numpy as np
import math
import os.path

from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.svm import SVR
from grakel.kernels import WeisfeilerLehman, VertexHistogram
from grakel import Graph

from util import *


if __name__ == "__main__":
    datasets_train = ['rand', 'product']#, 'docking']
    datasets_test = ['rand', 'product', 'docking', 'protein', 'dimacs', 'dense']

    model_name = 'svr_wl'
    log_scale = True
    test_size = 0.1

    # load and process data
    print('Loading and preprocessing data...')
    train_paths = ['../datasets/' + d + '_train.csv' for d in datasets_train]
    paths_tlimits = ['../datasets/' + d + '_train_tlimits.csv' for d in datasets_train]

    X = load_data(train_paths) # X is a list of graphs
    X = [Graph(g, node_labels={i:'a' for i in range(g.shape[0])}) for g in X]
    y = load_Tlimits(paths_tlimits)
    if log_scale: y = np.log10(np.array(y) + epsilon)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, shuffle=True)
    # create model, fit and predict data
    kernel = WeisfeilerLehman(base_graph_kernel=VertexHistogram, n_jobs=-1)
    K_train = kernel.fit_transform(X_train)
    K_test = kernel.transform(X_test)

    print("Fitting the model ...")
    model = SVR(kernel='precomputed', cache_size=7000, max_iter=100000)
    model.fit(K_train, y_train) # Fit on the train Kernel
    y_train_pred = model.predict(K_train)
    y_test_pred = model.predict(K_test)

    print('model name:', model_name)
    print('train RMSE:', mean_squared_error(y_train, y_train_pred, squared=False))
    print('test RMSE:', mean_squared_error(y_test, y_test_pred, squared=False))

    for dataset in datasets_test:
        graphs_path = '../datasets/{}_test.csv'.format(dataset)
        pred_path = '../pred/{}_{}.csv'.format(dataset, model_name)
        print('Making predictions ...', dataset)
        
        X = load_data([graphs_path]) # X is a list of graphs
        X = [Graph(g, node_labels={i:'a' for i in range(g.shape[0])}) for g in X]
        K_test = kernel.transform(X)
        preds = model.predict(K_test)

        for i in range(len(preds)):
            p = preds[i]
            if log_scale:
                if p < -6: preds[i] = -6
                if p > 0: preds[i] = 0
                preds[i] = math.pow(10, p)
            else:
                if p > 1.0: preds[i] = 1.0
                elif p < 0.0: preds[i] = 0.0
        save_predictions(pred_path, preds)