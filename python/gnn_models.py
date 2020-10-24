import numpy as np
import time
import math
from scipy.sparse import csr_matrix

import tensorflow as tf
import tensorflow.keras.backend as K
from tensorflow.keras.layers import Input, Dense, Dropout, BatchNormalization
from tensorflow.keras.losses import MeanSquaredError
from tensorflow.keras.metrics import RootMeanSquaredError
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam

from spektral.layers import GraphConvSkip, GlobalAvgPool, GraphConv, GINConv
from spektral.layers import ops
from spektral.layers.pooling import TopKPool
from spektral.utils.convolution import normalized_adjacency
from spektral.utils.data import batch_iterator, numpy_to_disjoint


class GNN_wrapper():

    def train(self, ):
        pass

    def predict(self, ):
        pass

    def create_model(model_name, load_model):
        pass
        model = None
        F = 1
        n_out = 1
        if model_name == 'gnn': model = create_model_gnn(F, n_out)
        elif model_name == 'gin': model = create_model_gin(F, n_out)
        elif model_name == 'gcn': model = create_model_gcn(F, n_out)
        else:
            print('No model with that name.')
            exit()
        if load_model: model.load_weights(model_name + '.h5')
        return model


class GCN(GNN_wrapper):

    def __init__(self, F, n_out):
        super().__init__()

        X_in = Input(shape=(F, ), name='X_in')
        A_in = Input(shape=(None,), sparse=True)
        I_in = Input(shape=(), name='segment_ids_in', dtype=tf.int32)

        X_1 = GraphConvSkip(256, activation='relu')([X_in, A_in])
        X_1, A_1, I_1 = TopKPool(ratio=0.5)([X_1, A_in, I_in])
        X_2 = GraphConvSkip(256, activation='relu')([X_1, A_1])
        X_2, A_2, I_2 = TopKPool(ratio=0.5)([X_2, A_1, I_1])
        X_3 = GraphConvSkip(256, activation='relu')([X_2, A_2])
        X_3 = GlobalAvgPool()([X_3, I_2])
        output = Dense(n_out)(X_3)
        self.model = Model(inputs=[X_in, A_in, I_in], outputs=output)
    
    def predict(self):
        pass



class GIN(GNN_Model):
    pass


class GAT(GNN_Model):
    pass