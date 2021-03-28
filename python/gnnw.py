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

from spektral.layers import GlobalAvgPool, GCNConv, GINConv, GATConv
from spektral.layers import ops
from spektral.layers.pooling import TopKPool
from spektral.utils.convolution import normalized_adjacency
from spektral.data.utils import batch_generator, to_disjoint

from util import *

def create_model(model_name, load_model):
    model = None
    F, n_out = 1, 1
    if model_name == 'gcn': model = create_model_gcn(F, n_out)
    elif model_name == 'gin': model = create_model_gin(F, n_out)
    elif model_name == 'gat': model = create_model_gat(F, n_out)
    else:
        print('No model with that name.')
        exit()
    if load_model: model.load_weights('saved_models/' + model_name + '.h5')
    return model

def create_model_gin(F, n_out):
    X_in = Input(shape=(F, ), name='X_in')
    A = Input(shape=(None,), sparse=True)
    I = Input(shape=(), name='segment_ids_in', dtype=tf.int32)

    X_1 = GINConv(128)([X_in, A])
    X_2 = GINConv(128)([X_1, A])
    X_3 = GlobalAvgPool()([X_2, I])
    X_4 = Dense(128, activation='relu')(X_3)
    output = Dense(n_out)(X_4)
    model = Model(inputs=[X_in, A, I], outputs=output)
    return model

def create_model_gcn(F, n_out):
    X_in = Input(shape=(F, ), name='X_in')
    A = Input(shape=(None,), sparse=True)
    I = Input(shape=(), name='segment_ids_in', dtype=tf.int32)

    X_1 = GCNConv(128, activation='relu')([X_in, A])
    X_2, A_1, I_1 = TopKPool(ratio=0.5)([X_1, A, I])
    X_3 = GCNConv(128, activation='relu')([X_2, A_1])
    X_4 = GlobalAvgPool()([X_3, I_1])
    X_5 = Dense(128, activation='relu')(X_4)
    output = Dense(n_out)(X_5)
    model = Model(inputs=[X_in, A, I], outputs=output)
    return model

def create_model_gat(F, n_out):
    X_in = Input(shape=(F, ), name='X_in')
    A = Input(shape=(None,), sparse=True)
    I = Input(shape=(), name='segment_ids_in', dtype=tf.int32)

    X_1 = GATConv(128, activation='relu')([X_in, A])
    X_2 = GATConv(128, activation='relu')([X_1, A])
    X_3 = GlobalAvgPool()([X_2, I])
    X_4 = Dense(128, activation='relu')(X_3)
    output = Dense(n_out)(X_4)
    model = Model(inputs=[X_in, A, I], outputs=output)
    return model

def make_prediction(A, X, model, log_scale):
    '''
    input -> g is a numpy adjacency matrix
    output -> Tlimit value
    '''
    X_, A_, I_ = to_disjoint([X], [A])
    A_ = ops.sp_matrix_to_sp_tensor(A_)
    pred = model([X_, A_, I_], training=False)
    pred = K.eval(pred[0][0]) # convert from tensor to float
    # Tlimit must be in range 0 <= Tlimit <= 1
    if log_scale:
        if pred < -6: pred = -6
        if pred > 0: pred = 0
        pred = math.pow(10, pred)
    else:
        if pred > 1.0: pred = 1.0
        elif pred < 0.0: pred = 0.0
    return pred

def predict_dataset(dataset, model_name, log_scale=False):
    print('Predicting ' + dataset + ' with ' + model_name)
    path_data = '../datasets/{}_test.csv'.format(dataset)
    path_out = '../pred/{}_{}.csv'.format(dataset, model_name)
    A_list, X_list = load_and_preprocess_test_data_weighted(path_data)
    model = create_model(model_name, True)
    Tlimits = []
    for i in range(len(A_list)):
        A = A_list[i]
        X = X_list[i]
        pred = make_prediction(A, X, model, log_scale)
        Tlimits.append(pred)
    save_predictions(path_out, Tlimits)

def train_eval(model_name):

    def evaluate(A_list, X_list, y_list, ops_list, batch_size):
        batches = batch_generator([X_list, A_list, y_list], batch_size=batch_size, epochs=1)
        output = []
        print('stuck in evaluate')
        for b in batches:
            X, A, I = to_disjoint(*b[:-1])
            A = ops.sp_matrix_to_sp_tensor(A)
            y = b[-1]
            pred = model([X, A, I], training=False)
            outs = [o(y, pred) for o in ops_list]
            #print("{} {}".format(y[0], pred[0]))
            output.append(outs)
        return np.mean(output, 0)

    # Parameters
    load_model = False
    train_model = True
    log_scale = True

    datasets_train = ['rand', 'product', 'docking']
    datasets_test = ['rand', 'product', 'docking', 'protein', 'dimacs', 'dense']
    train_paths = ['../datasets/' + d + '_train.csv' for d in datasets_train]
    paths_tlimits = ['../datasets/' + d + '_train_tlimits.csv' for d in datasets_train]
    #train_paths = ['../datasets/rand_data.csv']
    #paths_tlimits = ['../datasets/bugg.csv']


    learning_rate = 1e-4    # Learning rate for optimizer
    batch_size = 64         # Batch size
    epochs = 100             # Number of training epochs
    es_patience = 50        # Patience fot early stopping
    
    # Create model
    model = create_model(model_name, load_model)
    opt = Adam(lr=learning_rate)
    loss_fn = MeanSquaredError()
    acc_fn = RootMeanSquaredError()

    # Fit model
    if train_model:

        # Load and prepare train data
        print('Loading and preprocessing data...')
        A_train, A_val, A_test,\
        X_train, X_val, X_test,\
        y_train, y_val, y_test,\
        F, n_out = load_and_preprocess_train_data_weighted(train_paths, paths_tlimits, val_size=0.05, test_size=0.01)
        if log_scale:
            y_train = np.log10(y_train + epsilon)
            y_val = np.log10(y_val + epsilon)
            y_test = np.log10(y_test + epsilon)


        #print(X_train[0].dtype)
        #print(A_train[0].dtype)
        #print(type(y_train[0]))

        #@tf.function(
        #    input_signature=(tf.TensorSpec((None, F), dtype=tf.float32),
        #                    tf.SparseTensorSpec((None, None), dtype=tf.float32),
        #                    tf.TensorSpec((None,), dtype=tf.int32),
        #                    tf.TensorSpec((None, n_out), dtype=tf.float32)),
        #    experimental_relax_shapes=True)
        @tf.function(
            experimental_relax_shapes=True
        )
        def train_step(X_, A_, I_, y_):
            with tf.GradientTape() as tape:
                predictions = model([X_, A_, I_], training=True)
                loss = loss_fn(y_, predictions)
                loss += sum(model.losses)
                acc = acc_fn(y_, predictions)
            gradients = tape.gradient(loss, model.trainable_variables)
            opt.apply_gradients(zip(gradients, model.trainable_variables))
            return loss, acc

        
        current_batch = 0
        epoch = 0
        model_loss = 0
        model_acc = 0
        best_val_loss = np.inf
        best_weights = None
        patience = es_patience
        batches_in_epoch = np.ceil(y_train.shape[0] / batch_size)

        print('Fitting model')
        batches = batch_generator([X_train, A_train, y_train],
                                batch_size=batch_size, epochs=epochs)
        for b in batches:
            current_batch += 1

            X_, A_, I_ = to_disjoint(*b[:-1])
            A_ = ops.sp_matrix_to_sp_tensor(A_)
            y_ = b[-1]
            outs = train_step(X_, A_, I_, y_)
            model_loss += outs[0]
            model_acc += outs[1]
            if current_batch == batches_in_epoch:
                epoch += 1
                model_loss /= batches_in_epoch
                model_acc /= batches_in_epoch
                # Compute validation loss and accuracy
                val_loss, val_acc = evaluate(A_val, X_val, y_val, [loss_fn, acc_fn], batch_size=batch_size)
                print('Ep. {} - Loss: {:.6f} - RMSE: {:.4f} - Val loss: {:.6f} - Val RMSE: {:.4f}'
                    .format(epoch, model_loss, model_acc, val_loss, val_acc))
                # Check if loss improved for early stopping
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    patience = es_patience
                    print('New best val_loss {:.6f}'.format(val_loss))
                    best_weights = model.get_weights()
                else:
                    patience -= 1
                    if patience == 0:
                        print('Early stopping (best val_loss: {})'.format(best_val_loss))
                        break
                model_loss = 0
                model_acc = 0
                current_batch = 0

        print('Saving model to file')
        model.set_weights(best_weights)
        model.save_weights('saved_models/' + model_name + '.h5')

        # Evaluate model
        print('Testing model')
        test_loss, test_acc = evaluate(A_test, X_test, y_test, [loss_fn, acc_fn], batch_size=batch_size)
        print('Test loss: {:.4f}. Test RMSE: {:.4f}'.format(test_loss, test_acc))


    for dataset in datasets_test:
        predict_dataset(dataset, model_name, log_scale)


    '''
    # Predict single dimacs graph and measure inference time
    path_dimacs = '../data/dimacs/johnson8-2-4.clq'
    A = load_dimacs(path_dimacs)
    model = create_model('gnn', True)
    print('Making prediction for DIMACS graph')
    X = np.ones((A.shape[0],1))
    start = time.time()
    pred = make_prediction(A, X, model)
    end = time.time()
    print('Predicted Tlimit:', pred)
    print('Time for prediction:', end - start)
    # time for prediction: 0.02
    '''
    

if __name__ == "__main__":
    model_names = ['gcn', 'gat']
    for m in model_names:
        train_eval(m)
