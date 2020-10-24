import numpy as np
import math
from sklearn.svm import SVR
from sklearn.ensemble import AdaBoostRegressor, GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import os.path

from util import *

def preprocess_data(data, log_scale, test_size=0.3):
    a = data.to_numpy()
    X = a[:,:-1]
    y = a[:,-1]
    if log_scale: y = np.log10(y + epsilon)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size)
    return X_train, X_test, y_train, y_test

def predict_dataset(dataset, model_name, model, log_scale):
    features_path = '../data/{}_features.csv'.format(dataset)
    pred_path = '../data/pred/{}_{}.csv'.format(dataset, model_name)
    path_data = '../data/{}_data.csv'.format(dataset)
    if not os.path.isfile(features_path):
        print('Creating features')
        create_basic_data(path_data, features_path)
    print('Making predictions', dataset)
    X = load_basic_data(features_path).to_numpy()[:,:-1]
    pred = model.predict(X)
    if log_scale: pred = [math.pow(10, p) for p in pred]
    save_predictions(pred_path, pred)

def get_model(model_name):
    model = None
    if model_name == 'gbr': model = GradientBoostingRegressor()
    elif model_name == 'adaboost': model = AdaBoostRegressor(n_estimators=100)
    elif model_name == 'svr': model = SVR()
    return model

if __name__ == "__main__":
    path_in = '../data/train_data_3.csv'
    path_out = '../data/train_features.csv'

    model_name = 'gbr'
    log_scale = False

    #create_basic_data(path_in, path_out) # uncomment if the data is not yet created
    df = load_basic_data(path_out)
    X_train, X_test, y_train, y_test = preprocess_data(df, log_scale, test_size=0.1)

    model = get_model(model_name)
    model.fit(X_train, y_train)
    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)

    print('train RMSE:', mean_squared_error(y_train, y_train_pred, squared=False))
    print('test RMSE:', mean_squared_error(y_test, y_test_pred, squared=False))

    predict_dataset('rand', model_name, model, log_scale)
    predict_dataset('dimacs', model_name, model, log_scale)
    
    