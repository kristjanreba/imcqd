import numpy as np
import math
from sklearn.svm import SVR
from sklearn.ensemble import AdaBoostRegressor, GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import os.path
from xgboost import XGBRegressor

from util import *
    

def get_model(model_name):
    model = None
    if model_name == 'gbr': model = GradientBoostingRegressor()
    elif model_name == 'adaboost': model = AdaBoostRegressor(n_estimators=100)
    elif model_name == 'svr': model = SVR()
    elif model_name == 'xgb': model = XGBRegressor()
    return model



if __name__ == "__main__":
    datasets_train = ['rand', 'product', 'docking']
    datasets_test = ['rand', 'product', 'docking', 'protein', 'dimacs', 'dense']

    model_names = ['xgb', 'gbr', 'svr']
    log_scale = True
    test_size = 0.1

    # load and process data
    df = load_basic_data(['../datasets/' + d + '_train_features.csv' for d in datasets_train])
    y = load_Tlimits(['../datasets/' + d + '_train_tlimits.csv' for d in datasets_train])
    if log_scale: y = np.log10(np.array(y) + epsilon)
    a = df.to_numpy()
    X = a[:,:-1]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size)

    for model_name in model_names:
        # create model, fit and predict data
        model = get_model(model_name)
        model.fit(X_train, y_train)
        y_train_pred = model.predict(X_train)
        y_test_pred = model.predict(X_test)

        print('model name:', model_name)
        print('train RMSE:', mean_squared_error(y_train, y_train_pred, squared=False))
        print('test RMSE:', mean_squared_error(y_test, y_test_pred, squared=False))

        for dataset in datasets_test:
            features_path = '../datasets/{}_test_features.csv'.format(dataset)
            pred_path = '../pred/{}_{}.csv'.format(dataset, model_name)
            print('Making predictions', dataset)
            X = load_basic_data([features_path]).to_numpy()[:,:-1]
            preds = model.predict(X)
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
    
    
