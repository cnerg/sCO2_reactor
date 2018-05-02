import numpy as np
import pandas as pd
import itertools
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

from parse_outputs import filter_data, plot_results
import neutronic_sweeps as ns
import plotting as plot

data = pd.read_csv("depl_results.csv")
#data = filter_data([('keff', 'great', 1.0)], data)

def linear_regression(predictors, target, train_frac, form_predict=('lin', 'lin')):

    # apply linear or log shift to data
    ops = {'lin' : lambda x: x,
           'log' : lambda x: np.log(x)
          }

    X = ops[form_predict[0]](data[predictors])
    Y = ops[form_predict[1]](data[target])
    
    ntrain = int(len(X) * train_frac)
    X_train = X[:ntrain]
    X_test = X[ntrain:]

    Y_train = Y[:ntrain]
    Y_test = Y[ntrain:]

    regr = linear_model.LinearRegression()

    regr.fit(X_train, Y_train)

    y_pred = regr.predict(X_test)

    # The coefficients
    print('Coefficients: \n', regr.coef_)
    # The mean squared error
    print("Mean squared error: %.4f"
          % mean_squared_error(Y_test, y_pred))
    # Explained variance score: 1 is perfect prediction
    print('Variance score: %.4f' % r2_score(Y_test, y_pred))

    return Y_test, y_pred, regr

def compare_fits():
    """Try various forms for fitting the data, compare results.
    """
    print('Linear - Linear Fit\n-------------------\n')
    linear_regression(['core_r', 'PD', 'enrich'], 'keff', 0.5, ('lin', 'lin'))
    print('\nLinear - Log Fit\n----------------\n')
    linear_regression(['core_r', 'PD', 'enrich'], 'keff', 0.5, ('lin', 'log'))
    print('\nLog - Linear Fit\n----------------\n')
    linear_regression(['core_r', 'PD', 'enrich'], 'keff', 0.5, ('log', 'lin'))
    print('\nLog - Log Fit\n-------------\n')
    linear_regression(['core_r', 'PD', 'enrich'], 'keff', 0.5, ('log', 'log'))

def test_fit_config(vars, linlog):

    r = []
    mse = []
    train_frac = []

    for i in np.arange(0.01, 0.99, 0.01):
        Y_test, y_pred, regr =  linear_regression(vars,
                'keff', i, linlog)

        r.append(r2_score(Y_test, y_pred))
        mse.append(mean_squared_error(Y_test, y_pred))
        train_frac.append(i)

    return r, mse, train_frac

rs = []
tfrac = []


r, mse, train_frac = test_fit_config(['core_r'], ('log', 'lin'))
rs.append(r)
tfrac.append(train_frac)
r, mse, train_frac = test_fit_config(['core_r', 'PD'], ('log', 'lin'))
rs.append(r)
tfrac.append(train_frac)
r, mse, train_frac = test_fit_config(['core_r', 'PD', 'enrich'], ('log', 'lin'))
rs.append(r)
tfrac.append(train_frac)
r, mse, train_frac = test_fit_config(['core_r', 'PD', 'enrich', 'mass'], ('log', 'lin'))
rs.append(r)
tfrac.append(train_frac)
r, mse, train_frac = test_fit_config(['mass'], ('log', 'lin'))
rs.append(r)
tfrac.append(train_frac)
r, mse, train_frac = test_fit_config(['mass', 'PD'], ('log', 'lin'))
rs.append(r)
tfrac.append(train_frac)
r, mse, train_frac = test_fit_config(['mass', 'PD', 'enrich'], ('log', 'lin'))
rs.append(r)
tfrac.append(train_frac)


plot.plot_results(tfrac,rs, ['R', 'R(PD)', 'R(PD)E', 'R(PD)EM', 'M', 'M(PD)',
    'M(PD)E'])

compare_fits()
