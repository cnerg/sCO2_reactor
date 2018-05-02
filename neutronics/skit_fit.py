import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import numpy as np
import pandas as pd
import itertools
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

from parse_outputs import filter_data, plot_results
import neutronic_sweeps as ns


data = pd.read_csv("depl_results.csv")
#data = filter_data([('keff', 'great', 1.0)], data)

def linear_regression(predictors, target, ntrain, form_predict=('lin', 'lin')):

    # apply linear or log shift to data
    ops = {'lin' : lambda x: x,
           'log' : lambda x: np.log(x)
          }

    X = ops[form_predict[0]](data[predictors])
    Y = ops[form_predict[1]](data[target])

    X_train = X[:ntrain]
    X_test =  X[ntrain:]

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


def compare_fits():
    """Try various forms for fitting the data, compare results.
    """
    print('Linear - Linear Fit\n-------------------\n')
    linear_regression(['core_r', 'PD', 'enrich'], 'keff', 1500, ('lin', 'lin'))
    print('\nLinear - Log Fit\n----------------\n')
    linear_regression(['core_r', 'PD', 'enrich'], 'keff', 1500, ('lin', 'log'))
    print('\nLog - Linear Fit\n----------------\n')
    linear_regression(['core_r', 'PD', 'enrich'], 'keff', 1500, ('log', 'lin'))
    print('\nLog - Log Fit\n-------------\n')
    linear_regression(['core_r', 'PD', 'enrich'], 'keff', 1500, ('log', 'log'))

compare_fits()
