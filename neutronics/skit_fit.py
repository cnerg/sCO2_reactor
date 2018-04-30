import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

from parse_outputs import filter_data, plot_results
import neutronic_sweeps as ns


data = pd.read_csv("depl_results.csv")
data = filter_data([('keff', 'great', 1.0)], data)

X = data[['core_r', 'enrich', 'PD', 'cool_r', 'power']]
Y = data['keff']

midpoint = int(len(data) / 2)
X_train = np.log(X[:midpoint])
X_test = np.log(X[midpoint:])

Y_train = Y[:midpoint]
Y_test = Y[midpoint:]

regr = linear_model.LinearRegression()

regr.fit(X_train, Y_train)

y_pred = regr.predict(X_test)

# The coefficients
print('Coefficients: \n', regr.coef_)
# The mean squared error
print("Mean squared error: %.4f"
      % mean_squared_error(Y_test, y_pred))
# Explained variance score: 1 is perfect prediction
print('Variance score: %.3f' % r2_score(Y_test, y_pred))
print(regr.get_params())

def plot_results(ind, dep, colorplot=[]):
    """Generate Plots
    """
    # plot
    fig = plt.figure()
    if len(colorplot) > 0:
        plt.scatter(ind, dep, c=colorplot, s=6,
                cmap=plt.cm.get_cmap('plasma', len(set(colorplot))))
        plt.colorbar()
    else:
        plt.scatter(ind, dep, s=6)
    
    plt.show()

    return plt

err = np.abs(y_pred - Y_test)
plot_results(data['mass'][midpoint:], err, np.exp(X_test['enrich']))
