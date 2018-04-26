import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

from parse_outputs import filter_data, plot_results
import neutronic_sweeps as ns


data = pd.read_csv("depl_results.csv")
data = filter_data([('keff', 'great', 1.0)], data)


X = data[['core_r', 'enrich', 'PD']]
Y = data['keff']

X_train = np.log(X[:int(len(data)/2)])
X_test = np.log(X[int(len(data)/2):])

Y_train = np.log(Y[:int(len(data)/2)])
Y_test = np.log(Y[int(len(data)/2):])

regr = linear_model.LinearRegression()

regr.fit(X_train, Y_train)

y_pred = regr.predict(X_test)

# The coefficients
print('Coefficients: \n', regr.coef_)
# The mean squared error
print("Mean squared error: %.2f"
      % mean_squared_error(Y_test, y_pred))
# Explained variance score: 1 is perfect prediction
print('Variance score: %.2f' % r2_score(Y_test, y_pred))
print(regr.get_params())

def plot_results(ind, dep, colorplot):
    """Generate Plots
    """
    # plot
    fig = plt.figure()
    plt.scatter(ind, dep, c=colorplot, s=6,
                cmap=plt.cm.get_cmap('plasma', len(set(colorplot))))
    plt.colorbar()
    # titles and labels
    
    plt.show()

    return plt

err = np.abs(y_pred - Y_test)
plot_results(np.exp(X_test['PD']), err, np.exp(X_test['PD']))

