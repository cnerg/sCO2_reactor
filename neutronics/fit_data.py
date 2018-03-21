from sklearn.neural_network import MLPRegressor
import numpy as np
import neutronic_sweeps as ns

def load_from_csv(datafile="depl_results.csv"):
    """load the results data from a csv.
    """
    data = np.genfromtxt(datafile, delimiter=',',
            names=list(ns.parameters.keys()) + ['keff'])
    
    return data

dims = ['power', 'core_r', 'enrich', 'cool_r', 'AR', 'PD']

data = load_from_csv()
keff = data['keff']

xdata = []

for sample in data:
    sample_params = []
    for item in dims:
        sample_params.append(sample[item])
    xdata.append(sample_params)


nn = MLPRegressor()
nn.fit(xdata, keff)


xtest = [[142.417, 45.6, 0.822, 0.82378, 0.999, 1.3]]

y = nn.predict(xtest)

print(y)
