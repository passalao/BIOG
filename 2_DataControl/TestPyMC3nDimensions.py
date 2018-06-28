from pandas_datareader import data
import pandas as pd
from io import StringIO
import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt


data = np.genfromtxt("../../SourceData/Temperatures/LR04-stack.csv", delimiter=",")
print(data[:,1])
returns=pd.Series(data[:,1], index=data[:,0])#columns=['Age', 'Delta18O', 'Sigma'])
print(returns)

with pm.Model() as sp500_model:
    nu = pm.Exponential('nu', 1/10., testval=5.)
    sigma = pm.Exponential('sigma', 1/0.02, testval=.1)
    s = pm.GaussianRandomWalk('s', sd=sigma, shape=len(returns))
    volatility_process = pm.Deterministic('volatility_process', pm.math.exp(-2*s)**0.5)
    r = pm.StudentT('r', nu=nu, sd=volatility_process, observed=returns)

with sp500_model:
    trace = pm.sample(200)

pm.traceplot(trace, varnames=['nu', 'sigma']);

fig, ax = plt.subplots(figsize=(15, 8))
returns.plot(ax=ax)
ax.plot(returns.index, 1/np.exp(trace['s',::5].T), 'C3', alpha=.03);
ax.set(title='volatility_process', xlabel='time', ylabel='volatility');
ax.legend(['S&P500', 'stochastic volatility process']);

##########################
