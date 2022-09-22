
import numpy as np
import pandas as pd
from copy import deepcopy
from scipy import ndimage
# from sklearn.linear_model import ElasticNet
from sklearn.linear_model import ElasticNetCV
import tellurium as te
from ansatz import ansatz
import time
import matplotlib
import importlib
matplotlib.use('TkAgg')

start = time.time()


# todo: allow for variations of the ElasticNet algorithm
def rsindy(antstring, data_frame, rounding=2, rxn_types=None):

    ansatz_rxns, rates = ansatz(len(data_frame.columns)-1, rxn_types)

    ansatzfile = antstring.split('.')[0] + '_ansatz.py'
    with open(ansatzfile, 'w') as f:
        f.write(rates)

    ansatz_r = te.loada(ansatz_rxns)
    ansatz_mat = ansatz_r.getFullStoichiometryMatrix()

    t = data_frame['time']
    cols = list(data_frame)
    x = []
    for each in cols[1:]:
        x.append(data_frame[each])

    dx = []
    dc = []
    dg = []

    for each in x:
        dt = t[1] - t[0]
        dx.append(np.diff(each) / dt)
        dc.append(np.convolve(each, [1, -1]) / dt)
        dg.append(ndimage.gaussian_filter1d(each, sigma=1, order=1, mode='wrap') / dt)

    eval_x = []
    for each in x:
        eval_x.append(each[3:-3])
    eval_dx = []
    for each in dg:
        eval_dx.append(each[3:-3])

    ansatz_module = importlib.import_module(ansatzfile.split('.')[0])

    obs = []
    for each in cols[1:]:
        obs.append(np.vstack([ansatz_module.getAnsatzConcentrations(*val)*ansatz_mat[each] for val in zip(*eval_x)]))

    full_x = np.concatenate(tuple(obs), axis=0)
    full_dx = np.concatenate(tuple(eval_dx), axis=0)

    enet = ElasticNetCV(l1_ratio=1.0, eps=1e-6, tol=1e-6, max_iter=50000000, positive=True)
    enet_fit = enet.fit(full_x, full_dx)

    params = []
    rnfit = np.round(enet_fit.coef_, rounding)
    for i, each in enumerate(rnfit):
        params.append('k' + str(i) + ' = ' + str(each))

    return [enet_fit.n_iter_, enet_fit.alpha_, ansatz_rxns, params]
