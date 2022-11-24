
import numpy as np
import pandas as pd
from scipy import ndimage
# from sklearn.linear_model import ElasticNet
from sklearn.linear_model import ElasticNetCV
import tellurium as te
import matplotlib
import importlib
import sys
matplotlib.use('TkAgg')


# todo: allow for variations of the ElasticNet algorithm
def rsindy(antstring, t_span, steps, rounding=2):

    ant_rxns = []
    with open(antstring, 'r') as antfile:
        lines = antfile.readlines()
        for line in lines:
            if '->' in line:
                line_split = line.split(':')[1].split(';')[0].strip()
                line_split = line_split.split('->')
                line_split[0] = line_split[0].strip()
                line_split[1] = line_split[1].strip()
                if '+' in line_split[0]:
                    line_split[0] = line_split[0].split(' + ')
                    line_split[0][0] = int(line_split[0][0][1:])
                    line_split[0][1] = int(line_split[0][1][1:])
                    line_split[0].sort()
                    line_split[0][0] = 'S' + str(line_split[0][0])
                    line_split[0][1] = 'S' + str(line_split[0][1])
                    line_split[0] = ' + '.join(line_split[0])

                if '+' in line_split[1]:
                    line_split[1] = line_split[1].split(' + ')
                    line_split[1][0] = int(line_split[1][0][1:])
                    line_split[1][1] = int(line_split[1][1][1:])
                    line_split[1].sort()
                    line_split[1][0] = 'S' + str(line_split[1][0])
                    line_split[1][1] = 'S' + str(line_split[1][1])
                    line_split[1] = ' + '.join(line_split[1])

                line_split = ' -> '.join(line_split)

                ant_rxns.append(line_split.strip())

    r = te.loada(antstring)
    cols = r.getFloatingSpeciesIds() + r.getBoundarySpeciesIds()

    # print()
    # print('ANSATZ_RXNS---------------------')
    # print(ansatz_rxns)
    # print()
    # print('RATES---------------------------')
    # print(rates)
    # print()
    # print('ANSATZ_MAT----------------------')
    # print(ansatz_mat)

    cols.insert(0, 'time')

    sim = r.simulate(0, int(t_span), int(steps), selections=cols)
    sim_df = pd.DataFrame(sim, columns=cols)

    # print()
    # print('SIM_DF-DATA---------------------')
    # print(sim_df)

    t = sim_df['time']
    x = []
    for each in cols[1:]:
        x.append(sim_df[each])

    # print()
    # print(x[0])
    # print()
    # print(x[1])

    dx = []
    dc = []
    dg = []

    for each in x:
        dt = t[1] - t[0]
        dx.append(np.diff(each) / dt)
        dc.append(np.convolve(each, [1, -1]) / dt)
        dg.append(ndimage.gaussian_filter1d(each, sigma=1, order=1, mode='wrap') / dt)

    # print()
    # print('DG----------------------------')
    # print(dg)

    ansatz_function = antstring.split('.')[0] + '_ansatz.py'
    ansatz_rxns = antstring.split('.')[0] + '_ansatz.ant'

    ansatz_r = te.loada(ansatz_rxns)
    ansatz_mat = ansatz_r.getFullStoichiometryMatrix()

    eval_x = []
    for each in x:
        eval_x.append(each[3:-3])
    eval_dx = []
    for each in dg:
        eval_dx.append(each[3:-3])
    # print(ansatzfile)
    ansatz_module = importlib.import_module(ansatz_function.split('.')[0])

    # print()
    # print('-------------------------------')
    # print(eval_x[0])
    # print(eval_x[1])
    # print(list(zip(*eval_x)))

    obs = []
    for each in cols[1:]:
        obs.append(np.vstack([ansatz_module.getAnsatzConcentrations(*val)*ansatz_mat[each] for val in zip(*eval_x)]))

    # print()
    # print('OBS----------------------------')

    np.set_printoptions(threshold=sys.maxsize)
    # for each in obs:
    #     print(each)

    full_x = np.concatenate(tuple(obs), axis=0)
    full_dx = np.concatenate(tuple(eval_dx), axis=0)

    enet = ElasticNetCV(l1_ratio=1.0, eps=1e-6, tol=1e-4, max_iter=50000000, positive=True)
    enet_fit = enet.fit(full_x, full_dx)

    params = []
    rnfit = np.round(enet_fit.coef_, rounding)
    for i, each in enumerate(rnfit):
        params.append('k' + str(i) + ' = ' + str(each))

    return [enet_fit.n_iter_, enet_fit.alpha_, ansatz_rxns, params]
