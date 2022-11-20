
import numpy as np
import pandas as pd
from copy import deepcopy
from scipy import ndimage
# from sklearn.linear_model import ElasticNet
from sklearn.linear_model import ElasticNetCV
import tellurium as te
from random_ansatz import ansatz
import matplotlib
import importlib
import sys
matplotlib.use('TkAgg')


# todo: allow for variations of the ElasticNet algorithm
def generate_ansatz(antstring, anz_rxns_num, anz_num, rxn_types=None):

    if not rxn_types:
        rxn_types = ['syn', 'deg', 'uni-uni', 'bi-uni', 'uni-bi', 'bi-bi']

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

    species_names = deepcopy(cols)

    ansatz_sets = []
    rates_sets = []
    while len(ansatz_sets) < anz_num:
        ansatz_rxns, rates = ansatz(len(species_names), rxn_types, ant_rxns, anz_rxns_num)
        if ansatz_rxns not in ansatz_sets:
            ansatz_sets.append(ansatz_rxns)
            rates_sets.append(rates)

    # print()
    # print('ANSATZ_RXNS---------------------')
    # print(ansatz_rxns)
    # print()
    # print('RATES---------------------------')
    # print(rates)
    # quit()

    for i, each in enumerate(rates_sets):
        ansatzfile = antstring.split('.')[0] + '_' + str(i) + '_ansatz.py'
        with open(ansatzfile, 'w') as f:
            f.write(each)

    for i, each in enumerate(ansatz_sets):
        ansatz_ant = antstring.split('.')[0] + '_' + str(i) + '_ansatz.ant'
        with open(ansatz_ant, 'w') as f:
            f.write(each)

