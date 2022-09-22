
import rsindy
import numpy as np
import pandas as pd
from copy import deepcopy
import tellurium as te
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

if __name__ == "__main__":

    # todo: algorithm to set inference conditions and generate data
    modelname = 'lv.txt'
    model = te.loada(modelname)
    species = model.getFloatingSpeciesIds() + model.getBoundarySpeciesIds()
    cols = deepcopy(species)
    cols.insert(0, 'time')
    sim = model.simulate(0, 15, 200, selections=cols)
    sim_df = pd.DataFrame(sim, columns=cols)

    # It may be necessary to pass the maximum data set and instructions to subset it.
    # output = rsindy.rsindy('lv.txt')
    output = rsindy.rsindy(modelname, sim_df)

    print()
    print('number of iterations:', output[0])
    print()
    print('alpha parameter:', output[1])
    print()
    print('ansatz reactions')
    print(output[2])
    # print()
    print('inferred parameter values')
    for each in output[3]:
        print(each)
    print()
