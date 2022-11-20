
# creation of ansatz reactions

from itertools import combinations, permutations, product
import random


def ansatz(n_species, rxn_types, ant_rxns, rxn_num):

    species = []
    for i in range(n_species):
        species.append('S' + str(i))

    bi_species = None
    if 'bi-uni' in rxn_types or 'uni-bi' in rxn_types or 'bi-bi' in rxn_types:
        bi_species = list(combinations(species, 2))
        bi_species_ident = []
        for each in species:
            bi_species_ident.append((each, each))
        bi_species = bi_species_ident + bi_species

    rxns = []

    if 'syn' in rxn_types:
        for each in species:
            rxns.append(' -> ' + each + ';')

    if 'deg' in rxn_types:
        for each in species:
            rxns.append(each + ' -> ;')

    if 'deg2' in rxn_types:
        for each in bi_species:
            rxns.append(each[0] + ' + ' + each[1] + ' -> ;')

    if 'uni-uni' in rxn_types:
        uni_permutations = list(permutations(species, 2))
        for each in uni_permutations:
            rxns.append(each[0] + ' -> ' + each[1] + ';')

    tri_species = None
    if 'bi-uni' in rxn_types or 'uni-bi' in rxn_types:
        tri_species = list(product(species, bi_species))

    if 'bi-uni' in rxn_types:
        for each in tri_species:
            rxns.append(each[1][0] + ' + ' + each[1][1] + ' -> ' + each[0] + ';')

    if 'uni-bi' in rxn_types:
        for each in tri_species:
            rxns.append(each[0] + ' -> ' + each[1][0] + ' + ' + each[1][1] + ';')

    if 'bi-bi' in rxn_types:
        quad_species = list(permutations(bi_species, 2))
        for each in quad_species:
            rxns.append(each[0][0] + ' + ' + each[0][1] + ' -> ' + each[1][0] + ' + ' + each[1][1] + ';')

    full_ansatz = []
    new_ansatz = []
    pot_ansatz = []

    for i, each in enumerate(rxns):
        full_ansatz.append(each.split(';')[0].strip())
        if each.split(';')[0].strip() in ant_rxns:
            new_ansatz.append(each.split(';')[0].strip())
        else:
            pot_ansatz.append(each.split(';')[0].strip())

    if rxn_num < len(full_ansatz) and rxn_num > len(ant_rxns):
        selected_ansatz = random.sample(pot_ansatz, rxn_num - len(new_ansatz))
        new_ansatz.extend(selected_ansatz)
    if rxn_num >= len(full_ansatz):
        new_ansatz = full_ansatz
    if rxn_num <= len(ant_rxns):
        new_ansatz = ant_rxns

    new_ansatz2 = []
    for each in full_ansatz:
        if each in new_ansatz:
            new_ansatz2.append(each)

    ansatz_str = ''
    rxn_index = 0

    for each in new_ansatz2:
        ansatz_str += 'r' + str(rxn_index) + ': '
        if each[0] == '-':
            ansatz_str += ' ' + each + '; \n'
        elif each[-1] == '>':
            ansatz_str += each + ' ; \n'
        else:
            ansatz_str += each + '; \n'
        rxn_index += 1

    conc_func2 = 'import numpy as np\n\ndef getAnsatzConcentrations('
    for each in species[:-1]:
        conc_func2 += each + ', '
    conc_func2 += species[-1] + '): \n    return np.array([\n'
    for each in new_ansatz2:
        each_split = each.split('->')[0].strip()
        if each_split == '':
            conc_func2 += '        ' + str(1) + ', \n'
        elif '+' in each_split:
            each_split_split = each_split.split(' + ')
            conc_func2 += '        ' + each_split_split[0] + '*' + each_split_split[1] + ', \n'
        else:
            conc_func2 += '        ' + each_split + ', \n'
    conc_func2 = conc_func2[:-3]
    conc_func2 += '])\n'

    conc_func2 = conc_func2.replace('S', 's')

    return ansatz_str, conc_func2
