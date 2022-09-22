
# creation of ansatz reactions

from itertools import combinations, permutations, product


def ansatz(n_species, rxn_types=None):

    species = []
    for i in range(n_species):
        species.append('S' + str(i))

    if not rxn_types:
        rxn_types = ['syn', 'deg', 'uni-uni', 'bi-uni', 'uni-bi', 'bi-bi']

    bi_species = None
    if 'bi-uni' in rxn_types or 'uni-bi' in rxn_types or 'bi-bi' in rxn_types:
        bi_species = list(combinations(species, 2))
        bi_species_ident = []
        for each in species:
            bi_species_ident.append((each, each))
        bi_species = bi_species_ident + bi_species

    syn_r = []
    deg_r = []
    deg2_r = []
    uni_uni_r = []
    bi_uni_r = []
    uni_bi_r = []
    bi_bi_r = []

    conc_func = '\nimport numpy as np\n\n\ndef getAnsatzConcentrations('
    for each in species[:-1]:
        conc_func += each + ', '
    conc_func += species[-1] + '): \n    return np.array([\n'

    if 'syn' in rxn_types:
        for each in species:
            syn_r.append('-> ' + each + ';')
            conc_func += '        ' + str(1) + ', \n'

    if 'deg' in rxn_types:
        for each in species:
            deg_r.append(each + ' -> ;')
            conc_func += '        ' + each + ', \n'

    if 'deg2' in rxn_types:
        for each in bi_species:
            deg2_r.append(each[0] + ' + ' + each[1] + ' -> ;')
            conc_func += '        ' + each[0] + '*' + each[1] + ', \n'

    if 'uni-uni' in rxn_types:
        uni_permutations = list(permutations(species, 2))
        for each in uni_permutations:
            uni_uni_r.append(each[0] + ' -> ' + each[1] + ';')
            conc_func += '        ' + each[0] + ', \n'

    tri_species = None
    if 'bi-uni' in rxn_types or 'uni-bi' in rxn_types:
        tri_species = list(product(species, bi_species))

    if 'bi-uni' in rxn_types:
        for each in tri_species:
            bi_uni_r.append(each[1][0] + ' + ' + each[1][1] + ' -> ' + each[0] + ';')
            conc_func += '        ' + each[1][0] + '*' + each[1][1] + ', \n'

    if 'uni-bi' in rxn_types:
        for each in tri_species:
            uni_bi_r.append(each[0] + ' -> ' + each[1][0] + ' + ' + each[1][1] + ';')
            conc_func += '        ' + each[0] + ', \n'

    if 'bi-bi' in rxn_types:
        quad_species = list(permutations(bi_species, 2))
        for each in quad_species:
            bi_bi_r.append(each[0][0] + ' + ' + each[0][1] + ' -> ' + each[1][0] + ' + ' + each[1][1] + ';')
            conc_func += '        ' + each[0][0] + '*' + each[0][1] + ', \n'
    conc_func = conc_func[:-3]
    conc_func += '])\n'

    ansatz_str = ''

    rxn_index = 0
    for each in syn_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    for each in deg_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    for each in deg2_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    for each in uni_uni_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    for each in bi_uni_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    for each in uni_bi_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    for each in bi_bi_r:
        ansatz_str += 'r' + str(rxn_index) + ': ' + each + ' \n'
        rxn_index += 1

    conc_func = conc_func.replace('S', 's')

    return ansatz_str, conc_func
