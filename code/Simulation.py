import cobra
from cobra.io import read_sbml_model
from cobra import Reaction, Metabolite
import numpy as np
import pandas as pd
import copy


def simulateGrowth(model, alpha):
    tmpmodel = model.copy()

    # max growth
    tmpmodel.objective = 'BIOMASS_SCO_tRNA'
    sol = tmpmodel.optimize()

    # fix growth and max product
    tmpgrow = sol.objective_value * alpha
    tmpmodel.reactions.get_by_id('BIOMASS_SCO_tRNA').bounds = (tmpgrow, tmpgrow)
    tmpmodel.objective = 'EX_spiA_e'
    heme_max = tmpmodel.optimize()
    print(heme_max.objective_value)

    # fix product and perform pFBA
    tmpmodel.reactions.get_by_id('EX_spiA_e').bounds = (heme_max.objective_value * 0.999,
                                                        heme_max.objective_value * 0.999)
    sol_pfba = cobra.flux_analysis.pfba(tmpmodel)
    return sol_pfba


def div_(start_, end_, num):
    l = []
    interv = (end_ - start_) / num
    co = 0
    while co < num + 1:
        l.append(start_ + co * interv)
        co += 1
    return l


def compare_substrate(model):
    flux_WT = simulateGrowth(model, 1)
    aex = flux_WT.fluxes.loc['BIOMASS_SCO_tRNA']
    alpha = div_(0.2 * aex, aex * 0.8, 15)
    v_matrix = []
    k_matrix = []
    for a in alpha:
        tmpflux = simulateGrowth(model, a)
        try:
            v_matrix.append(tmpflux.fluxes)
            k_matrix.append(np.asarray(tmpflux.fluxes) /
                            np.asarray(flux_WT.fluxes))
        except AttributeError:
            v_matrix.append(tmpflux)
            k_matrix.append(tmpflux)
        print(a)
        print(alpha.index(a))
    return v_matrix, k_matrix, alpha, flux_WT


model = read_sbml_model('../data/sco-GEM.xml')
model.reactions.get_by_id('BIOMASS_SCO').bounds = (0, 0)

# add C4 rxns and mets
forosamine_c = Metabolite('forosamine_c', compartment='c')
model.add_metabolites(forosamine_c)

agl_c = Metabolite('agl_c', compartment='c')
model.add_metabolites(agl_c)

psa_c = Metabolite('psa_c', compartment='c')
model.add_metabolites(psa_c)

spiA_c = Metabolite('spiA_c', compartment='c')
model.add_metabolites(spiA_c)

spiA_e = Metabolite('spiA_e', compartment='e')
model.add_metabolites(spiA_e)


forosamine_syth = Reaction('forosamine_syth')
forosamine_syth.add_metabolites({
    model.metabolites.get_by_id('g1p_c'): -1,
    model.metabolites.get_by_id('dtdp_c'): -1,
    model.metabolites.get_by_id('amet_c'): -2,
    model.metabolites.get_by_id('forosamine_c'): 1,
    model.metabolites.get_by_id('ahcys_c'): 2})
forosamine_syth.gene_reaction_rule = 'forosamine_syth'
model.add_reactions([forosamine_syth])


agl_syth = Reaction('agl_syth')
agl_syth.add_metabolites({
    model.metabolites.get_by_id('malcoa_c'): -2,
    model.metabolites.get_by_id('mmcoa__R_c'): -9,
    model.metabolites.get_by_id('agl_c'): 1,
    model.metabolites.get_by_id('coa_c'): 11,
    model.metabolites.get_by_id('co2_c'): 9,
    model.metabolites.get_by_id('h2o_c'): 2})
agl_syth.gene_reaction_rule = 'agl_syth'
model.add_reactions([agl_syth])


psa_syth = Reaction('psa_syth')
psa_syth.add_metabolites({
    model.metabolites.get_by_id('agl_c'): -1,
    model.metabolites.get_by_id('dtdprmn_c'): -1,
    model.metabolites.get_by_id('amet_c'): -3,
    model.metabolites.get_by_id('ahcys_c'): 3,
    model.metabolites.get_by_id('psa_c'): 1,
    model.metabolites.get_by_id('h2o_c'): 1,
    model.metabolites.get_by_id('dtdp_c'): 2})
psa_syth.gene_reaction_rule = 'psa_syth'
model.add_reactions([psa_syth])


spiA_syth = Reaction('spiA_syth')
spiA_syth.add_metabolites({
    model.metabolites.get_by_id('psa_c'): -1,
    model.metabolites.get_by_id('forosamine_c'): -1,
    model.metabolites.get_by_id('spiA_c'): 1,
    model.metabolites.get_by_id('h2o_c'): 1,
    model.metabolites.get_by_id('dtdp_c'): 1})
spiA_syth.gene_reaction_rule = 'spiA_syth'
model.add_reactions([spiA_syth])


spiA_exchange = Reaction('spiA_exchange')
spiA_exchange.bounds = (-1000, 1000)
spiA_exchange.add_metabolites({
    model.metabolites.get_by_id('spiA_c'): -1,
    model.metabolites.get_by_id('spiA_e'): 1})
spiA_exchange.gene_reaction_rule = 'spiA_exchange'
model.add_reactions([spiA_exchange])

spiA_transport = Reaction('EX_spiA_e')
spiA_transport.bounds = (0, 1000)
spiA_transport.add_metabolites({
    model.metabolites.get_by_id('spiA_e'): -1})
#spiA_transport.gene_reaction_rule = 'spiA_transport'

model.add_reactions([spiA_transport])

model.objective = model.reactions.get_by_id(spiA_transport.id)
solu = model.optimize()

model.objective = model.reactions.get_by_id('BIOMASS_SCO_tRNA')
solu2 = model.optimize()

# Follow Olena P. Ishchuk's paper (https://doi.org/10.1073/pnas.2108245119)
v_matrix, k_matrix, alpha, flux_WT = compare_substrate(model)
all_flux = pd.DataFrame(v_matrix)
all_flux = all_flux.T
all_k = pd.DataFrame(k_matrix)
all_k = all_k.T
all_k.index = all_flux.index

f1 = list()
for h in all_flux.index:
    if bool(model.reactions.get_by_id(h).gene_reaction_rule == ''):
        f1.append(h)

all_flux.drop(index=f1, inplace=True)
all_k.drop(index=f1, inplace=True)
print()
print('delete all nan')
all_flux.dropna(axis=0, how='all', inplace=True)
all_k.dropna(axis=0, how='all', inplace=True)
print(len(all_k.index))

print('nan --> 1')
all_flux.fillna(1, inplace=True)
all_k.fillna(1, inplace=True)
print(len(all_k.index))
# inf --> 1000
all_k.replace([np.inf, -np.inf], 1000, inplace=True)

print('delete inconsistent rxn')
f5 = []
for i in all_k.index:
    big = any(all_k.loc[i, :] > 1)
    small = any(all_k.loc[i, :] < 1)
    if big == small == True:
        f5.append(i)
all_k.drop(index=f5, inplace=True)
print(len(all_k.index))


# 删除最大化生产小于0.1的k
# f7 = []
# for i in all_k.index:
#     if 0.02 < np.abs(all_k.loc[i, 15]) < 500:
#         f7.append(i)
# all_k.drop(index=f7, inplace=True)
# print(len(all_k.index))

# order k(descend)
orderd_k = copy.deepcopy(all_k)
orderd_k['meank'] = orderd_k.mean(axis=1)
orderd_k = orderd_k.sort_values(by=["meank"], ascending=False)
print('get gene')
cons_g = {}
single_g = {}
for ge in model.genes:
    # 取基因的k值
    ge_rxn = [r.id for r in ge.reactions if r.id in orderd_k.index]
    ge_k = [orderd_k.loc[r.id, 'meank'] for r in ge.reactions if r.id in orderd_k.index]
    # delete inconsistent gene
    if ge_k != []:
        gk = np.asarray(ge_k)
        if all(gk > 1) == True or all(gk < 1) == True:
            cons_g[ge.id] = np.mean(ge_k)
print(len(cons_g))
print(len(single_g))

cons_g_f = {}
for k, v in cons_g.items():
    if (v < 0.002) or (v > 500):
        cons_g_f[k] = v


cons_g_f = pd.Series(cons_g_f)
cons_g_f.sort_values(ascending=False, inplace=True)
cons_g_f.to_excel('../output/targets.xlsx')
