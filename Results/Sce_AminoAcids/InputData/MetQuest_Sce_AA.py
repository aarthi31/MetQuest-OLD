# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 10:18:39 2016
This function identifies sub-networks from E. coli genome-scale metabolic
model for the specified set of source and target molecules, for a given cutoff of 15.

@author: aarthi
"""


from __future__ import division
from networkx import read_gpickle, get_node_attributes
import os
from pickle import dump
import itertools
import time
from generate_partitions import generate_partitions
import cPickle
import gzip
from operator import mul
from forward_pass import forward_pass
import math


tic = time.clock()
def calculate_biomass_components(G, src, tar1, alwaysexist1, cutoff, filenames):
    H = G.reverse()
    nattr = get_node_attributes(G, 'bipartite')
    nattr_inv = {}
    for k, v in nattr.iteritems():
        nattr_inv[v] = nattr_inv.get(v, [])
        nattr_inv[v].append(k)  # 0 metabolites , 1 are reactions
    for rxns1 in nattr_inv[1]:
        if len(set(pred(rxns1)) - alwaysexist1) >= 5:
            G.remove_node(rxns1)
    tar= tar1
    lowerbound, lowerboundreaction, status_dictfp, me, alwaysexist_first = \
        forward_pass(G, alwaysexist1, src, 0)
    lowerboundrev, lowerboundreactionrev, status_dictrev, me1, ae1 = \
        forward_pass(H, me, tar, 1, src)
    status_dict = list(set(status_dictfp).intersection(set(status_dictrev)))
    dagsfound2 = {}
    for item in list(alwaysexist1):
        dagsfound2[item] = {0: ''}
    for rxns in status_dict:
        if set(pred(rxns)).issubset(alwaysexist1):
            for succmets in succ(rxns):
                if succmets not in dagsfound2 and succmets not in alwaysexist1:
                    dagsfound2[succmets] = {1:[]}
            for succmets in succ(rxns):
                if succmets not in alwaysexist1:
                    dagsfound2[succmets][1].append(set([rxns]))
    for rxns in status_dict:
        if set(pred(rxns)).issubset(alwaysexist1):
            for succmets in succ(rxns):
                if succmets not in alwaysexist1:
                    dagsfound2[succmets] = {1: [set([rxns])]}
    for i in range(2, cutoff+1):
        for rxns in status_dict:
            if set(pred(rxns)).issubset(dagsfound2):
                metsrequired = list(set(pred(rxns)) - alwaysexist1)
                shortestpathlist = []
                if metsrequired:# and goahead == 'Yes':
                    for predmets in metsrequired:
                        shortestpathlist.append(min(lowerbound[predmets]))
    # ========================================================================
    # First loop
    # ========================================================================
                    for val in range(i-1, len(metsrequired)*(i-2)+1):
                        temp = int(math.floor(val/((i-1))))
                        for idx in range(1, temp+1):
                            for combivalue in itertools.combinations(metsrequired, idx):
                                variablemets = list(set(metsrequired) - set(combivalue))
                                temprxnlist = []
                                templen = {}

                                for metabolites in list(combivalue):
                                    if metabolites in dagsfound2:
                                        if i-1 in dagsfound2[metabolites]:
                                            temprxnlist.append([[rxns]])
                                            templen[metabolites]=len(dagsfound2[metabolites][i-1])
                                    else:
                                        break
                                flag = ''
                                if templen:
                                    for mets in list(variablemets):
                                        templen[mets] = len(dagsfound2[mets])
                                    if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]): #, set(succ(rxns)).isdisjoint(set(tar))
                                        flag = 'NA'
                                    else:
                                        for mets in combivalue:
                                            if i-1 in dagsfound2[mets]:
                                                temprxnlist.append(dagsfound2[mets][i-1])
                                            else:
                                                flag = 'NA'
                                                break
                                    if flag != 'NA':
                                        tempshortestpath = []
                                        for varmet in list(variablemets):
                                            tempshortestpath.append(min(lowerbound[varmet]))
                                        par = generate_partitions(i, tempshortestpath, val-((i-1)*idx))
                                        for partitions in par:
                                            counter = 0
                                            for varmetidx in range(len(variablemets)):
                                                if variablemets[varmetidx] in dagsfound2:
                                                    if partitions[varmetidx] in dagsfound2[variablemets[varmetidx]]:
                                                        counter += 1
                                                        templen[variablemets[varmetidx]] = len(dagsfound2[variablemets[varmetidx]][partitions[varmetidx]])
                                            if counter == len(variablemets): #because we ae appending the value to templen, we need another flag to check if the other metabolite is present or not
                                                if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]): #, set(succ(rxns)).isdisjoint(set(tar))
                                                    break #or pass
                                                else:
                                                    temprxnlist1 = temprxnlist[:]
                                                    for varmetidx in range(len(variablemets)):
                                                        temprxnlist1.append(dagsfound2[variablemets[varmetidx]][partitions[varmetidx]])
                                                    for prodvec in itertools.product(*temprxnlist1):
                                                        rxncomb = set([])
                                                        reactants = set([]) #This is to check for cycles
                                                        for rxnentry in prodvec:
                                                            for individualele in rxnentry:
                                                                if individualele not in rxncomb:
                                                                    rxncomb.add(individualele)
                                                                    for metsreq in pred(individualele):
                                                                        reactants.add(metsreq)
                                                        for succmets in succ(rxns):
#                                                            if succmets in reactants:
#                                                                pdb.set_trace()
                                                            if succmets not in alwaysexist1 and succmets not in reactants:
                                                                if succmets in dagsfound2:
                                                                    if len(rxncomb) in dagsfound2[succmets]:
                                                                        try:
                                                                            if dagsfound2[succmets][len(rxncomb)].index(rxncomb):
                                                                                pass
                                                                        except:
                                                                                dagsfound2[succmets][len(rxncomb)].append(rxncomb)
                                                                    else:
                                                                        dagsfound2[succmets].update({len(rxncomb): [rxncomb]})
                                                                else:
                                                                    dagsfound2[succmets] = {len(rxncomb): [rxncomb]}
    #==============================================================================
    # Second loop
    #==============================================================================
                    for val1 in range(len(metsrequired)*(i-2)+1, len(metsrequired)*(i-1)+1):
                        par1 = generate_partitions(i, shortestpathlist,val1)
                        for partitions in par1:
                            temprxnlist2 = []
                            templen = {}
                            flag1 = ''
                            for item in range(len(metsrequired)):
                                if metsrequired[item] in dagsfound2:
                                    if partitions[item] in dagsfound2[metsrequired[item]]:
                                        templen[metsrequired[item]]=len(dagsfound2[metsrequired[item]][partitions[item]])
                            if templen:
                                if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]): #, set(succ(rxns)).isdisjoint(set(tar))
                                    flag1 = 'NA'
                                else:
                                    for item in range(len(metsrequired)):
                                        if partitions[item] in dagsfound2[metsrequired[item]]:
                                            temprxnlist2.append(dagsfound2[metsrequired[item]][partitions[item]])
                                        else:
                                            flag1 = 'NA'
                                if flag1!='NA':
                                    temprxnlist2.append([[rxns]])
                                    for prodvec in itertools.product(*temprxnlist2):
                                        rxncomb1 = set([])
                                        for rxnentry in prodvec:
                                            for individualele in rxnentry:
                                                if individualele not in rxncomb1:
                                                    rxncomb1.add(individualele)
                                        for succmets in succ(rxns):
                                            if succmets not in alwaysexist1:
                                                if succmets in dagsfound2:
                                                    if len(rxncomb1) in dagsfound2[succmets]:
                                                            try:
                                                                if dagsfound2[succmets][len(rxncomb1)].index(rxncomb1):
                                                                    pass
                                                            except:
                                                                dagsfound2[succmets][len(rxncomb1)].append(rxncomb1)
                                                    else:
                                                        dagsfound2[succmets].update({len(rxncomb1): [rxncomb1]})
                                                else:
                                                    dagsfound2[succmets] = {len(rxncomb1): [rxncomb1]}

        print 'Column value', i, 'Cumulative number of metabolites', len(dagsfound2.keys())
    toc = time.clock()
    return lowerbound, lowerboundreaction, dagsfound2, toc-tic,me
cutoff1 = 30
gname1 = ['iMM904.gpickle']
G = read_gpickle(gname1[0])
filenames = gname1[0].split('/')[-1].split('.')[0].split('_')
src = ['iMM904 glc_DASH_D_c']
alwaysexist1 = set(['iMM904 glc_DASH_D_c', 'iMM904 accoa_c','iMM904 accoa_m','iMM904 ACP_c', 'iMM904 ACP_m', 'iMM904 ACP_p', 'iMM904 adp_c', 'iMM904 adp_m', 'iMM904 adp_p', 'iMM904 nh4_c','iMM904 atp_c', 'iMM904 atp_m', 'iMM904 atp_p', 'iMM904 co2_c', 'iMM904 co2_m', 'iMM904 co2_p', 'iMM904 coa_c', 'iMM904 coa_m', 'iMM904 coa_p', 'iMM904 h2o_m', 'iMM904 h2o_p', 'iMM904 h_c', 'iMM904 h_m', 'iMM904 h_p', 'iMM904 nad_c', 'iMM904 nad_m', 'iMM904 nad_p', 'iMM904 nadh_c', 'iMM904 nadh_m', 'iMM904 nadh_p', 'iMM904 nadp_c', 'iMM904 nadp_m', 'iMM904 nadp_p', 'iMM904 nadph_c', 'iMM904 nadph_m', 'iMM904 nadph_p', 'iMM904 pi_c', 'iMM904 pi_m', 'iMM904 pi_p', 'iMM904 h2o_c', 'iMM904 ppi_c', 'iMM904 ppi_m', 'iMM904 ppi_p'])
G = read_gpickle(gname1[0])
print 'Cutoff', cutoff1
taritems = ['iMM904 phe_DASH_L_c']
print taritems[0]
succ = G.successors
pred = G.predecessors
lb, lbr, dagsfound, timetaken,me = calculate_biomass_components(G, src, taritems, alwaysexist1, cutoff1, filenames)#, gname2)
print 'Time taken', timetaken
onlysrctotar = []
for idx in dagsfound[taritems[0]]:
    for items in dagsfound[taritems[0]][idx]:
        if set(succ(src[0])).intersection(items):
            onlysrctotar.append(list(items))
if not os.path.exists('Results'):
    os.makedirs('Results')
datafilename = 'Results/'+filenames[0] + '_' + src[0].split(' ')[1] + '_' + taritems[0].split(' ')[1] + '_' + str(cutoff1) +'.gz'
fp=gzip.open(datafilename,'wb')
cPickle.dump([dagsfound, alwaysexist1,onlysrctotar],fp)
fp.close()

with open('iMM904_namemap.pickle','r') as f:
    names = cPickle.load(f)

for k in dagsfound[taritems[0]].keys():
    fnames = 'Results/' + filenames[0] + '_' + src[0].split(' ')[1] + '_' + taritems[0].split(' ')[1] + '_' + str(k) +'.txt'
    count = 0
    with open(fnames,'w') as f:
        for item in onlysrctotar:
            if len(item) == k:
                count += 1
                print>>f, "Path" + " " + str(count)
                for i in item:
                    print>>f, names[i],' ', pred(i), ' -> ', succ(i)
                print>>f, '\n'