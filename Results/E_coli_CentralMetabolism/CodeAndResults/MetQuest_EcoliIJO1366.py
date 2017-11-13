# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 10:18:39 2016
This function identifies sub-networks from E. coli genome-scale metabolic
model for the specified set of source and target molecules, for a given cutoff of 15.

@author: aarthi
"""


from __future__ import division
from networkx import read_gpickle, get_node_attributes
from pickle import dump
from cPickle import load
import itertools
import time
from generate_partitions import generate_partitions
from operator import mul
import os
from forward_pass import forward_pass
import math

tic = time.clock()

def calculate_biomass_components(G, src, tar1, alwaysexist1, cutoff, filenames):
    H = G.reverse()
    nattr = get_node_attributes(G, 'bipartite')
    nattr_inv = {}
    for k, v in nattr.iteritems():
        nattr_inv[v] = nattr_inv.get(v, [])
        nattr_inv[v].append(k)  #0 metabolites , 1 are reactions

    for rxns1 in nattr_inv[1]:
        if len(set(pred(rxns1)) - alwaysexist1) >= 5:
            G.remove_node(rxns1)
    tar= tar1
    lowerbound, lowerboundreaction, status_dictfp, me, alwaysexist_first = forward_pass(G, alwaysexist1, src, 0)
    lowerboundrev, lowerboundreactionrev, status_dictrev, me1, ae2 = forward_pass(H, me, tar, 1, src)
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
                    dagsfound2[succmets][1].append(set([rxns]))# = {1: [set([rxns])]}
    for i in range(2, cutoff+1):
        for rxns in status_dict: 
            if set(pred(rxns)).issubset(dagsfound2):
                metsrequired = list(set(pred(rxns)) - alwaysexist1) # This variable avoids the generation of partitions pertaining to seed metabolites
                shortestpathlist = []
                if metsrequired:
                    for predmets in metsrequired:
                        shortestpathlist.append(min(lowerbound[predmets]))
    ##==============================================================================
    ## First loop
    ##==============================================================================
                    for val in range(i-1, len(metsrequired)*(i-2)+1): #range does not include the end value
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
                                    if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]): 
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
                                            if counter == len(variablemets): 
                                                if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]): 
                                                    break 
                                                else:
                                                    temprxnlist1 = temprxnlist[:]
                                                    for varmetidx in range(len(variablemets)):
                                                        temprxnlist1.append(dagsfound2[variablemets[varmetidx]][partitions[varmetidx]])
                                                    for prodvec in itertools.product(*temprxnlist1):
                                                        rxncomb = set([])
                                                        reactants = set([])
                                                        for rxnentry in prodvec:
                                                            for individualele in rxnentry:
                                                                if individualele not in rxncomb:
                                                                    rxncomb.add(individualele)
                                                                    for metsreq in pred(individualele):
                                                                        reactants.add(metsreq)
                                                        for succmets in succ(rxns):
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
                                    if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]): 
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
                                            reactants1 = set([]) 
                                            for rxnentry in prodvec:
                                                for individualele in rxnentry:
                                                    if individualele not in rxncomb1:
                                                        rxncomb1.add(individualele)
                                                        for metsreq in pred(individualele):
                                                            reactants1.add(metsreq)
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
    return lowerbound, lowerboundreaction, dagsfound2, toc-tic,me, status_dict
if not os.path.exists('Results'):
    os.makedirs('Results')
cutoff1 = 15
gname1 = ['iJO1366.gpickle']
src = ['glc_DASH_D_e']
taritems = ['iJO1366 pyr_c']
G = read_gpickle(gname1[0])
filenames = gname1[0].split('/')[-1].split('.')[0].split('_')
print 'cutoff', cutoff1
alwaysexist1 = set(['glc_DASH_D_e','iJO1366 ACP_c', 'iJO1366 adp_c', 'iJO1366 amp_c', 'iJO1366 atp_c', 'iJO1366 co2_c', 'iJO1366 coa_c', 'iJO1366 h2o_c', 'iJO1366 h2o_p', 'iJO1366 h_c', 'iJO1366 h_p', 'iJO1366 nad_c', 'iJO1366 nadh_c', 'iJO1366 nadp_c', 'iJO1366 nadph_c', 'iJO1366 pi_c', 'iJO1366 pi_p', 'iJO1366 ppi_c'])
succ = G.successors
pred = G.predecessors
lb, lbr, dagsfound, timetaken,me,status_dict = calculate_biomass_components(G, src, taritems, alwaysexist1, cutoff1, filenames)
print 'Time taken', timetaken
with open('iJO1366_namemap.pickle', 'r') as f:
    namemap = load(f)
fname1 = 'Results/' + src[0] + '_' + taritems[0] + '_' + str(cutoff1) + '.txt'
pathnumcount = 0
with open(fname1,'w') as f:
   for idx in dagsfound[taritems[0]]:
       print>>f, idx, '\n'
       for items in dagsfound[taritems[0]][idx]:
           pathnumcount += 1
           print>>f, pathnumcount
           for elements in list(items):
               print>>f, namemap[elements], '\t', pred(elements),'-->', succ(elements)
           print>>f, '\n'

onlysrctotar = []
for idx in dagsfound[taritems[0]]:
    for items in dagsfound[taritems[0]][idx]:
        if set(succ(src[0])).intersection(items):
            onlysrctotar.append(list(items))
kegg_ids_paths = {}
for x,elements1 in enumerate(onlysrctotar):
    currentpath = set([])
    for elements in elements1:
        currentpath.add(namemap[elements])
    if len(currentpath) in kegg_ids_paths:
        kegg_ids_paths[len(currentpath)].append(currentpath)
    else:
        kegg_ids_paths[len(currentpath)] = [currentpath]

fname2 = 'Results/' + 'only' + src[0] + '_' + taritems[0] + '_' + str(cutoff1) + '.txt'
with open(fname2,'w') as f:
   for x,elements1 in enumerate(onlysrctotar):
       print>>f, x+1
       for elements in elements1:
           print>>f, namemap[elements], '\t', pred(elements),'-->',succ(elements)
       print>>f,'\n'
