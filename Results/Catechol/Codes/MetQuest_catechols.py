# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 10:18:39 2016
This function identifies sub-networks from Pseudomonas putida
genome scale metabolic model for the specified set of source
and target molecules, for a given cutoff of 25.

@author: aarthi
"""

from __future__ import division
from networkx import read_gpickle, get_node_attributes
from cPickle import load
import itertools
import time
from generate_partitions import generate_partitions
from operator import mul
import os
from forward_pass import forward_pass
import math

tic = time.clock()

def calculate_biomass_components(G, src, tar1, alwaysexist1, cutoff, filenames):#, gname2):
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
        for rxns in status_dict: #Can you parallalize this?
            if set(pred(rxns)).issubset(dagsfound2):
                metsrequired = list(set(pred(rxns)) - alwaysexist1) # This variable avoids the generation of partitions pertaining to seed metabolites
                shortestpathlist = []
                if metsrequired:# and goahead == 'Yes':
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
                                    if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]): #, set(succ(rxns)).isdisjoint(set(tar))
                                        flag = 'NA'
                                        #print templen, 'Loop1'
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
                                            #print partitions,mets
                                            #num += 1
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
                                #print partitions,metsrequired
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
                                        #print templen
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
                                            reactants1 = set([]) #This is to check for cycles
                                            for rxnentry in prodvec:
                                                for individualele in rxnentry:
                                                    if individualele not in rxncomb1:
                                                        rxncomb1.add(individualele)
                                                        for metsreq in pred(individualele):
                                                            reactants1.add(metsreq)
                                            for succmets in succ(rxns):
                                                if succmets not in alwaysexist1:# and succmets not in reactants1:
                                                    if succmets in dagsfound2:
                                                        #pdb.set_trace()
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
cutoff1 = 25
gname1 = ['iJN746.gpickle']
src = ['catechol_e']
taritems = ['iJN746 fum_c']
G = read_gpickle(gname1[0])
filenames = gname1[0].split('/')[-1].split('.')[0].split('_')
print 'cutoff', cutoff1
#alwaysexist1 = set(['iMM904 ACP_c', 'iMM904 adp_c', 'iMM904 amp_c', 'iMM904 atp_c', 'iMM904 co2_c', 'iMM904 coa_c', 'iMM904 h2o_c',  'iMM904 h2o_p', 'iMM904 h_c', 'iMM904 h_p', 'iMM904 nad_c', 'iMM904 nadh_c', 'iMM904 nadp_c', 'iMM904 nadph_c', 'iMM904 pi_c', 'iMM904 pi_p', 'iMM904 ppi_c'])
alwaysexist1 = set(['catechol_e','iJN746 adp_c', 'iJN746 amp_c', 'iJN746 atp_c', 'iJN746 co2_c', 'iJN746 coa_c', 'iJN746 h2o_c', 'iJN746 h2o_p', 'iJN746 h_c', 'iJN746 h_p', 'iJN746 nad_c', 'iJN746 nadh_c', 'iJN746 nadp_c', 'iJN746 nadph_c', 'iJN746 pi_c', 'iJN746 pi_p', 'iJN746 ppi_c', 'iJN746 o2_e', 'iJN746 co2_e', 'iJN746 succoa_c' ])
for sourcemets in src:
    alwaysexist1.add(sourcemets)
G = read_gpickle(gname1[0])
succ = G.successors
pred = G.predecessors
lb, lbr, dagsfound, timetaken,me,status_dict = calculate_biomass_components(G, src, taritems, alwaysexist1, cutoff1, filenames)#, gname2)
print 'Time taken', timetaken
if taritems[0] in dagsfound:
   if cutoff1 in dagsfound[taritems[0]]:
       print 'Number of sub-networks found for target metabolites from seed set inclu. source for a cutoff', str(cutoff1), '-', len(dagsfound[taritems[0]][cutoff1])
       with open('iJN746_namemap.pickle','r') as g:
           namemap = load(g)
       if not os.path.exists('Results/'+taritems[0].split(' ')[1]):
           os.makedirs('Results/'+taritems[0].split(' ')[1])


       foldername1 = 'Results/' + taritems[0].split(' ')[1] + '/' + src[0] + '_' + 'Cutoff' + str(cutoff1) + '.txt'
       with open(foldername1, 'w') as f:
           print >> f, alwaysexist1
           for x,ele in enumerate(dagsfound[taritems[0]][cutoff1]):
               print>>f, str(x+1)
               for rxnname in list(ele):
                   print>>f, namemap[rxnname],'\t', pred(rxnname),'-->', succ(rxnname)
               print>>f, '\n'
   else:
       print taritems[0] , 'not found in ', cutoff1, 'steps'

elif taritems[0] in me:
   print taritems[0], 'found in seed, but no paths of length', str(cutoff1)
onlysrctotar = []
for idx in dagsfound[taritems[0]]:
   for items in dagsfound[taritems[0]][idx]:
       if set(succ(src[0])).intersection(items):
           onlysrctotar.append(list(items))
foldername2 = 'Results/' + taritems[0].split(' ')[1] + '/' + src[0] + '_' + 'Cutoff' +'_' +'All' + str(cutoff1) +'Onlysrctotar'+ '.txt'
with open(foldername2, 'w') as f:
   print >> f, alwaysexist1
   for x,ele in enumerate(onlysrctotar):
       print>>f, str(x+1)
       for rxnname in list(ele):
           print>>f, namemap[rxnname],'\t', pred(rxnname),'-->', succ(rxnname)
       print>>f, '\n'
