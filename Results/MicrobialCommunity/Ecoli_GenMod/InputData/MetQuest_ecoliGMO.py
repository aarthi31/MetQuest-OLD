# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 10:18:39 2016
This file will generate the sub-networks between two different E. coli strains
having genetic manipulations,i.e., a single gene deletion, and will find the
metabolic exchanges happening between these two organisms.

@author: aarthi
"""


from __future__ import division
from networkx import read_gpickle, get_node_attributes
import pickle
import itertools
from generate_partitions import generate_partitions
from operator import mul
import os
from forward_pass import forward_pass
import gzip
import cPickle
import math

def calculate_biomass_components(G, src, tar1, alwaysexist1, cutoff, filenames):
    H = G.reverse()
    succ = G.successors
    pred = G.predecessors
    nattr = get_node_attributes(G, 'bipartite')
    nattr_inv = {}
    for k, v in nattr.iteritems():
        nattr_inv[v] = nattr_inv.get(v, [])
        nattr_inv[v].append(k)  #0 metabolites , 1 are reactions

    for rxns1 in nattr_inv[1]:
        if len(set(pred(rxns1)) - alwaysexist1) >= 5:
            G.remove_node(rxns1)
    lowerbound, lowerboundreaction, status_dict, me = forward_pass(G, alwaysexist1, src, 0)
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
                                    if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]):#, set(succ(rxns)).isdisjoint(set(tar))]):
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
                                                if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]):#, set(succ(rxns)).isdisjoint(set(tar))]):
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
                                    if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2))]):#, set(succ(rxns)).isdisjoint(set(tar))
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
                                                if succmets not in alwaysexist1 and succmets not in reactants1:
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
    tarname = os.path.join(filenames[0],"".join(src).replace(" ","") +'241'+'.pickle.gz')
    picfname = gzip.open(tarname, 'wb')
    cPickle.dump([dagsfound2,cutoff, alwaysexist1],  picfname)
    picfname.close()
    return lowerbound, lowerboundreaction, dagsfound2

cutoff1 = 20
gname1 = ['b2276_b3708.gpickle']
G = read_gpickle(gname1[0])
filenames = gname1[0].split('/')[-1].split('.')[0].split('_')
if not os.path.exists(filenames[0]):#+'_'+filenames[1]):
    os.mkdir(filenames[0])#+'_'+filenames[1])
alwaysexist1 = set([])
deg = G.degree()
degsort = {}
nattr = get_node_attributes(G, 'bipartite')
nattr_inv = {}
alwaysexist = set(['b2276 q8_c','b2276 glyclt_c', 'b3708 glyclt_c', 'b3708 ACP_c', 'b3708 atp_c', 'b3708 accoa_c', 'b2276 accoa_c', 'b3708 h2o_c', 'b2276 ppi_c', 'b2276 nadh_c', 'b2276 coa_c', 'b3708 coa_c', 'b2276 h2o_p', 'b2276 h_c', 'b2276 pi_c', 'b2276 adp_c', 'b2276 pi_p', 'b2276 h2o_c',
                   'b2276 ACP_c', 'b2276 atp_c', 'b3708 ppi_c', 'b3708 pi_p', 'b3708 nadp_c','b2276 ACP_c' ,'b2276 adp_c' ,'b2276 adp_c' , 'b2276 adp_c' , 'b2276 atp_c' , 'b2276 coa_c' , 'b2276 h2o_c' , 'b2276 h2o_p', 'b2276 h2o_p' , 'b2276 h_c' , 'b2276 h_p' , 'b2276 nadh_c' , 'b2276 nadp_c' , 'b2276 nadph_c' , 'b2276 pi_c' , 'b2276 pi_p' , 'b2276 ppi_c', 'b2276 nadp_c', 'b3708 nadph_c', 'b2276 h_p', 'b3708 adp_c', 'b3708 pi_c', 'b3708 nadh_c', 'b3708 h_p', 'b2276 nad_c', 'b3708 o2_c', 'b2276 nadph_c',  'b3708 h_c', 'b2276 co2_c', 'b3708 h2o_p'])
src = ['b2276 glyclt_c', 'b3708 glyclt_c']
print 'with', src, 'without co2'
taritems = []
print gname1[0]
lb, lbr, dagsfound = calculate_biomass_components(G, src, taritems, alwaysexist, cutoff1, filenames)
pred = G.predecessors
succ = G.successors
excinpaths = {}
sourcerxns = set([])
for ele in src:
    for succrxns in succ(ele):
        sourcerxns.add(succrxns)
onlysrctotar = []
for mets in dagsfound:
    excinpaths[mets] = []
    for i in dagsfound[mets]:
        for items in dagsfound[mets][i]:
            if sourcerxns.intersection(items):
                onlysrctotar.append(list(items))
                for rxns in items:
                    if 'ER' in rxns:
                        excinpaths[mets].append(rxns)
tar1 = 'b3708 ins_c'
fname_2 = 'b2276_acetate_involving_pathway' + tar1 + '.txt'
pathcount = 1

with open('b2276_b3708_namemap.pickle','r') as f:
    namemap = pickle.load(f)

with open(fname_2,'w') as f:
    for ky in dagsfound[tar1]:
        for rxns in dagsfound[tar1][ky]:
            if 'Org_b2276 ER34' in rxns:
                print>>f, 'Sub-network', pathcount
                pathcount += 1
                for item in rxns:
                    try:
                        print>>f, namemap[item],'\t', pred(item),'\t', succ(item)
                    except:
                        print>>f, item, '\t', pred(item), '\t', succ(item)
                print>>f, '\n'
metsforwhichexchappen = []
for item in excinpaths:
    if excinpaths[item]:
        metsforwhichexchappen.append(item)

uniqueexcinapths ={}
for ele in metsforwhichexchappen:
    uniqueexcinapths[ele] = set(excinpaths[ele])
fname_1 = 'exchanges_with_' + '_'.join(src)+tar1 +str(cutoff1) +'withoutco2' + '.txt'
with open(fname_1,'w') as f:
    print>>f, alwaysexist
    print>>f, 'Target Metabolite', '\t', 'Exchange reaction ID', '\t', 'Reactants', '\t', 'Products'
    for item in uniqueexcinapths:
        for rxns in uniqueexcinapths[item]:
            print>>f, item, '\t', rxns, '\t', pred(rxns), '\t', succ(rxns)
        print>>f,'\n'

for lth in dagsfound['b3708 acald_c']:
    for rxns in dagsfound['b3708 acald_c'][lth]:
        if sourcerxns.intersection(rxns):
            if 'Org_b2276 ER34' in rxns:
                print rxns
