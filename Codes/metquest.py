# -*- coding: utf-8 -*-
"""

@author: aarthi
"""
from __future__ import division
from forward_pass import forward_pass
from generate_partitions import generate_partitions
from operator import mul
import time
import math
import itertools
from networkx import get_node_attributes

tic = time.clock()
def find_pathways(G, src, target, seedmets, cutoff):
    H = G.reverse()
    succ = G.successors
    pred = G.predecessors
    nattr = get_node_attributes(G, 'bipartite')
    nattr_inv = {}
    for k, v in nattr.iteritems():
        nattr_inv[v] = nattr_inv.get(v, [])
        nattr_inv[v].append(k)  #0 metabolites , 1 are reactions
    for rxns1 in nattr_inv[1]:
        if len(set(pred(rxns1)) - seedmets) >= 5:
            G.remove_node(rxns1)
    lowerbound, lowerboundreaction, status_dictfp, scope = forward_pass(G, seedmets, src, 0)
    #if target[0] not in G:
    #    status_dict = status_dictfp.copy()
    #    print target, "Since target is not found in G, only the forward pass is used"
    #else:
    lowerboundrev, lowerboundreactionrev, status_dictrev,scope_rev = forward_pass(H,scope, target, 1, src)
    status_dict = list(set(status_dictfp).intersection(set(status_dictrev)))
    pathway_table = {}
    for item in list(seedmets):
        pathway_table[item] = {0: ''}
    for rxns in status_dict:
        if set(pred(rxns)).issubset(seedmets):
            for succmets in succ(rxns):
                if succmets not in pathway_table and succmets not in seedmets:
                    pathway_table[succmets] = {1:[]}
            for succmets in succ(rxns):
                if succmets not in seedmets:
                    pathway_table[succmets][1].append(set([rxns]))
    for jdx in range(2, cutoff+1):
        for rxns in status_dict:
            if set(pred(rxns)).issubset(pathway_table):
                metsrequired = list(set(pred(rxns)) - seedmets)
                shortestpathlist = []
                if metsrequired:
                    for predmets in metsrequired:
                        shortestpathlist.append(min(lowerbound[predmets]))
                    for val in range(jdx-1, len(metsrequired)*(jdx-2)+1):
                        temp = int(math.floor(val/((jdx-1))))
                        for idx in range(1, temp+1):
                            for combined_value in itertools.combinations(metsrequired, idx):
                                mets_participant = list(set(metsrequired) - set(combined_value))
                                temprxnlist = []
                                templen = {}

                                for metabolites in list(combined_value):
                                    if metabolites in pathway_table:
                                        if jdx-1 in pathway_table[metabolites]:
                                            temprxnlist.append([[rxns]])
                                            templen[metabolites] = len(pathway_table[metabolites][jdx-1])
                                    else:
                                        break
                                flag = ''
                                if templen:
                                    for mets in list(mets_participant):
                                        templen[mets] = len(pathway_table[mets])
                                    if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(pathway_table))]):
                                        flag = 'NA'

                                    else:
                                        for mets in combined_value:
                                            if jdx-1 in pathway_table[mets]:
                                                temprxnlist.append(pathway_table[mets][jdx-1])
                                            else:
                                                flag = 'NA'
                                                break
                                    if flag != 'NA':
                                        tempshortestpath = []
                                        for varmet in list(mets_participant):
                                            tempshortestpath.append(min(lowerbound[varmet]))

                                        all_partitions_1 = generate_partitions(jdx, tempshortestpath, val-((jdx-1)*idx))

                                        for partitions in all_partitions_1:
                                            counter = 0
                                            for varmetidx in range(len(mets_participant)):
                                                if mets_participant[varmetidx] in pathway_table:
                                                    if partitions[varmetidx] in pathway_table[mets_participant[varmetidx]]:
                                                        counter += 1
                                                        templen[mets_participant[varmetidx]] = len(pathway_table[mets_participant[varmetidx]][partitions[varmetidx]])
                                            if counter == len(mets_participant):
                                                if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(pathway_table))]): #, set(succ(rxns)).isdisjoint(set(tar))
                                                    break
                                                else:
                                                    temprxnlist1 = temprxnlist[:]
                                                    for varmetidx in range(len(mets_participant)):
                                                        temprxnlist1.append(pathway_table[mets_participant[varmetidx]][partitions[varmetidx]])
                                                    for rxnunion in itertools.product(*temprxnlist1):
                                                        reaction_combntn = set([])
                                                        reactants = set([]) #This is to check for cycles
                                                        for rxnentry in rxnunion:
                                                            for individualele in rxnentry:
                                                                if individualele not in reaction_combntn:
                                                                    reaction_combntn.add(individualele)
                                                                    for metsreq in pred(individualele):
                                                                        reactants.add(metsreq)
                                                        for succmets in succ(rxns):
                                                            if succmets not in seedmets:
                                                                if succmets in pathway_table:
                                                                    if len(reaction_combntn) in pathway_table[succmets]:
                                                                        try:
                                                                            if pathway_table[succmets][len(reaction_combntn)].index(reaction_combntn):
                                                                                pass
                                                                        except:
                                                                                pathway_table[succmets][len(reaction_combntn)].append(reaction_combntn)
                                                                    else:
                                                                        pathway_table[succmets].update({len(reaction_combntn): [reaction_combntn]})
                                                                else:
                                                                    pathway_table[succmets] = {len(reaction_combntn): [reaction_combntn]}
                    for val1 in range(len(metsrequired)*(jdx-2)+1, len(metsrequired)*(jdx-1)+1):
                        all_partitions_2 = generate_partitions(jdx, shortestpathlist, val1)
                        for partitions in all_partitions_2:
                                #print partitions,metsrequired
                                temprxnlist2 = []
                                templen = {}
                                flag1 = ''
                                for item in range(len(metsrequired)):
                                    if metsrequired[item] in pathway_table:
                                        if partitions[item] in pathway_table[metsrequired[item]]:
                                            templen[metsrequired[item]]=len(pathway_table[metsrequired[item]][partitions[item]])
                                if templen:
                                    if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(pathway_table))]): #, set(succ(rxns)).isdisjoint(set(tar))
                                        flag1 = 'NA'
                                        #print templen
                                    else:
                                        for item in range(len(metsrequired)):
                                            if partitions[item] in pathway_table[metsrequired[item]]:
                                                temprxnlist2.append(pathway_table[metsrequired[item]][partitions[item]])
                                            else:
                                                flag1 = 'NA'
                                    if flag1!='NA':
                                        temprxnlist2.append([[rxns]])
                                        for rxnunion in itertools.product(*temprxnlist2):
                                            reaction_combntn = set([])
                                            reactants1 = set([]) #This is to check for cycles
                                            for rxnentry in rxnunion:
                                                for individualele in rxnentry:
                                                    if individualele not in reaction_combntn:
                                                        reaction_combntn.add(individualele)
                                                        for metsreq in pred(individualele):
                                                            reactants1.add(metsreq)
                                            for succmets in succ(rxns):
                                                if succmets not in seedmets:# and succmets not in reactants1:
                                                    if succmets in pathway_table:
                                                        #pdb.set_trace()
                                                        if len(reaction_combntn) in pathway_table[succmets]:
                                                                try:
                                                                    if pathway_table[succmets][len(reaction_combntn)].index(reaction_combntn):
                                                                        pass
                                                                except:
                                                                    pathway_table[succmets][len(reaction_combntn)].append(reaction_combntn)
                                                        else:
                                                            pathway_table[succmets].update({len(reaction_combntn): [reaction_combntn]})
                                                    else:
                                                        pathway_table[succmets] = {len(reaction_combntn): [reaction_combntn]}

        print 'Pathway length', jdx
        print 'Number of metabolites for which pathways have been found', len(pathway_table.keys())

    toc = time.clock()
    timetaken = toc - tic
    return lowerbound, lowerboundreaction, pathway_table, timetaken, scope, status_dict
