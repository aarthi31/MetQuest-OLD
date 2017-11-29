# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:46:52 2017

This code constructs the bipartite graph from the input SBML models.
@author: Aarthi Ravikrishnan
"""
import networkx as nx
import fetchreactions
import itertools
import os
import pdb
import sys
from cPickle import dump

def create_graph_with_internal_reaction(organismsdata):
    orgkeys1 = organismsdata.keys()
    G = nx.DiGraph()
    for j in range(len(orgkeys1)):
        G.add_nodes_from(organismsdata[orgkeys1[j]]['exchno'],bipartite=1)
        G.add_nodes_from(organismsdata[orgkeys1[j]]['irrrxnno'],bipartite=1)
        G.add_nodes_from(organismsdata[orgkeys1[j]]['revrxnno'],bipartite = 1)
        G.add_nodes_from(organismsdata[orgkeys1[j]]['revbacrxno'],bipartite = 1)
        irrevlhs = list(set([item for sublist in organismsdata[orgkeys1[j]]['irrevlhsnodes'] for item in sublist]))
        irrevrhs = list(set([item for sublist in organismsdata[orgkeys1[j]]['irrevrhsnodes'] for item in sublist]))
        revlhs = list(set([item for sublist in organismsdata[orgkeys1[j]]['revlhsnodes'] for item in sublist]))
        revrhs = list(set([item for sublist in organismsdata[orgkeys1[j]]['revrhsnodes'] for item in sublist]))
        G.add_nodes_from(irrevlhs,bipartite = 0)
        G.add_nodes_from(irrevrhs,bipartite = 0)
        G.add_nodes_from(revlhs,bipartite = 0)
        G.add_nodes_from(revrhs,bipartite = 0)
        for idx in range(len(organismsdata[orgkeys1[j]]['irrrxnno'])):
            for idx2 in range(len(organismsdata[orgkeys1[j]]['irrevlhsnodes'][idx])):
                G.add_edges_from([(organismsdata[orgkeys1[j]]['irrevlhsnodes'][idx][idx2],organismsdata[orgkeys1[j]]['irrrxnno'][idx])])
            for idx3 in range(len(organismsdata[orgkeys1[j]]['irrevrhsnodes'][idx])):
                G.add_edges_from([(organismsdata[orgkeys1[j]]['irrrxnno'][idx],organismsdata[orgkeys1[j]]['irrevrhsnodes'][idx][idx3])])
        for idx in range(len(organismsdata[orgkeys1[j]]['revrxnno'])):
             for idx2 in range(len(organismsdata[orgkeys1[j]]['revlhsnodes'][idx])):
                 G.add_edges_from([(organismsdata[orgkeys1[j]]['revlhsnodes'][idx][idx2],organismsdata[orgkeys1[j]]['revrxnno'][idx])])
                 G.add_edges_from([(organismsdata[orgkeys1[j]]['revbacrxno'][idx],organismsdata[orgkeys1[j]]['revlhsnodes'][idx][idx2])])
             for idx3 in range(len(organismsdata[orgkeys1[j]]['revrhsnodes'][idx])):
                 G.add_edges_from([(organismsdata[orgkeys1[j]]['revrxnno'][idx],organismsdata[orgkeys1[j]]['revrhsnodes'][idx][idx3])])
                 G.add_edges_from([(organismsdata[orgkeys1[j]]['revrhsnodes'][idx][idx3],organismsdata[orgkeys1[j]]['revbacrxno'][idx])])
    return G

def create_graph_with_exchange_reactions(G, orgs, namemap):
    orgkeys1 = orgs.keys()
    mexcids = []
    for j in range(len(orgkeys1)):
        ec = orgs[orgkeys1[j]]['exchange']
        mexcids.append(ec)
    commonexc = list(set.intersection(*map(set,mexcids))) #Common exchange metabolites in different organisms
    for j in range(len(orgkeys1)):
        renamedexc = [orgkeys1[j] + ' ' + s for s in commonexc]
        exclen = range(0,len(commonexc))
        excrxnno = ['Org_%s ER' %orgkeys1[j] +  str(t+1) for t in exclen]
        excrxnnor = ['Org_%s ERR' %orgkeys1[j] + str(t+1) for t in exclen]
        G.add_nodes_from(excrxnno, bipartite=1)
        G.add_nodes_from(excrxnnor, bipartite=1)
        G.add_nodes_from(commonexc, bipartite=0)
        G.add_nodes_from(renamedexc, bipartite=0)
        for k in range(len(renamedexc)):
            namemap[excrxnno[k]] = commonexc[k]
            namemap[excrxnnor[k]] = commonexc[k]
            G.add_edges_from([(renamedexc[k],excrxnno[k])])
            G.add_edges_from([(excrxnno[k],commonexc[k])])
            G.add_edges_from([(commonexc[k],excrxnnor[k])])
            G.add_edges_from([(excrxnnor[k],renamedexc[k])])
    for j in range(len(orgkeys1)):
        metitems = orgs[orgkeys1[j]]['exchange']
        noncommonexc = list(set(metitems) - set(commonexc))
        ncrenamedexc = [orgkeys1[j] + ' ' + s for s in noncommonexc]
        ncexclen = range(0,len(noncommonexc))
        ncexcrxnno = ['Org_%s NCER' %orgkeys1[j] +  str(t+1) for t in ncexclen]
        ncexcrxnnor = ['Org_%s NCERR' %orgkeys1[j] + str(t+1) for t in ncexclen]
        G.add_nodes_from(ncexcrxnno, bipartite=1)
        G.add_nodes_from(ncexcrxnnor, bipartite=1)
        G.add_nodes_from(noncommonexc, bipartite=0)
        G.add_nodes_from(ncrenamedexc, bipartite=0)
        for k in range(len(ncrenamedexc)):
            namemap[ncexcrxnno[k]] = noncommonexc[k]
            namemap[ncexcrxnnor[k]] = noncommonexc[k]
            G.add_edges_from([(ncrenamedexc[k],ncexcrxnno[k])])
            G.add_edges_from([(ncexcrxnno[k],noncommonexc[k])])
            G.add_edges_from([(noncommonexc[k],ncexcrxnnor[k])])
            G.add_edges_from([(ncexcrxnnor[k],ncrenamedexc[k])])
    return G, excrxnno, excrxnnor,noncommonexc, ncexclen,ncrenamedexc,commonexc,namemap


print 'Enter full path where the .xml files are'
pname = raw_input()
excrxnno =[]
excrxnnor =[]
organismsdata, allnamemap, modelids = fetchreactions.segregate_reactions_from_models(pname)
if organismsdata:
    print 'How many organisms should be in the community?'
    no_of_orgs = input()
    orgkeys = organismsdata.keys()
    combolist = list(itertools.combinations(range(len(orgkeys)),int(no_of_orgs)))
    if combolist:
        for m in range(len(combolist)):
            fnm = ''
            tempdict = {}
            for n in range(len(combolist[m])):
                tempdict[orgkeys[combolist[m][n]]] = organismsdata[orgkeys[combolist[m][n]]]
                fnm = fnm +  orgkeys[combolist[m][n]] + '_'
            H = create_graph_with_internal_reaction(tempdict)
            G, excnos1, excnos2,noncommonexc,B,C,commonexc,namemap = create_graph_with_exchange_reactions(H, tempdict, allnamemap)
            nx.write_gpickle(H , fnm+ '.gpickle')
            print 'Graph created'
            print 'Number of edges', len(G.edges())
            print 'Number of nodes', len(G.nodes())
            with open(fnm + 'namemap' + '.pickle', 'w') as f:
                dump(namemap, f)
        sys.path.append(pname)
    else:
        print 'Number of organisms for creating a consortium graph is more than the models given'
        print 'Program will now exit'
        sys.exit()
else:
    print "Cannot create graph"

