# -*- coding: utf-8 -*-
"""

This program will check if the prerequisites for running MetQuest are met.
If not, the program will automatically exist.
This is an interactive program, input to this is the pathname of the folder
which contains the sbml files of the organisms, the source, seed and the
target metabolites in a separate text file named as "source_mets.txt",
"seed_mets.txt" and "target_mets.txt". All the entries in this file has to
be separated by newline.

This file can be run on command line as
python exec_metquest.py
Follow the instructions therein.

@author: aarthi
"""
from __future__ import division
from sys import exit
try:
    from networkx import read_gpickle, get_node_attributes
except ImportError:
    print "Networkx not found. Please install networkx package"
    print "Program will now exit"
    exit()

try:
    from cPickle import load
except ImportError:
    print "cPickle not found. Please install cPickle package"
    print "Program will now exit"
    exit()

try:
    import cobra
except ImportError:
    print "cobrapy not found. Please install cobrapy package"
    print "Program will now exit"
    exit()

try:
    from collections import deque, defaultdict
except ImportError:
    print "collections not found. Please install collections package"
    print "Program will now exit"
    exit()

import os
import metquest
import construct_graph

print "Enter the path length"
cutoff = int(input())
for files in os.listdir(construct_graph.pname):
    if files.endswith(".txt"):
        if files.startswith("seed"):
            with open(files, 'r') as f:
                seedmetslist = f.read().splitlines()
            seedmetabolites = set(seedmetslist)
        elif files.startswith("source"):
            with open(files, 'r') as f:
                src = f.read().splitlines()
        elif files.startswith("target"):
            with open(files, 'r') as f:
                targetmetabolites = f.read().splitlines()
alwaysexist2 = seedmetabolites.copy()
for compounds in alwaysexist2:
    if compounds not in construct_graph.G:
        print compounds, "not in G"
        seedmetabolites.remove(compounds)
for sourcele in src:
    if sourcele not in construct_graph.G:
        sourcenotfound = 'Y'
        print sourcele, "not in G, traversal will use only the seed metabolites"
    else:
        seedmetabolites.add(sourcele)

print "Total length of seed metabolite set and the source is" , len(seedmetabolites)

for targetele in targetmetabolites:
    if targetele not in construct_graph.G:
        print targetele, "not in G, MetQuest will only find the pathways to all the metabolites in the scope of seed"
        tarnotfound = 'Y'
filenames = '_'.join(construct_graph.modelids)
succ = construct_graph.G.successors
pred = construct_graph.G.predecessors
lowerbound, lowerboundreaction, pathway_table, timetaken, scope, status_dict = metquest.find_pathways(construct_graph.G, src, targetmetabolites, seedmetabolites, cutoff)
print 'Time taken', timetaken
if not os.path.exists('Results'):
    os.makedirs('Results')
for itemtar in targetmetabolites:
    if itemtar not in pathway_table:
        print "Target could not be found. Consider changing the cut-off or the seed metabolite set"
    else:
        if cutoff not in pathway_table[itemtar]:
            print "No pathways of length", str(cutoff), "from seed to target"
        else:
            for sourcemets in src:
                if sourcemets not in construct_graph.G:
                    print "Source not found in graph"
                    seedpaths = raw_input("Do you want to print the paths arising from seed metabolites? Y or N \n")
                    if seedpaths == 'Y':
                        fname1 = 'Results/' + itemtar + '_allpathwaysfromseed' + '.txt'
                        pathnumcount = 0
                        with open(fname1,'w') as f:
                            for idx in pathway_table[itemtar]:
                                print>>f, idx, '\n'
                                for items in pathway_table[itemtar][idx]:
                                    pathnumcount += 1
                                    print>>f, pathnumcount
                                    for entities in list(items):
                                        print>>f, construct_graph.namemap[entities],\
                                           '\t', ' + '.join(pred(entities)),\
                                           '->', ' + '.join(succ(entities))
                                    print>>f, '\n'
                else:
                    fname1 = 'Results/' + itemtar + '_allpathwaysfromsourceandseed' + '.txt'
                    fname2 = 'Results/' + itemtar + '_allpathwaysfromsource' + '.txt'

                    pathnumcount = 0
                    print_flag = raw_input("Do you want to print all pathways of length less than or equal to pathways? Y or N \n")
                    if print_flag == 'Y':
                        onlysrctotar = []
                        for sourcemets in src:
                            for idx in pathway_table[itemtar]:
                                for items in pathway_table[itemtar][idx]:
                                    if set(succ(sourcemets)).intersection(items):
                                        onlysrctotar.append(list(items))
                        with open(fname1,'w') as f:
                            for idx in pathway_table[itemtar]:
                                print>>f, idx, '\n'
                                for items in pathway_table[itemtar][idx]:
                                    pathnumcount += 1
                                    print>>f, pathnumcount
                                    for entities in list(items):
                                        print>>f, construct_graph.namemap[entities], \
                                       '\t', ' + '.join(pred(entities)), '->',\
                                      ' + '.join(succ(entities))
                                    print>>f, '\n'
                        with open(fname2, 'w') as f:
                            for x,listentries in enumerate(onlysrctotar):
                                print>>f, x+1
                                for entities in listentries:
                                        print>>f, construct_graph.namemap[entities], \
                                       '\t', ' + '.join(pred(entities)), '->',\
                                      ' + '.join(succ(entities))

                                print>>f,'\n'
                    else:
                        onlysrctotar = []
                        for sourcemets in src:
                            for items in pathway_table[itemtar][cutoff]:
                                if set(succ(sourcemets)).intersection(items):
                                    onlysrctotar.append(list(items))

                        with open(fname1, 'w') as f:
                            for items in pathway_table[itemtar][cutoff]:
                                    pathnumcount += 1
                                    print>>f, pathnumcount
                                    for entities in list(items):
                                        print>>f, construct_graph.namemap[entities], \
                                               '\t', ' + '.join(pred(entities)), '->',\
                                                ' + '.join(succ(entities))
                                    print>>f, '\n'

                        with open(fname2, 'w') as f:
                            for x,listentries in enumerate(onlysrctotar):
                                print>>f, x+1
                                for entities in listentries:
                                        print>>f, construct_graph.namemap[entities], \
                                       '\t', ' + '.join(pred(entities)), '->',\
                                      ' + '.join(succ(entities))
                                print>>f,'\n'
