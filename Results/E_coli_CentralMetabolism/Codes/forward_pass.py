# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 09:58:43 2016

This code performs the Phase 1 of MetQuest.

Input to this code:
1. Python NetworkX graph object (DiGraph)
2. set of seed metabolites
3. src - Set of source molecules
4. Counter value
        0 - If the Guided BFS has to be carried out on a 
 	   forward graph
	1 - If the Guided BFS has to be carried out on a 
 	   reverse graph
@author: Aarthi Ravikrishnan
"""

from collections import deque, defaultdict

def forward_pass(G, alwaysexist2, src, countertoadd, *args): 
    targets1 = set([])
    pred = G.predecessors
    succ = G.successors
    alwaysexist = alwaysexist2.copy()
    for ele in args:
        for targets in ele:
            targets1.add(targets)
            if targets in alwaysexist:
                alwaysexist.remove(targets)
    if countertoadd == 1:
        for srcs in src:
            for tarrxns in succ(srcs):
                for tarreqmets in pred(tarrxns):
                    alwaysexist.add(tarreqmets)
    lowerboundmetabolite = defaultdict(list)
    lowerboundreaction = defaultdict(list)

    for item in alwaysexist:
        lowerboundmetabolite[item].append(0)
    stage = 1
    queue = deque([])
    status_dict = defaultdict(str)
    alwaysexist.add(src[0])
    mediacomp = alwaysexist.copy()
    alwaysexist_first = alwaysexist.copy()
    startnode = []
    cannotbetriggered = defaultdict(str)
    for startingnodes in alwaysexist:
        if startingnodes in G:
            for startingrxns in succ(startingnodes):
                if set(pred(startingrxns)).issubset(alwaysexist):
                    startnode.append(startingrxns)
                    for metsprod in succ(startingrxns):
                        mediacomp.add(metsprod)
                        alwaysexist_first.add(metsprod)

                        if stage not in lowerboundmetabolite[metsprod]:
                            lowerboundmetabolite[metsprod].append(stage)
                    if stage not in lowerboundreaction[startingrxns]:
                        lowerboundreaction[startingrxns].append(stage)
    for startrxn in startnode:
        temp = succ(startrxn)
        for item in temp:
            for rxns in succ(item):
                if set(pred(rxns)).issubset(mediacomp):
                    queue.append(rxns)
                else:
                    cannotbetriggered[rxns] = 'Y'
        status_dict[startrxn] = 'V'

    while queue:
        stage += 1
        for parent in list(queue):

            if status_dict[parent] == '':
                temp = succ(parent)
                if stage not in lowerboundreaction[parent]:
                    lowerboundreaction[parent].append(stage)
                for item in temp:
                    mediacomp.add(item)
                    if stage not in lowerboundmetabolite[item]:
                        lowerboundmetabolite[item].append(stage)
                    for progeny in succ(item):
                        if set(pred(progeny)).issubset(mediacomp):
                            if status_dict[progeny] != 'V':
                                queue.append(progeny)
                                status_dict[parent] = 'V'
                        else:
                            cannotbetriggered[progeny] = 'Y'
            elif status_dict[parent] == 'V':
                for item2 in succ(parent):
                    if stage not in lowerboundmetabolite[item2]:
                        lowerboundmetabolite[item2].append(stage)
            queue.popleft()
    return lowerboundmetabolite, lowerboundreaction, status_dict, mediacomp, alwaysexist_first
