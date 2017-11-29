# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:46:52 2017

This code carries out the Guided Breadth First Search as explained in 
the main manuscript.

@author: Aarthi Ravikrishnan
"""


def forward_pass(G, alwaysexist2, src, countertoadd, *args): #tar,
    from collections import deque, defaultdict
    global deque
    global defaultdict
    targets1 = set([])
    pred = G.predecessors
    succ = G.successors
    alwaysexist = alwaysexist2.copy()
    for sourcemets in src:
        alwaysexist.add(sourcemets)
    for ele in args:
        for targets in ele:
            targets1.add(targets)
            if targets in alwaysexist:
                alwaysexist.remove(targets)
    if countertoadd == 1:
        for srcs in src:
            if srcs in G:
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

    scope = alwaysexist.copy()
    alwaysexist_first = alwaysexist.copy()
    startnode = []
    cannotbetriggered = defaultdict(str)
    for startingnodes in alwaysexist:
        if startingnodes in G:
            for startingrxns in succ(startingnodes):
                if set(pred(startingrxns)).issubset(alwaysexist):
                    startnode.append(startingrxns)
                    for metsprod in succ(startingrxns):
                        scope.add(metsprod)
                        alwaysexist_first.add(metsprod)

                        if stage not in lowerboundmetabolite[metsprod]:
                            lowerboundmetabolite[metsprod].append(stage)
                    if stage not in lowerboundreaction[startingrxns]:
                        lowerboundreaction[startingrxns].append(stage)
    for startrxn in startnode:
        temp = succ(startrxn)
        for item in temp:
            for rxns in succ(item):
                if set(pred(rxns)).issubset(scope):
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
                    scope.add(item)
                    if stage not in lowerboundmetabolite[item]:
                        lowerboundmetabolite[item].append(stage)
                    for progeny in succ(item):
                        if set(pred(progeny)).issubset(scope):# and set(pred(progeny)).isdisjoint(targets1):
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
    return lowerboundmetabolite, lowerboundreaction, status_dict, scope
