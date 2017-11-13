# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 16:24:38 2017

@author: aarthi
"""
from __future__ import division
from itertools import combinations
from matplotlib import pyplot as plt
import numpy as np
import collections
import pickle
import matplotlib.font_manager as fm
import matplotlib as mpl

def find_jaccard_between_paths(onlysrctotar):
    jaccard_values = {}
    flag = 0
    pathsanalysed = {}
    for reactionlists in combinations(onlysrctotar, 2):
        flag+= 1
        j_value = len(set(reactionlists[0]).intersection(set(reactionlists[1])))/len(set(reactionlists[0]).union(set(reactionlists[1])))
        jaccard_values[flag] =j_value
        pathsanalysed[flag] = reactionlists
    return jaccard_values, pathsanalysed

if __name__ == '__main__':
    print 'Load List of lists of subnetworks'
    [jaccard_values, pathsanalysed] = find_jaccard_between_paths(onlysrctotar)
    with open('jaccard-values.pickle','w') as f:
        pickle.dump(jaccard_values,f)
    csfont = {'fontname':'Comic Sans MS'}
    hfont = {'fontname':'Helvetica'}
    paths_diff_len = []
    for item in onlysrctotar:
        paths_diff_len.append(len(item))
    numbers = collections.Counter(paths_diff_len)
mostdiffent = [330148, 30150, 330287, 330289]#335221,335223,335506,335508,338392,338394,338551,338553] #Have least Jaccard value
for i in mostdiffent:
    diffones = pathsanalysed[i]
    for x,item in enumerate(diffones):
        fname = 'Results/' + str(i) + str(x) + '.txt'
        with open(fname, 'w') as g:
            for item2 in item:
                print >> g, namemap[item2],'\t', pred(item2),'\t', succ(item2)

for i in mostdiffent:
    diffones = pathsanalysed[i]
    for x,item in enumerate(diffones):
        fname = 'Results/' + str(i) + '_forgraph' + str(x) + '.txt'
        with open(fname, 'w') as g:
            for item2 in item:
                for predele in pred(item2):
                    if len(predele.split(' ')) > 1:
                        print >> g, predele.split(' ')[1], '\t', namemap[item2]
                    else:
                        print >> g, predele, '\t', namemap[item2]
                for succele in succ(item2):
                    if len(succele.split(' ')) > 1:
                        print >> g, namemap[item2], '\t', succele.split(' ')[1]
                    else:
                        print >> g, namemap[item2], '\t', succele
            print>>g, '\n'

