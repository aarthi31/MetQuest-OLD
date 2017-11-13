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
    jaccard_values = []
    for reactionlists in combinations(onlysrctotar, 2):
        j_value = len(set(reactionlists[0]).intersection(set(reactionlists[1])))/len(set(reactionlists[0]).union(set(reactionlists[1])))
        jaccard_values.append(j_value)
    return jaccard_values

if __name__ == '__main__':
    path = '/home/aarthi/Dropbox/Fonts/Helvetica.ttf'
    prop = fm.FontProperties(fname=path)
    #prop.set_weight = 'light'
    mpl.rcParams['font.family'] = prop.get_name()
    mpl.rcParams['font.weight'] = 'light'

    print 'Load List of lists of subnetworks'
    jaccard_values = find_jaccard_between_paths(onlysrctotar)
    with open('jaccard-values.pickle','w') as f:
        pickle.dump(jaccard_values,f)
    csfont = {'fontname':'Comic Sans MS'}
    hfont = {'fontname':'Helvetica'}
    prop = fm.FontProperties(fname=path)
    plt.rc('xtick',labelsize=36)
    plt.rc('ytick',labelsize=36)
    plt.hist(jaccard_values, color = 'black', bins=np.arange(min(jaccard_values), max(jaccard_values) +0.025,0.025), alpha = 0.5)
    plt.xlabel("Jaccard values", fontsize=36,fontproperties=prop)
    plt.ylabel("Count",fontsize=36,fontproperties=prop)
    #plt.title("Histogram of jaccard values of the paths",fontsize=14,fontproperties=prop)
  #  plt.savefig("Histogram-jaccard.jpeg")
    paths_diff_len = []
    for item in onlysrctotar:
        paths_diff_len.append(len(item))
    numbers = collections.Counter(paths_diff_len)
    plt.bar(numbers.keys(), numbers.values())
