# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 21:26:46 2016

This code takes as input the columnvalue (j), values of the shortest path of
each of the metabolites (given as a list) and the sum that has to be obtained
using these combination of numbers

For instance, if the column value is 7, the number of imputs is 2,
and the shortest path of the metabolites is 4,3 respectively, and the maximum
sum that has to be obtained is 8, then

>>> generate_partitions(7,[4,3],8)
[(4, 4), (5, 3)]

>>> generate_partitions(4, [2,1,1], 5)
[(2, 1, 2), (2, 2, 1), (3, 1, 1)]

@author: Aarthi Ravikrishnan
"""
from itertools import product

def generate_partitions(columnvalue,shortestpathvalue, checkvalue):
    allcombinations= []
    for item in shortestpathvalue:
        allcombinations.append(range(item, columnvalue))
    #print allcombinations
    allpartitions = []
    for partitions in product(*allcombinations):
        if sum(partitions) == checkvalue:
            allpartitions.append(partitions)
    return allpartitions




