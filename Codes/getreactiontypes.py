# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:38:56 2017

@author: Aarthi Ravikrishnan
"""

import numpy as np

def find_exchange_reactions(x,model,temp,mets,fname):
    xdim,ydim = np.shape(x)
    Amat=[]
    Bmat=[]
    Cmat=[]
    Amat_len=[]
    Bmat_len=[]
    Cmat_len=[]
    ExcIdx=[]
    rxns=[]
    rxnNames=[]
    for j in model.reactions:
        rxns.append(j.reaction)
    for j in range(len(model.reactions)):
        reactions_in_model = model.reactions[j]
        rxnNames.append(reactions_in_model.id)
    for i in range(xdim):
        Amat.append(np.where(x[i]==-1))
        Bmat.append(np.where(x[i]!=0))
        Cmat.append(np.where(x[i]==1))
        Amat_len.append(len(Amat[i][0]))
        Bmat_len.append(len(Bmat[i][0]))
        Cmat_len.append(len(Cmat[i][0]))
        #Case 1 - Presence of bulk metabolites in the medium
        if rxns[i][-1] == 'b': #Assuming the bulk metabolites end in 'b'
            if Amat_len[i]==1 and Cmat_len[i]==1:
                ExcIdx.append(i)
        #Case 2 - Presence of exchange metabolites
        elif Amat_len[i]==1 and Bmat_len[i]==1:
            ExcIdx.append(i)
        elif Cmat_len[i]==1 and Bmat_len[i]==1:
            ExcIdx.append(i)
    tempmets = []
    excnodes = []
    excmetids =[]
    exchids=[]
    excnodestemp =[]
    for l in ExcIdx:
        exchids.append(rxnNames[l])
        if rxns[l][-1] == 'b':
            excnodes.append(mets[np.nonzero(x[l])[0][0]]) # This will generate a list of tuples
        else:
            excmetids.append(np.nonzero(x[l]))
    if excmetids:
        for k in range(len(excmetids)):
            dummy2 =excmetids[k][0].tolist() # To convert it to an array
            excnodestemp.append(dummy2)
        excnodestemp1=[item4 for sublist in excnodestemp for item4 in sublist]
        for m in excnodestemp1:
            excnodes.append(mets[m])
    allrxnids = []
    for i in range(len(rxns)):
        allrxnids.append(i)

    InternalRxns=list(set(allrxnids)^set(ExcIdx))
    ReversibleRxns = []
    IrreversibleRxns = []
    lb = temp.lower_bounds
    ub = temp.upper_bounds
    for i in InternalRxns: #Changed to InternalRxns from range(len(InternalRxns))
        if lb[i] < 0 and ub[i] >= 0:
            ReversibleRxns.append(i)
        elif lb[i] >=0 and ub[i] >=0:
            IrreversibleRxns.append(i)
    #Irreversible part
    irrevlhstemp1 = []
    irrevlhstemp =[]
    irrevlhsnodes=[]
    tempmets =[]
    irrevrxnids = []
    for j in IrreversibleRxns:
        irrevrxnids.append(rxnNames[j])
        irrevlhstemp1.append(np.where(x[j]<0))
    for k in range(len(irrevlhstemp1)):
        dummy =irrevlhstemp1[k][0].tolist()
        irrevlhstemp.append(dummy)
    for l in range(len(irrevlhstemp)):
        for m in irrevlhstemp[l]:
            metch=fname + ' ' + mets[m]
            tempmets.append(metch)
        irrevlhsnodes.append(tempmets)
        tempmets =[]
    tempmets=[]
    irrevrhsnodes=[]
    irrevrhstemp1 = []
    irrevrhstemp =[]
    for j in IrreversibleRxns:
        irrevrhstemp1.append(np.where(x[j]>0)) #Empty are the ones where the stoichiometric matrix itself has nothing
    for k in range(len(irrevrhstemp1)):
        dummy =irrevrhstemp1[k][0].tolist()
        irrevrhstemp.append(dummy)
    for l in range(len(irrevrhstemp)):
        for m in irrevrhstemp[l]:
            metch=fname + ' ' + mets[m]
            tempmets.append(metch)
        irrevrhsnodes.append(tempmets)
        tempmets =[]
    #Reversible part
    revlhsnodes =[]
    revlhstemp1 = []
    revlhstemp =[]
    revrxnids = []
    for j in ReversibleRxns:
        revrxnids.append(rxnNames[j])
        revlhstemp1.append(np.where(x[j]<0))
    for k in range(len(revlhstemp1)):
        dummy =revlhstemp1[k][0].tolist()
        revlhstemp.append(dummy)
    for l in range(len(revlhstemp)):
        for m in revlhstemp[l]:
            metch=fname + ' ' + mets[m]
            tempmets.append(metch)
        revlhsnodes.append(tempmets)
        tempmets =[]
    revrhsnodes =[]
    revrhstemp1 = []
    revrhstemp =[]
    for j in ReversibleRxns:
        revrhstemp1.append(np.where(x[j]>0))
    for k in range(len(revrhstemp1)):
        dummy =revrhstemp1[k][0].tolist()
        revrhstemp.append(dummy)
    for l in range(len(revrhstemp)):
        for m in revrhstemp[l]:
            metch=fname + ' ' + mets[m]
            tempmets.append(metch)
        revrhsnodes.append(tempmets)
        tempmets =[]
    return ReversibleRxns, IrreversibleRxns, ExcIdx, excnodes, irrevlhsnodes, \
        irrevrhsnodes, revlhsnodes, revrhsnodes, exchids, irrevrxnids, revrxnids
