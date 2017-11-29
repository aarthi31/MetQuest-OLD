# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:38:56 2017

@author: Aarthi Ravikrishnan
"""

from getreactiontypes import find_exchange_reactions
import cobra
import os
import glob
import copy
import pdb
organismsdata = {}
namemap = {}
allmodelids = []
def segregate_reactions_from_models(pathname):
    os.chdir(pathname)
    fnames = glob.glob('*.xml')
    if len(fnames) == 0:
        print "There are no .xml files. Please check the path"

    organismsdata = {}
    print 'Filenames', fnames

    for m in range(len(fnames)):
        #pdb.set_trace()
        model = cobra.io.read_sbml_model(fnames[m])
        #try:
        #    len
        #if len(model) > 1:
        #    [model, badrxns] = cobra.io.read_sbml_model(fnames[m])

        #print model.id
        #print type(model)
        #try:
        #if len(model) > 1:

        if isinstance(model, tuple):
            temp = cobra.core.ArrayBasedModel(model[0])
            if model[0].id:
                org = model[0].id

            else:
                print "Model ID not found; using file name instead"
                org =  fnames[m].split('.')[0]
            model = model[0]

        else:
            temp = cobra.core.ArrayBasedModel(model)
            if model.id:
                org = model.id
            else:
                print "Model ID not found; using file name instead"
                org =  fnames[m].split('.')[0]

        allmodelids.append(org)
        orgclasdata={org:{'exchange':[],'irrevlhsnodes':[],'irrevrhsnodes':[],'revrhsnodes':[],'revlhsnodes':[],'irrrxnno':[],'revrxnno':[],'exchno':[],'totalnodes':[],'modelrxns':[],'metabolites':[],'exchname':[],'irrevname':[],'revname':[]}}
        metstmp= []
        rxnids =[]
        mets =[]
        rxnstmp =[]
        #pdb.set_trace()
        for metnames in model.metabolites:
            metstmp.append(metnames.id)
            #print metstmp
        for i in range(len(metstmp)):
            mets.append(metstmp[i])
        for rxns in model.reactions:
            rxnstmp.append(rxns)
        for i in range(len(rxnstmp)):
            rxnids.append(rxnstmp[i])
        stoi = temp.S
        stoi_matrix = stoi.toarray()
        stoi_matrix = stoi.T
        #stoi_copy = copy.deepcopy(stoi_matrix)
        x = stoi_matrix.toarray()
        revrxn, irrevrxn, excrxn, exc1, irr1, irr2, rev1, rev2, excname, irrevname, revname = find_exchange_reactions(x,model,temp,mets, org)
        orgclasdata[org]['exchange'] = exc1
        orgclasdata[org]['irrevlhsnodes'] = irr1
        orgclasdata[org]['irrevrhsnodes'] = irr2
        orgclasdata[org]['revlhsnodes'] = rev1
        orgclasdata[org]['revrhsnodes'] = rev2
        orgclasdata[org]['exchname'] = excname
        orgclasdata[org]['irrevname'] = irrevname
        orgclasdata[org]['revname'] = revname
        irrevrxno = []
        for num in range(len(irr1)):
            rno='Org_%s IR' %org + str(num+1)
            irrevrxno.append(rno)
            namemap[rno]= irrevname[num]
        orgclasdata[org]['irrrxnno'] = irrevrxno
        revrxno = []
        for num in range(len(rev1)):
            rno='Org_%s RR' %org + str(num+1)
            revrxno.append(rno)
            namemap[rno]= revname[num]
        orgclasdata[org]['revrxnno'] = revrxno
        revbacrxno = []
        for num in range(len(rev1)):
            rno='Org_%s RevBR' %org + str(num+1)
            revbacrxno.append(rno)
            namemap[rno]= revname[num]
        orgclasdata[org]['totalnodes'] = len(exc1)+len(irr1)+len(rev1)
        orgclasdata[org]['modelrxns']=rxnids
        orgclasdata[org]['revbacrxno'] = revbacrxno
        orgclasdata[org]['metabolites'] = mets
        organismsdata.update(orgclasdata)
    return organismsdata, namemap, allmodelids
