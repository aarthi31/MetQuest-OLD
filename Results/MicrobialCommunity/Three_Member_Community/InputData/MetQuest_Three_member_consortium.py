# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 10:18:39 2016
This code consists of 2 functions. The first function will report the 
sub-networks on a 3-member-consortium. The seconf function will report 
the metabolic exchanges happening between the organisms and the 
metabolites for which acetate and ethanol are being exchanged between
Clostridium cellulolyticum and other organisms.
@author: aarthi
"""


from __future__ import division
from networkx import read_gpickle, get_node_attributes
import cPickle
import itertools
from generate_partitions import generate_partitions
from operator import mul
import os
from collections import Counter,defaultdict
from forward_pass import forward_pass
import gzip
import cPickle
import math

def calculate_biomass_components(G, src, tar1, alwaysexist1, cutoff, filenames):
    succ = G.successors
    pred = G.predecessors
    nattr = get_node_attributes(G, 'bipartite')
    nattr_inv = {}
    for k, v in nattr.iteritems():
        nattr_inv[v] = nattr_inv.get(v, [])
        nattr_inv[v].append(k)  #0 metabolites , 1 are reactions
    for rxns1 in nattr_inv[1]:
        if len(set(pred(rxns1)) - alwaysexist1) >= 5:
            G.remove_node(rxns1)
    #pdb.set_trace()
    if tar1 not in alwaysexist1:
        print tar1
        tar= [tar1]
        lowerbound, lowerboundreaction, status_dict, me,alwaysexist_first  = forward_pass(G, alwaysexist1, src, 0)
        dagsfound2 = {}
        for item in list(alwaysexist1):
            dagsfound2[item] = {0: ''}
        for rxns in status_dict:
            if set(pred(rxns)).issubset(alwaysexist1):
                for succmets in succ(rxns):
                    if succmets not in dagsfound2 and succmets not in alwaysexist1:
                        dagsfound2[succmets] = {1:[]}
                for succmets in succ(rxns):
                    if succmets not in alwaysexist1:
                        dagsfound2[succmets][1].append(set([rxns]))
        for i in range(2, cutoff+1):
            #pdb.set_trace()
            for rxns in status_dict: #Can you parallalize this?
                #if i == 25 and 'S1 pro_DASH_L_c' in succ(rxns):
                 #   pdb.set_trace()
                if set(pred(rxns)).issubset(dagsfound2):

                    metsrequired = list(set(pred(rxns)) - alwaysexist1) # This variable avoids the generation of partitions pertaining to seed metabolites
                    shortestpathlist = []
                    if metsrequired:# and goahead == 'Yes':
                        for predmets in metsrequired:
                            shortestpathlist.append(min(lowerbound[predmets]))
        ##==============================================================================
        ## First loop
        ##==============================================================================
                        for val in range(i-1, len(metsrequired)*(i-2)+1): #range does not include the end value
                            temp = int(math.floor(val/((i-1))))
                            for idx in range(1, temp+1):
                                for combivalue in itertools.combinations(metsrequired, idx):
                                    variablemets = list(set(metsrequired) - set(combivalue))
                                    temprxnlist = []
                                    templen = {}

                                    for metabolites in list(combivalue):
                                        if metabolites in dagsfound2:
                                            if i-1 in dagsfound2[metabolites]:
                                                temprxnlist.append([[rxns]])
                                                templen[metabolites]=len(dagsfound2[metabolites][i-1])
                                        else:
                                            break
                                    flag = ''
                                    if templen:
                                        for mets in list(variablemets):
                                            templen[mets] = len(dagsfound2[mets])
                                        if all([reduce(mul,templen.values()) > 500, set(succ(rxns)).issubset(set(dagsfound2)), set(succ(rxns)).isdisjoint(set(tar))]):
                                            flag = 'NA'
                                            #print templen, 'Loop1'
                                        else:
                                            for mets in combivalue:
                                                if i-1 in dagsfound2[mets]:
                                                    temprxnlist.append(dagsfound2[mets][i-1])
                                                else:
                                                    flag = 'NA'
                                                    break
                                        if flag != 'NA':
                                            tempshortestpath = []
                                            for varmet in list(variablemets):
                                                tempshortestpath.append(min(lowerbound[varmet]))

                                            par = generate_partitions(i, tempshortestpath, val-((i-1)*idx))

                                            for partitions in par:
                                                #print partitions,mets
                                                #num += 1
                                                counter = 0
                                                for varmetidx in range(len(variablemets)):
                                                    if variablemets[varmetidx] in dagsfound2:
                                                        if partitions[varmetidx] in dagsfound2[variablemets[varmetidx]]:
                                                            counter += 1
                                                            templen[variablemets[varmetidx]] = len(dagsfound2[variablemets[varmetidx]][partitions[varmetidx]])
                                                if counter == len(variablemets): #because we ae appending the value to templen, we need another flag to check if the other metabolite is present or not
                                                    if all([reduce(mul,templen.values()) > 1000, set(succ(rxns)).issubset(set(dagsfound2)), set(succ(rxns)).isdisjoint(set(tar))]):
                                                        break #or pass
                                                    else:
                                                        temprxnlist1 = temprxnlist[:]
                                                        for varmetidx in range(len(variablemets)):
                                                            temprxnlist1.append(dagsfound2[variablemets[varmetidx]][partitions[varmetidx]])
                                                        for prodvec in itertools.product(*temprxnlist1):
                                                            rxncomb = set([])
                                                            reactants = set([]) #This is to check for cycles
                                                            for rxnentry in prodvec:
                                                                for individualele in rxnentry:
                                                                    if individualele not in rxncomb:
                                                                        rxncomb.add(individualele)
                                                                        for metsreq in pred(individualele):
                                                                            reactants.add(metsreq)
                                                            for succmets in succ(rxns):
                                                                if succmets not in alwaysexist1 and succmets not in reactants:
                                                                    if succmets in dagsfound2:
                                                                        if len(rxncomb) in dagsfound2[succmets]:
                                                                            try:
                                                                                if dagsfound2[succmets][len(rxncomb)].index(rxncomb):
                                                                                    pass
                                                                            except:
                                                                                    dagsfound2[succmets][len(rxncomb)].append(rxncomb)
                                                                        else:
                                                                            dagsfound2[succmets].update({len(rxncomb): [rxncomb]})
                                                                    else:
                                                                        dagsfound2[succmets] = {len(rxncomb): [rxncomb]}
        #==============================================================================
        # Second loop
        #==============================================================================
                        for val1 in range(len(metsrequired)*(i-2)+1, len(metsrequired)*(i-1)+1):
                            par1 = generate_partitions(i, shortestpathlist,val1)
                            for partitions in par1:
                                    #print partitions,metsrequired
                                    temprxnlist2 = []
                                    templen = {}
                                    flag1 = ''
                                    for item in range(len(metsrequired)):
                                        if metsrequired[item] in dagsfound2:
                                            if partitions[item] in dagsfound2[metsrequired[item]]:
                                                templen[metsrequired[item]]=len(dagsfound2[metsrequired[item]][partitions[item]])
                                    if templen:
                                        if all([reduce(mul,templen.values()) > 5000, set(succ(rxns)).issubset(set(dagsfound2)), set(succ(rxns)).isdisjoint(set(tar))]):
                                            flag1 = 'NA'
                                            #print templen
                                        else:
                                            for item in range(len(metsrequired)):
                                                if partitions[item] in dagsfound2[metsrequired[item]]:
                                                    temprxnlist2.append(dagsfound2[metsrequired[item]][partitions[item]])
                                                else:
                                                    flag1 = 'NA'
                                        if flag1!='NA':
                                            temprxnlist2.append([[rxns]])
                                            for prodvec in itertools.product(*temprxnlist2):
                                                rxncomb1 = set([])
                                                reactants1 = set([]) #This is to check for cycles
                                                for rxnentry in prodvec:
                                                    for individualele in rxnentry:
                                                        if individualele not in rxncomb1:
                                                            rxncomb1.add(individualele)
                                                            for metsreq in pred(individualele):
                                                                reactants1.add(metsreq)
                                                for succmets in succ(rxns):
                                                    if succmets not in alwaysexist1 and succmets not in reactants1:
                                                        if succmets in dagsfound2:
                                                            if len(rxncomb1) in dagsfound2[succmets]:
                                                                    try:
                                                                        if dagsfound2[succmets][len(rxncomb1)].index(rxncomb1):
                                                                            pass
                                                                    except:
                                                                        dagsfound2[succmets][len(rxncomb1)].append(rxncomb1)
                                                            else:
                                                                dagsfound2[succmets].update({len(rxncomb1): [rxncomb1]})
                                                        else:
                                                            dagsfound2[succmets] = {len(rxncomb1): [rxncomb1]}

            print 'Column value', i, 'Cumulative number of metabolites', len(dagsfound2.keys())
        return lowerbound, lowerboundreaction, dagsfound2

def find_exc_in_paths(dagsfound):
    exc_ethanol_rxnnames = []
    exc_ethanol_mets = defaultdict(list)
    exc_acetate_rxnnames = []
    exc_acetate_mets = defaultdict(list)

    for item in dagsfound:
            for pathlen in dagsfound[item]:
                for paths in dagsfound[item][pathlen]:
                    for rxnnames in paths:
                        if 'Org_BMID000000140414 ER231'  in rxnnames: #Ethanol exchange from Clostridium cellulolyticum
                            exc_ethanol_rxnnames.append(rxnnames)
                            exc_ethanol_mets[item].append(rxnnames)
                        elif 'Org_BMID000000140414 ER477' in rxnnames: #Acetate exchange from Clostridium cellulolyticum
                            exc_acetate_rxnnames.append(rxnnames)
                            exc_acetate_mets[item].append(rxnnames)


    return exc_ethanol_rxnnames, exc_ethanol_mets,exc_acetate_rxnnames, exc_acetate_mets
if __name__ == '__main__':
    gname1 = ['/home/aarthi/Dropbox/Consortia_Manuscript/ConsortiaExamples/Simulations/CommunitySImulations/Miller2010Examples/THree/BMID000000142315_BMID000000140414_BMID000000141529.gpickle']
    cutoff1 = 20
    G = read_gpickle(gname1[0])
    filenames = gname1[0].split('/')[-1].split('.')[0].split('_')
    if not os.path.exists('_'.join(filenames)):
        os.mkdir('_'.join(filenames))
    src = ['MNXM776_i']#'MNXM99_i','MNXM105_i']
    alwaysexist1 = set([])
    alwaysexistdec = ['bigg_h_i', 'bigg_h2o_i', 'bigg_atp_i', 'bigg_pi_i',
                      'bigg_nad_i', 'bigg_nadh_i', 'bigg_nadp_i', 'bigg_co2_i',
                      'bigg_k_i', 'bigg_na1_i','bigg_nadph_i', 'bigg_adp_i',
                      'bigg_ppi_i', 'bigg_co2_i', 'bigg_coa_i', 'bigg_amp_i',
                      'bigg_amp_bm', 'bigg_nh3_i', 'bigg_o2_i', 'bigg_glu_L_i',
                      'bigg_gtp_i', 'bigg_ump_i', 'bigg_udp_i',
                      'bigg_gdp_i', 'bigg_cmp_i', 'bigg_ctp_i', 'bigg_utp_i',
                      'bigg_fad_i', 'bigg_gmp_i', 'bigg_dna_i', 'bigg_cdp_i',
                      'bigg_so3_i', 'bigg_dgtp_i', 'MNXM96063_i', 'bigg_no2_i', 'bigg_h2s_i',
                      'bigg_imp_i', 'bigg_dadp_i', 'bigg_dgdp_i', 'bigg_itp_i',
                      'bigg_dcdp_i', 'bigg_datp_i', 'bigg_dctp_i', 'bigg_dtmp_i','bigg_dtmp_bm',
                      'bigg_dutp_i', 'bigg_dudp_i', 'bigg_cl_i', 'bigg_dgmp_i','bigg_gmp_bm',
                      'bigg_idp_i', 'bigg_dcmp_i', 'bigg_dump_i', 'bigg_so4_i','bigg_dgmp_bm',
                      'bigg_pppi_i', 'bigg_damp_i', 'bigg_fmnh2_i', 'MNXM537_i',
                      'bigg_fe2_i', 'bigg_h2_i', 'bigg_ditp_i', 'bigg_ca2_i',
                      'bigg_cobalt2_i', 'metacyc_NAD_P_H_i', 'MNXM24_i',
                      'MNXM90227_i', 'bigg_fadh2_i', 'bigg_nac_i',
                      'bigg_ncam_i', 'bigg_damp_bm','bigg_dtmp_bm',
                      'bigg_ump_bm', 'bigg_dcmp_bm', 'bigg_atp_bm', 'metacyc_NAD_P__i',
                      'bigg_glutrna_i', 'bigg_trnalys_i', 'bigg_asntrna_i', 'bigg_asptrna_i',
                      'bigg_trnamet_i', 'bigg_mettrna_i', 'bigg_trnaglu_i',
                      'bigg_trnaasn_i', 'bigg_glntrna_i', 'bigg_sertrna_sec__i',
                      'bigg_fmettrna_i', 'bigg_leutrna_i',
                      'bigg_trnagln_i', 'bigg_argtrna_i',
                      'bigg_trnaarg_i', 'bigg_alatrna_i', 'bigg_trnaala_i',
                      'bigg_cystrna_i', 'bigg_trnacys_i', 'bigg_glytrna_i',
                      'bigg_trnagly_i', 'bigg_sertrna_i', 'bigg_trnaser_i',
                      'bigg_trnatyr_i', 'bigg_protrna_i', 'bigg_trnapro_i',
                      'bigg_trnatrp_i', 'bigg_trnaval_i', 'bigg_trptrna_i',
                      'bigg_tyrtrna_i', 'bigg_valtrna_i', 'bigg_phetrna_i',
                      'bigg_trnaphe_i', 'bigg_trnaasp_i',
                      'bigg_iletrna_i', 'bigg_trnaile_i',
                      'bigg_trnaleu_i', 'bigg_thrtrna_i', 'bigg_trnathr_i',
                      'bigg_histrna_i', 'bigg_trnahis_i', 'bigg_lystrna_i',
                      'bigg_sectrna_i','bigg_NPmehis_i','MNXM2623_i',
                      'MNXM11736_i','MNXM11739_i','MNXM96096_i','bigg_argtrna_i',
                      'MNXM101_i','bigg_3mop_i','MNXM24_i']
    for alwaysitem in alwaysexistdec:
        for modelnames in filenames:
            alwaysexist1.add(modelnames +' ' +alwaysitem)
    alwaysexist1.add('MNXM776_i')
    lowerbound, lowerboundreaction, status_dictfp, me,ae1 = forward_pass(G, alwaysexist1, src, 0)
    taritems = ''
    lb, lbr, dagsfound = calculate_biomass_components(G, src, taritems, alwaysexist1, cutoff1, filenames)
    [exc_ethanol_rxnnames, exc_ethanol_mets,exc_acetate_rxnnames, exc_acetate_mets] = find_exc_in_paths(dagsfound)
