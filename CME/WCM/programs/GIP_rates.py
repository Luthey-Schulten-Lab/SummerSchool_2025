"""
Author: Zane Thornburg

Modified by Enguang Fu in 2024

calculate the rates for genetic information processes (GIP) in CME

"""
# 

### Rate Constant Calculations #####


from collections import defaultdict, OrderedDict


#########################################################################################
def replicationRate(sim_properties,locusNum):
    

    locusNumtoGeneSeq = sim_properties['locusNumtoGeneSeq']
    geneSeq = locusNumtoGeneSeq[locusNum]

    # K0rep = 0.26e-3    # mM Binding of DNAP with gene
    KDrep = 0.001      # mM Binding of monomer with the replication complex
    kcatrep = 100      # bp/s
    
    countsDic = sim_properties['counts']

    datp = partTomM(countsDic['M_datp_c'][-1], sim_properties) # 0.009 # mM
    dttp = partTomM(countsDic['M_dttp_c'][-1], sim_properties) #0.011 # mM
    dctp = partTomM(countsDic['M_dctp_c'][-1], sim_properties) #0.006 # mM
    dgtp = partTomM(countsDic['M_dgtp_c'][-1], sim_properties) #0.0035 # mM


    DNApol3 = partTomM(35, sim_properties) # 35*countToMiliMol # mM 0.00173 mM Concentration of DNAP; Not explicitly shown in the reactions but in the rate law


    baseCount = defaultdict(int)
    for base in set(geneSeq):
        baseCount[base] = geneSeq.count(base)

    n_tot = sum(list(baseCount.values()))

    subsystem = 'replication'
    if checkMonomer(sim_properties, subsystem):
        return n_tot, 0

    NMono_A = baseCount["A"]

    NMono_C = baseCount["C"]

    NMono_G = baseCount["G"]

    NMono_T = baseCount["T"]

    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_T]
    

    NMonoSum = NMono_A*KDrep/datp + NMono_C*KDrep/dctp + NMono_T*KDrep/dttp + NMono_G*KDrep/dgtp

    k_gene_rep = kcatrep/((KDrep**2/datp/dttp) + NMonoSum + n_tot - 1)

    # k_gene_rep = kcatrep/((1+K0rep/DNApol3)*(KDrep**2/datp/dttp) + NMonoSum + n_tot - 1)

    return n_tot, k_gene_rep
#########################################################################################



#########################################################################################
def TranscriptionRate(sim_properties, locusNum):
    subsystem = 'trsc_poly'

    if checkMonomer(sim_properties, subsystem):
        return 0

    genome = sim_properties['genome']
    countsDic = sim_properties['counts']

    locusTag = 'JCVISYN3A_' + locusNum
    rnasequence = genome[locusTag]['RNAsequence']

    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
        
    PromoterStrength = sim_properties['promoters'][locusNum]
    
    rnaPolKcat = 20 # nt/s

    kcat_mod = min(rnaPolKcat*(PromoterStrength),85)

    kcat_mod = max(10,kcat_mod)

    # Add total number of monomers to parameter dict
    baseMap = OrderedDict({ "A":"M_atp_c", "U":"M_utp_c", "G":"M_gtp_c", "C":"M_ctp_c" })
    # baseMapToMonoP = OrderedDict({ "A":"M_amp_c", "U":"M_ump_c", "G":"M_gmp_c", "C":"M_cmp_c" })

    Mono1 = baseMap[ rnasequence[0] ]
    CMono1 = partTomM(countsDic[Mono1][-1], sim_properties)

    Mono2 = baseMap[ rnasequence[1] ]
    CMono2 = partTomM(countsDic[Mono2][-1], sim_properties)

    n_tot = sum(list(baseCount.values()))

    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["U"]
    
    NMono_C = baseCount["C"]
    
    NMono_G = baseCount["G"]

    atp = partTomM(countsDic['M_atp_c'][-1], sim_properties)
    ctp = partTomM(countsDic['M_ctp_c'][-1], sim_properties)
    gtp = partTomM(countsDic['M_gtp_c'][-1], sim_properties)
    utp = partTomM(countsDic['M_utp_c'][-1], sim_properties)

    rnaPolKd = 0.1 #mM
    NMonoSum = NMono_A*rnaPolKd/atp + NMono_C*rnaPolKd/ctp + NMono_U*rnaPolKd/utp + NMono_G*rnaPolKd/gtp
    
    k_transcription = kcat_mod / ((rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
    
    return k_transcription
#########################################################################################


#########################################################################################
# Define how to calculate translation rate constants as in equation 3 for translation reactions.

def TranslationRate(sim_properties,locusNum, aasequence):
    """
    Called by addGeneticInformationProcess

    Description: calculate the rates of translation reactions based on charged tRNA counts
    """

    # Check if gtp count == 0 
    subsystem = 'tran_poly'

    if checkMonomer(sim_properties, subsystem):
        return 0
    
    
    # Considers amino acids up to the first stop codon.
#     aasequence = aasequence[0:aasequence.find("*")]
    
    # Check that we know all residues used in the sequence
#     if ( set(aasequence) - set(aaMap.keys()) ):
#         raise Exception("Unknown residue(s) in Protein sequence {}".format(set(aasequence) - set(aaMap.keys())) )
    

    aaCostMap = sim_properties['aaCostMap']
    tRNAmap = sim_properties['trna_map']
    
    currenttime_second = sim_properties['time_second'][-1]

    riboKcat = 12  # /s
    riboKd = 0.001 # mM
    kcat_mod = riboKcat #*ribo_num #*0.4

    # genome = sim_properties['genome']
    # countsDic = sim_properties['counts']
    # locusTag = 'JCVISYN3A_' + locusNum
    # aasequence = genome[locusTag]['AAsequence']

    countsDic = sim_properties['counts']

    # Count how many times each aa residue is used
    aaCount = defaultdict(int)
    for aa in set(aasequence):
        aaCount[aa] = aasequence.count(aa)
    

    NStop = aaCount["*"]

    if NStop > 1:
        print("EXTRA STOP CODON: MISTAKE IN TRANSLATION")
    
    # the first two amino acids for each protein
    monomer1 = aasequence[0]
    monomer2 = aasequence[1]
    # if monomer1 == monomer2:
    #     print(monomer1,monomer2)
    NMonoSum = 0

    if currenttime_second == 0:
        # The assumed counts of aa:tRNA at t = 0 second are 160
        monomer1conc = partTomM(160, sim_properties)
        monomer2conc = partTomM(160, sim_properties)

        for aa in aaCount.keys():
            if aa != '*':
                aaNum = aaCount[aa]
                # aaName: ALA, AUG. ...
                aaName = aaCostMap[aa].split('_')[0]
                tRNAlist = tRNAmap[aaName]

                tRNAcounts = len(tRNAlist)*160
                tRNAconc = partTomM(tRNAcounts, sim_properties)

                NMonoSum +=  aaNum*riboKd/tRNAconc

    else:

        for aa in aaCount.keys():
            if aa != '*':
                aaNum = aaCount[aa]
                # aaName: ALA, AUG. ...
                aaName = aaCostMap[aa].split('_')[0]
                tRNAlist = tRNAmap[aaName]
                # one aa can have multiple tRNA carriers
                tRNAcounts = 0
                for tRNA in tRNAlist:
                    tRNAcounts += countsDic[tRNA + '_ch'][-1]
                
                tRNAconc = partTomM(tRNAcounts, sim_properties)
                if tRNAconc == 0:
                    if locusNum == '0001': # Avoid print WARNING repetitively for all 452 proteins
                        print("WARNING: charged tRNA {0} is 0 at time {1}".format(tRNA, currenttime_second))
                    return 0 # return 0
                NMonoSum +=  aaNum*riboKd/tRNAconc

                if aa == monomer1:
                    monomer1conc = tRNAconc
                if aa == monomer2:
                    monomer2conc = tRNAconc
            

    n_tot = sum(list(aaCount.values()))

    
    k_translation = kcat_mod / ((riboKd**2)/(monomer1conc*monomer2conc) + NMonoSum + n_tot - 1)

    
    return k_translation


# def translationRate_polysome(sim_properties, aasequence):
#     """
    
#     Called by addGeneticInformationProcess 

#     Description: calculate the rates of translation reactions based on charged tRNA counts
#                 consider polysome as a prefactor in the kinectic constant

#     """
#     aaCostMap = {"A":"ALA_cost", "R":"ARG_cost", 
#         "N":"ASN_cost", "D":"ASP_cost", "C":"CYS_cost", "E":"GLU_cost", "Q":"GLN_cost", "G":"GLY_cost", 
#             "H":"HIS_cost", "I":"ILE_cost", "L":"LEU_cost", "K":"LYS_cost", "M":"MET_cost", "F":"PHE_cost", 
#         "P":"PRO_cost", "S":"SER_cost", "T":"THR_cost", "W":"TRP_cost", "Y":"TYR_cost", "V":"VAL_cost"}
    
#     tRNAmap = sim_properties['trna_map']
    
#     currenttime_second= sim_properties['time_second'][-1]

#     riboKcat = 12  # /s
#     riboKd = 0.001 # mM
#     kcat_mod = riboKcat #*ribo_num #*0.4
    
#     countsDic = sim_properties['counts']

#     # Count how many times each aa residue is used
#     aaCount = defaultdict(int)
#     for aa in set(aasequence):
#         aaCount[aa] = aasequence.count(aa)
    

#     NStop = aaCount["*"]

#     if NStop > 1:
#         print("EXTRA STOP CODON: MISTAKE IN TRANSLATION")
    
#     # the first two amino acids for each protein
#     monomer1 = aasequence[0]
#     monomer2 = aasequence[1]
#     # if monomer1 == monomer2:
#     #     print(monomer1,monomer2)
#     NMonoSum = 0

#     if currenttime_second == 0:
#         # The assumed counts of aa:tRNA at t = 0 second are 160
#         monomer1conc = partTomM(160, sim_properties)
#         monomer2conc = partTomM(160, sim_properties)

#         for aa in aaCount.keys():
#             if aa != '*':
#                 aaNum = aaCount[aa]
#                 # aaName: ALA, AUG. ...
#                 aaName = aaCostMap[aa].split('_')[0]
#                 tRNAlist = tRNAmap[aaName]

#                 tRNAcounts = len(tRNAlist)*160
#                 tRNAconc = partTomM(tRNAcounts, sim_properties)

#                 NMonoSum +=  aaNum*riboKd/tRNAconc

#     else:

#         for aa in aaCount.keys():
#             if aa != '*':
#                 aaNum = aaCount[aa]
#                 # aaName: ALA, AUG. ...
#                 aaName = aaCostMap[aa].split('_')[0]
#                 tRNAlist = tRNAmap[aaName]
#                 # one aa can have multiple tRNA carriers
#                 tRNAcounts = 0
#                 for tRNA in tRNAlist:
#                     tRNAcounts += countsDic[tRNA + '_ch'][-1]
                
#                 tRNAconc = partTomM(tRNAcounts, sim_properties)

#                 NMonoSum +=  aaNum*riboKd/tRNAconc

#                 if aa == monomer1:
#                     monomer1conc = tRNAconc
#                 if aa == monomer2:
#                     monomer2conc = tRNAconc
            

#     n_tot = sum(list(aaCount.values()))

    
#     k_translation = kcat_mod / ((riboKd**2)/(monomer1conc*monomer2conc) + NMonoSum + n_tot - 1)

    

#     return k_translation

#########################################################################################
def mrnaDegradationRate(sim_properties, locusNum):
    # indenpendent with metabolism

    # Count how many times each base is used
    subsystem = 'deg_depoly'

    if checkMonomer(sim_properties, subsystem):
        return 0
    
    genome = sim_properties['genome']

    locusTag = 'JCVISYN3A_' + locusNum

    rnasequence = genome[locusTag]["RNAsequence"]

    baseCount = defaultdict(int)

    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
    
    kcat = 88 #1/s

    n_tot = sum(list(baseCount.values()))

    k_deg = kcat / n_tot 
    # # The degradation rates of mRNAs vary from 0.016 to 1.1.
    return k_deg
#########################################################################################


#########################################################################################
def TranslocationRate(sim_properties, locusNum):
    # Independent with metabolism
    # Count how many times each residue is used
    subsystem = 'translocation'
    if checkMonomer(sim_properties, subsystem):
        return 0
    
    genome = sim_properties['genome']

    locusTag = 'JCVISYN3A_' + locusNum

    aasequence = genome[locusTag]["AAsequence"]

    aaCount = defaultdict(int)
    for aa in set(aasequence):
        aaCount[aa] = aasequence.count(aa)
    
    ptnLen = sum(list(aaCount.values()))
    
    k_transloc = 50/ptnLen #secyNum*
    
    return k_transloc
#########################################################################################

def checkMonomer(sim_properties, subsytem):
    """
    Return: True or False, True means shortage, False means no shortage

    Check the availabilty of monomers in subsystems including replication, trancription, translation, mRNA degradation, and translocation 
    """
    countsDic = sim_properties['counts']

    atp = countsDic['M_atp_c'][-1]

    if subsytem == 'replication':
        datp = countsDic['M_datp_c'][-1]; dttp = countsDic['M_dttp_c'][-1]
        dctp = countsDic['M_dctp_c'][-1]; dgtp = countsDic['M_dgtp_c'][-1]
        if datp == 0 or dttp == 0 or dctp == 0 or dgtp == 0 or atp == 0:
            return True
        
    elif subsytem == 'trsc_poly':
        atp = countsDic['M_atp_c'][-1]; utp = countsDic['M_utp_c'][-1]
        ctp = countsDic['M_ctp_c'][-1]; gtp = countsDic['M_gtp_c'][-1]
        if atp == 0 or utp == 0 or ctp == 0 or gtp == 0:
            return True
    
    elif subsytem == 'tran_poly':
        gtp = countsDic['M_gtp_c'][-1]
        if gtp == 0:
            return True
        
    elif subsytem == 'translocation':
        if atp == 0:
            return True
    elif subsytem == 'deg_depoly':
        if atp == 0:
            return True


    rxns_prefix = {'initiation':'init', 'replication':'rep', 'trsc_binding':'trsc_binding', 'trsc_poly':'trsc_poly',
                   'tran_binding':'tran_binding', 'tran_poly':'tran_poly', 'translocation':'translocation','deg_binding':'deg_binding', 'deg_depoly': 'deg_depoly', 
                   'ribo_biogenesis':'ribo', 'tRNACharging':'tRNA'}
    
#########################################################################################
def partTomM(particles, sim_properties):
    """
    Convert particle counts to mM concentrations for the ODE Solver

    Parameters:
    particles (int): The number of particles for a given chemical species

    Returns:
    conc (float): The concentration of the chemical species in mM
    """

    ### Constants
    NA = 6.022e23 # Avogadro's

    conc = (particles*1000.0)/(NA*sim_properties['volume_L'][-1])

    return conc
#########################################################################################



