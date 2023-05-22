'''

 3DS-soil-model, a 3D Soil model adapted from STICS : plant functions
 ******************************
 Authors: G. Louarn 
 

 

'''
from scipy import *
import numpy as np
from copy import deepcopy
from soil3ds.miscel_functions import * #soil3ds miscellaneous soil functions


########## diverse plant fonctions - soil water balance


def Transpi_NC(Et0, ls_epsi, ls_FTSW, leafAlbedo=0.15, FTSWThreshold=0.4):
    """
    
    """
    ls_transp = []
    for i in range(len(ls_epsi)):

        ##Riou
        #coeffTV = ls_epsi[i]/ (1-leafAlbedo)   
        #potentialVineTranspiration = coeffTV * Et0

        ##adaptation eq 7.1 p127 STICS et 7.8 p131 book
        k_eq = 0.7 #LAI equivalent pour un k=0.7; pause k car les deux k et LAI sont inconnu
        LAIeq = -log(1-ls_epsi[i])/k_eq 
        potentialTranspiration = Et0 * (1 - exp(-(k_eq-0.2)*LAIeq))#complement a potentialSoilEvaporation
        #rq: la bonne variable a recuperer a la place de espi serait le taux de couverture ou le transmis vertical pour lequel le correctif 0.2 n'est pas necessaire (Eq. 7.2)
        # faut introduire crop coefficient KMAXp cf eq. 7.7 -> pour le moment limite a Et0 (reference gazon)

      
        if (ls_FTSW[i] > FTSWThreshold):#previousTSW/TTSW
            # la transpiration de la plante est faite a son potentielle
            Ks = 1.
        else:
            # regulation
            # la quantite d'eau presente dans le sol n'est pas suffisante pour
            # que la transpiration de plante se fasse a son maximum   
            Ks = ls_FTSW[i]/FTSWThreshold

        ls_transp.append(Ks*potentialTranspiration)

    return ls_transp

#Transpi_NC(2., [0.4], [0.8], leafAlbedo=0.15, FTSWThreshold=0.4)







########## diverse plant fonctions - root distribution



## gestion des racines
##rq: density above 0.5 cm.cm-3 are not taken into account for absorption of water and N in STICS (p90)
##rq2: valeur de 0.1 cm.cm-3 pour definir profondeur d'enracinement dans L2SU(p138)

def vert_roots(dxyz, lvopt):
    """
    
    """
    #""" pour generer un profil vertical de densite racinaire a partir d'une liste - (homogene en x,y)"""
    
    m_roots = []
    for z in range(len(dxyz[2])):
        v = []
        for x in range(len(dxyz[0])):
            vv = []
            for y in range(len(dxyz[1])):
                vv.append(lvopt[z])

            v.append(vv)
        m_roots.append(v)

    return array(m_roots)

#R1 = vert_roots(dxyz, [0.5]*nstrate[2])
#R2 = vert_roots(dxyz, [0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05])


def effective_root_lengths(ls_roots, tresh = 0.5):
    """
    
    """
    #""" for multiple root systems in competition
    #H0 : perfect symetry locally (in a voxel) for ressource aquisition
    #treshold of effective density similar to single species
    #fraction of effective density to each species/plant = proportion of total root length density"""
    
    m = deepcopy(ls_roots)
    nb_r = len(ls_roots)
    tot_root_dens = sum_mat(ls_roots) #sum bug
    for rt in range(nb_r):
        #frac_root = divide(ls_roots[rt], tot_root_dens)## rq: gere bien les voxels vides
        m[rt] = where(tot_root_dens>tresh, tresh*ls_roots[rt]/tot_root_dens, ls_roots[rt])
        
    return m

#def effective_root_length(m_root, tresh = 0.5):#faire avec une liste de root systems ls_root_syst pour competition
#    """ for a single root system """
#    m = deepcopy(m_root)
#    for z in range (len(m_root)):
#        for x in range (len(m_root[z])):
#            for y in range (len(m_root[z][x])):
#                if m_root[z][x][y]>tresh:
#                    m[z][x][y] = tresh
#                else:
#                    m[z][x][y] = m_root[z][x][y]
#    return m

#effective_root_length(R2, 0.5)

def build_ls_roots(RLprofil, S):
    """
    
    """
    #""" version 1 root system :a modifier pour prendre en compte une liste de RLprofil"""
    
    idz = list(RLprofil.keys())
    idz.sort()
    RLprofil_ls = []
    for i in idz: RLprofil_ls.append(RLprofil[i])

    R1 = vert_roots(S.dxyz, RLprofil_ls)*100. /  (S.m_soil_vol*100.*100*100.)#profil densite racinaire en cm.cm-3
    ls_roots = [R1]#[R1, R2]#[R3]#
    return ls_roots

#RLprofil = {0: 0.12450386886407872, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0}#sortie l-system
#RLprofil = [{0: 0.12450386886407872, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0}, {0: 0.145, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0}...]#sortie l-system

def build_ls_roots_mult(RLprofil, S):
    """
    
    """
    #""" version 2: prends en compte une liste de RLprofil"""
    
    nbp = len(RLprofil)
    idz = list(RLprofil[0].keys())
    idz.sort()
    ls_roots = []
    for p in range(nbp):
        RLprofil_ls = []
        for i in idz: RLprofil_ls.append(RLprofil[p][i])
        R = vert_roots(S.dxyz, RLprofil_ls)*100. /  (S.m_soil_vol*100.*100*100.)#profil densite racinaire en cm.cm-3
        ls_roots.append(R)# = [R1]#[R1, R2]#[R3]#

    return ls_roots


def RLprof_t(t, ncouches_sol):
    """
    
    """
    #""" pour generer un profil de densite qui varie au cours du temps / phase de test """
    
    RL_profil_max = [0.64]*31+[0.64, 0.6, 0.5, 0.32, 0.16, 0.08, 0.04, 0.02, 0.01]
    n0 = ncouches_sol-t
    if t<ncouches_sol:
        res  = RL_profil_max[-1-t:-1]+[0]*n0
    else:
        res = RL_profil_max
    return res

#for i in range(50):print RLprof_t(i, ncouches_sol)
    



########## diverse plant fonctions - soil nitrogen balance



def critN (MS, a=4.8, b=-0.33):
    """ courbe critique de dilution de l'N - marche aussi pour array"""
    # MS = array od MS values (T.ha-1)
    vals = a*MS**b #en %
    if vals.size>1:
        for i in range(len(vals)): vals[i]=min(a, vals[i])
    else:
        vals = min(a, vals)
    return vals 


def demandeNdefaut(MSp,dMSp,Npc, surfsolref, a=4.8, b=-0.33):
    """ demande N pour parties aerienne - suppose meme courbe critique pour tout le monde - base sur N crit de la biomasse totale """
    #MSp = array des MSp (g.plant-1)
    #dMSp = array des dMSp (g.plant-1)
    #Npc = array des Npc plante (%)
    #surfsol sur laquelle sont les plantes #m2

    QN = MSp * Npc/100. #gN.plant-1
    MStot = array(sum(MSp+dMSp))/(surfsolref*100.)#MS new (T.ha-1)
    NcritTot = critN (MStot, a, b)#N crit de MS new
    PotN = (MSp + dMSp) * NcritTot/100. #gN.plant-1
    ls_demandeN = PotN-QN
    ls_demandeN[ls_demandeN<0.]=0.#gN.plant-1
    return ls_demandeN


#surfsolref = 0.05
#MSp = array([1.,1.2, 2.])
#dMSp = array([0.1,0.15,0.2])
#Npc = array([4., 3., 2.])
#demandeNdefaut(MSp,dMSp,Npc, surfsolref)

def demandeNdefaut2(MSp,dMSp,Npc, surfsolref, a=4.8, b1=-0.1 ,b2=-0.33):
    """ demande N pour parties aerienne - suppose meme courbe critique pour tout le monde - base sur N crit de la biomasse totale """
    #MSp = array des MSp (g.plant-1)
    #dMSp = array des dMSp (g.plant-1)
    #Npc = array des Npc plante (%)
    #surfsol sur laquelle sont les plantes #m2

    QN = MSp * Npc/100. #gN.plant-1
    MStot = array(sum(MSp+dMSp))/(surfsolref*100.)#MS new (T.ha-1)
    if MStot>=1.:
        NcritTot = a*MStot**b2#critN (MStot, a, b2)#N crit de MS new dense
    else:
        NcritTot = a*MStot**b1#critN (MStot, a, b1)#N crit de MS new isole

    #filtre NcritTot
    NcritTot[NcritTot>9.]=9.#gN.plant-1

    PotN = (MSp + dMSp) * NcritTot/100. #gN.plant-1
    ls_demandeN = PotN-QN
    ls_demandeN[ls_demandeN<0.]=0.#gN.plant-1
    return ls_demandeN, NcritTot, MStot #renvoie aussi Ncrit et MStot


def demandeNroot(MSpiv,dMSpiv,Npcpiv, surfsolref, Noptpiv):
    """ demande N pour parties racinaire - suppose N critique constant - s'applique aux racines et aux pivots """
    #MSp = array des MSp (g.plant-1)
    #dMSp = array des dMSp (g.plant-1)
    #Npc = array des Npc plante (%)
    #surfsol sur laquelle sont les plantes #m2

    QNpiv = MSpiv * Npcpiv/100. #gN.plant-1
    PotNpiv = (MSpiv + dMSpiv) * Noptpiv/100. #gN.plant-1
    ls_demandeN = PotNpiv-QNpiv
    ls_demandeN[ls_demandeN<0.]=0.#gN.plant-1
    return ls_demandeN




