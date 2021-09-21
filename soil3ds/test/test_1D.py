#############
# GL - test en 1D
# choix : ongletIn = 'Lusignan30_1' # pour sol avec 1 unique voxel/compartiment
# choix : ongletIn = 'Lusignan30' # pour sol avec 30 voxels verticaux
#############

from scipy import *
from soil3ds import soil_moduleN as solN
import soil3ds

import os

path_ = os.path.dirname(os.path.abspath(soil3ds.__file__))  # path ou trouver les inputs
path_leg = os.path.join(path_, 'test',
                        'inputs')  # r'C:\devel\l-egume\l-egume\input'#r'C:\devel\grassland'#r'H:\devel\grassland\grassland\L-gume' #r'C:\devel\grassland'

import IOxls


def init_sol(inis, meteo_j, par_sol, par_SN, Lsol, discret_solXY, dz_sol, pattern8, opt_residu, obstarac=None):
    """ soil initialisation comme dans L-py"""
    # vecteurs d'initialisation du sol
    Tsol = meteo_j['Tsol']  # 15. #degresC
    num_nb = list(map(int, inis['num_nb']))  # [6,6,18] #nbr de couche de chaque num de sol
    vsoilnumbers = [1] * num_nb[0] + [2] * num_nb[1] + [3] * num_nb[2]  # convention autorise 3 types d'horizon max
    # vDA = [par_SN['DA'][0]]*num_nb[0] + [par_SN['DA'][1]]*num_nb[1] + [par_SN['DA'][2]]*num_nb[2] #densite apparente de sol
    vCN = [par_SN['CN0_30']] * num_nb[0] + [par_SN['CN30_60']] * num_nb[1] + [par_SN['CN60_90']] * num_nb[
        2]  # maxi 3 horizons
    vMO = [par_SN['MO0_30']] * num_nb[0] + [par_SN['MO30_60']] * num_nb[1] + [par_SN['MO60_90']] * num_nb[
        2]  # maxi 3 horizons
    vARGIs = [par_SN['ARGIs0_30']] * num_nb[0] + [par_SN['ARGIs30_60']] * num_nb[1] + [par_SN['ARGIs60_90']] * num_nb[2]
    vCALCs = [par_SN['CALCs']] * ncouches_sol
    vNH4 = inis['NH4']  # [2.]*ncouches_sol # #!! kg d'N.ha-1 (entree de STICS)
    vNO3 = inis['NO3']  # [0.]*ncouches_sol
    HRpinit = inis['HRp']  # []
    if type(HRpinit) != type([0]):  # si 1 seule valure/Scalaire
        vNH4 = [inis['NH4']]
        vNO3 = [inis['NO3']]
        HRpinit = [inis['HRp']]

    if min(HRpinit) < 0:  # code -1 pour pas d'initialisation
        HRpinit = []

    vDA = []
    for i in vsoilnumbers:
        vDA.append(par_sol[str(i)]['DA'])

    # vsoilnumbers = [1]+[2]*3+[3]*13+[4]*13 #numeros de sol du profil -> mesures acsyd11
    # vDA = [1.81]+[1.31]*3+[1.37]*13+[1.42]*13 #densite apparente de sol (mesure pesees initial aschyd11)
    # vCN = [par_SN['CN0_30']]*ncouches_sol #maxi 90cm en strates de 5cm
    # vMO = [par_SN['MO0_30']]*ncouches_sol #maxi 90cm en strates de 5cm
    # vARGIs = [par_SN['ARGIs']]*ncouches_sol #maxi 90cm
    # vCALCs = [par_SN['CALCs']]*ncouches_sol
    # vNH4 = [2.]*ncouches_sol # #!! kg d'N.ha-1 (entree de STICS)
    # coeff = 0.#0.09#coeff perte ressuyage -> a ajuster pour avoir environ 600 kg N.ha-1
    # vNO3 = [91.*coeff]*ncouches_sol # kg d'N.ha-1 (entree de STICS)
    # vNO3 = array([16.96, 16.07, 15.17, 33.92, 33.92, 33.92, 33.92, 62.49, 82.13, 89.27, 76.77, 107.13, 124.98, 142.84, 124.98, 142.84, 160.69, 151.76, 151.76, 142.84, 178.55, 133.91, 98.20, 89.27, 83.92, 89.27, 73.20, 89.27, 87.45, 62.49])*coeff #issu du profil en sol nu
    # HRpinit = [25.5,26.,25.,25.5,26.,26.,26.,26.5,26.5,27.,27.,27.,27.5,27.5,27.5,27.5,27.5,29,29,29,29,29,29,29,29,30,30,30,30,30]#-> mesures ahscyd au jour 195 (140711) -> init sol nu

    ## soil initialisation
    S = solN.SoilN(par_sol, par_SN, soil_number=vsoilnumbers,
                   dxyz=[[Lsol / discret_solXY[0]] * discret_solXY[0], [largsol / discret_solXY[1]] * discret_solXY[1],
                         [dz_sol / 100.] * ncouches_sol], vDA=vDA, vCN=vCN, vMO=vMO, vARGIs=vARGIs, vNO3=vNO3,
                   vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol, pH=par_SN['pH'], ZESX=par_SN['ZESX'], CFES=par_SN['CFES'],
                   obstarac=obstarac, pattern8=pattern8)

    if HRpinit != []:  # initialise humidite si un vecteur est fourni
        S.init_asw(HRp_init=HRpinit)

        # lims_sol = rtd.lims_soil(pattern8, dxyz=[[Lsol], [largsol], [dz_sol/100.]*ncouches_sol])

    if opt_residu == 1:  # initialisatio de residus
        S.init_residues(vCNRESt, vAmount, vProps, vWC, vCC)

    # print 'sol', sum(S.m_NO3), sum(S.m_NH4), sum(S.m_QH20fc)-sum(S.m_QH20wp)

    # Uval = 0.9*2.61#(epaisseur de sol* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
    Uval = par_SN['q0'] * 0.1 * sum(S.m_QH20fc[0]) * surfsolref / (S.dxyz[2][
                                                                       0] * 100.)  # (epaisseur de sol (cm)* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
    stateEV = [0., 0., 0.]  # pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
    b_ = solN.bEV(par_SN['ACLIMc'], par_SN['ARGIs'], HXs=0.261)  # 1.#valeur empirique tres proche#0.1#0.63#0.63
    # !!!

    return S, Tsol, Uval, stateEV, b_


def critN(MS, a=4.8, b=-0.33):
    """ courbe critique de dilution de l'N """
    return min(6.5, a * MS ** b)  # en %


## 1) lecture fichier initialisation
meteo_path = os.path.join(path_leg,
                          'meteo_exemple.xls')  # 'meteo_exemple_debugL_gl.xls')#r'H:\devel\grassland\grassland\L-gume\meteo_exemple2.xls'
ongletM = 'Lusignan30'  # 'Lusignan302ans'#'DivLeg15'#'morpholeg15'#'combileg15'#'combileg16'#'Avignon30'#'exemple'#'morpholeg15'#'testJLD'#'competiluz'#
meteo = IOxls.read_met_file(meteo_path, ongletM)

## lecture fichier management
mn_path = os.path.join(path_leg,
                       'management_exemple.xls')  # 'management_exemple3_debugL_gl.xls')#r'H:\devel\grassland\grassland\L-gume\management_exemple.xls'
ongletMn = 'Lusignan30IrrN2'  # 'Lusignan30IrrN2ans'#'DivLeg15'#'Lusignan30IrrN'#'illimite-sanscoupe'#'combileg15-irrigajusteeLUZTVMIN'#'combileg16-irrigajusteeMIN'#'Lusignan30'#'Avignon30IrrN'#'Avignon30'#
mng = IOxls.read_met_file(mn_path, ongletMn)

inis_path = os.path.join(path_leg, 'Init_sol_exemple.xls')  # 'Initialisation_sol_exemple.xls')
ongletIn = 'Lusignan30'#'Lusignan30_1'  #
inis = IOxls.read_plant_param(inis_path, ongletIn)

# lecture des parametres du sol
path_sol = os.path.join(path_leg, 'Parametres_sol_exemple.xls')  # 'Parametres_sol_exemple2_debugL_glbis.xls')#
ongletS = 'lusignan99'  # 'morpholeg'#'combileg2015vshallow'#'combileg16vshallow'#'ASCHYD11'#
par_SN, par_sol = IOxls.read_sol_param(path_sol, ongletS)

# Param Plante
plant_path = os.path.join(path_leg, 'Parametres_plante_exemple.xls')  # 'Initialisation_sol_exemple.xls')
ongletP = 'Orca'  #
ParamP = IOxls.read_plant_param(plant_path, ongletP)

# 2) definition du pattern et discretisation sol
cote = 100.
pattern8 = [[0, 0], [cote, cote]]
Lsol = max((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
largsol = min((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
surfsolref = Lsol * largsol  # m2
dz_sol = inis['dz_sol']  # 4.#5. #cm
ncouches_sol = int(inis['ncouches_sol'])  # 4#10#30
prof_sol_max = ncouches_sol * dz_sol  # 80.

discret_solXY = list(map(int, inis['discret_solXY']))  # [10,10]# nb de discretisation du sol en X et en Y
# lims_sol = rtd.lims_soil(pattern8, dxyz=[[Lsol / discret_solXY[0]] * discret_solXY[0],
#                                         [largsol / discret_solXY[1]] * discret_solXY[1],
#                                         [dz_sol / 100.] * ncouches_sol])

opt_residu = 0

# debut, fin de simulation
DOY_deb, DOY_fin = 100, 300  # 239,623

# initialisation sol
meteo_j = IOxls.extract_dataframe(meteo, ['TmoyDay', 'I0', 'Et0', 'Precip', 'Tsol'], 'DOY', val=DOY_deb)
S, Tsol, Uval, stateEV, b_ = init_sol(inis, meteo_j, par_sol, par_SN, Lsol, discret_solXY, dz_sol, pattern8, opt_residu,
                                      obstarac=None)

# simulation d'un sol 1D
##vegetation avec racine, LAI et MSaerien constants
R1 = S.m_1 * 200.  # vert_roots(S.dxyz, [0.000000001,0.,0.,0.]) #pas zero sinon buf FTSW
# R1[0,:,:] = R1[0,:,:]+0.000000001 #ajoute epsilon ds 1er horizon
ls_roots = [R1]
ls_epsi = [0.2]

# initialise teneur en N des plantes
Npc = 2.  # %
MSa = 1.5  # T.ha-1
QN = MSa * Npc / 100. * 1000  # kg N.ha-1 #%N libre

# initialisation de variables de sorties
cumEV, cumET0, cumPP, cumD, profH20, cumTransp, vlix, azomes = [], [], [], [], [], [], [], []

##boucle journaliere couplage sol-plante
for DOY in range(DOY_deb, DOY_fin):

    # MAJ meteo / mng
    # meteo_j = IOxls.extract_dataframe(meteo, ['TmoyDay','I0','Et0','Precip','Irrig','Coupe','FertNO3','FertNH4','Tsol'], 'DOY', val=DOY)
    meteo_j = IOxls.extract_dataframe(meteo, ['TmoyDay', 'I0', 'Et0', 'Precip', 'Tsol'], 'DOY', val=DOY)
    mng_j = IOxls.extract_dataframe(mng, ['Coupe', 'Irrig', 'FertNO3', 'FertNH4', 'Hcut'], 'DOY', val=DOY)
    print(DOY)
    for k in list(meteo_j.keys()): meteo_j[k] = meteo_j[k][0]
    for k in list(mng_j.keys()): mng_j[k] = mng_j[k][0]

    # entrees eau
    # Precip = meteo_j['Precip']+meteo_j['Irrig']
    Rain = meteo_j['Precip']
    Irrig = mng_j['Irrig']

    # entrees N
    # map_N = 0.*S.m_1[0,:,:]
    mapN_Rain = 1. * S.m_1[0, :, :] * Rain * par_SN['concrr']  # Nmin de la pluie
    mapN_Irrig = 1. * S.m_1[0, :, :] * Irrig * par_SN['concrr']  # Nmin de l'eau d'irrigation
    mapN_fertNO3 = 1. * S.m_1[0, :, :] * mng_j['FertNO3'] * S.m_vox_surf[0, :, :] / 10000.  # kg N par voxel
    mapN_fertNH4 = 1. * S.m_1[0, :, :] * mng_j['FertNH4'] * S.m_vox_surf[0, :, :] / 10000.  # kg N par voxel

    # entre Tsol (lecture)
    S.updateTsol(meteo_j['Tsol'])  # (meteo_j['TmoyDay']) #Tsol forcee comme dans STICS

    # demande N plante pour 1 couvert

    PotN = MSa * critN(MSa) / 100. * 1000  # kg N.ha-1
    demande_N_plt = max(PotN - QN, 0.)  # kg N.ha-1
    ls_demandeN = [sum(demande_N_plt) / 10000.]  # kg N.surface de sol

    # Calcul du bilan hydrique
    ls_transp, evapo_tot, Drainage, stateEV, m_frac_transpi, m_frac_evap, ls_ftsw = S.stepWBmc(
        meteo_j['Et0'] * surfsolref, ls_roots, ls_epsi, Rain * surfsolref, Irrig * surfsolref, stateEV,
        ZESX=par_SN['ZESX'], leafAlbedo=0.15, U=Uval, b=b_, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)

    # Calcul du bilan N

    # Mineralisation / Nitif /Inflitration N
    S.stepNB(par_SN)
    S.stepNitrif(par_SN)
    S.stepNINFILT(mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, Drainage, opt=1)

    # uptake plante
    ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N, idmin = S.stepNuptakePlt(par_SN, [ParamP], ls_roots, m_frac_transpi,
                                                                     ls_demandeN)

    # update QN et Npc
    QN += sum(ActUpNtot) * 10000  # kg N.ha-1
    Npc = QN / (MSa * 10)

    # sorties
    print(DOY, 'tsw_t: ', S.tsw_t[0, 0, 0], 'evapotot: ', evapo_tot)  # sum3(S.tsw_t)
    print('Npc', Npc)
    # print('QN', QN)
    HAx = S.HRp()
    cumEV.append(evapo_tot)
    cumTransp.append(sum(ls_transp))
    cumET0.append(meteo_j['Et0'] * surfsolref)
    cumPP.append(meteo_j['Precip'] * surfsolref)
    cumD.append(Drainage[-1][0][0])
    profH20.append([DOY] + HAx[:, 0, 0].tolist() + [evapo_tot, Drainage[-1][0][0]])

    vlix.append(S.lixiNO3 * 10000)
    azomes.append(S.m_NH4.sum() + S.m_NO3.sum() * 10000)

##termes du bilan hydrique global
S.CloseWbalance()  # -> equilibre
S.CloseCbalance()  # -> equilibre
S.CloseNbalance()  # -> equilibre










