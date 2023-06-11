#############
# GL - test default soil functions
# choix : ongletIn = 'Lusignan30_1' # pour sol avec 1 unique voxel/compartiment
# choix : ongletIn = 'Lusignan30' # pour sol avec 30 voxels verticaux
#############


import os
import soil3ds
from soil3ds import soil_moduleW as solW
from soil3ds import soil_moduleN as solN
from scipy import *

from legume import initialisation # require legume package for 'init_sol_fromLpy' function
from legume import IOxls
#from soil3ds import IOxls

path_ = os.path.dirname(os.path.abspath(soil3ds.__file__))  # path ou trouver les inputs
path_leg = os.path.join(path_, 'test', 'inputs')


## 1) lecture fichier initialisation
meteo_path = os.path.join(path_leg, 'meteo_exemple.xls')  # 'meteo_exemple_debugL_gl.xls')#r'H:\devel\grassland\grassland\L-gume\meteo_exemple2.xls'
ongletM = 'Lusignan30'  # 'Lusignan302ans'#'DivLeg15'#'morpholeg15'#'combileg15'#'combileg16'#'Avignon30'#'exemple'#'morpholeg15'#'testJLD'#'competiluz'#
meteo = IOxls.read_met_file(meteo_path, ongletM)

## lecture fichier management
mn_path = os.path.join(path_leg, 'management_exemple.xls')  # 'management_exemple3_debugL_gl.xls')#r'H:\devel\grassland\grassland\L-gume\management_exemple.xls'
ongletMn = 'Lusignan30IrrN2'  # 'Lusignan30IrrN2ans'#'DivLeg15'#'Lusignan30IrrN'#'illimite-sanscoupe'#'combileg15-irrigajusteeLUZTVMIN'#'combileg16-irrigajusteeMIN'#'Lusignan30'#'Avignon30IrrN'#'Avignon30'#
mng = IOxls.read_met_file(mn_path, ongletMn)

# inis_path = os.path.join(path_leg, 'Init_sol_exemple.xls')  # 'Initialisation_sol_exemple.xls')
# ongletIn = 'Lusignan30'#'Lusignan30_1'  #
# inis = IOxls.read_plant_param(inis_path, ongletIn)

# lecture des parametres du sol par defaut
par_sol = solW.default_par_sol()
par_SN = solN.default_parSN()

# Param Plante defaut

ParamP =  {'Vmax1': 0.0018,
 'Kmax1': 50.0,
 'Vmax2': 0.05,
 'Kmax2': 25000.0}

# faire fonction par defaut




# 2) definition du pattern et discretisation sol
cote = 100.
pattern8 = [[0, 0], [cote, cote]]
Lsol = max((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
largsol = min((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
surfsolref = Lsol * largsol  # m2
dz_sol = 5. #inis['dz_sol']  #cm
ncouches_sol = 30 #int(inis['ncouches_sol'])
prof_sol_max = ncouches_sol * dz_sol
discret_solXY = [1, 1] #list(map(int, inis['discret_solXY']))


# vecteurs profils qui devraient venir de inis
vsoilnumbers = [1]*ncouches_sol
vDA = [1.42]*ncouches_sol
vCN = [9.52]*ncouches_sol
vMO = [18.02]*ncouches_sol
vARGIs = [18.3]*ncouches_sol
vCALCs = [0.2]*ncouches_sol
vNH4 = [0.2]*ncouches_sol #inis['NH4'] #
vNO3 = [1.]*ncouches_sol #inis['NO3'] #


opt_residu = 0

# debut, fin de simulation
DOY_deb, DOY_fin = 100, 300  # 239,623

# initialisation sol
meteo_j = IOxls.extract_dataframe(meteo, ['TmoyDay', 'I0', 'Et0', 'Precip', 'Tsol'], 'DOY', val=DOY_deb)

S = solN.SoilN(par_sol, par_SN, soil_number=vsoilnumbers,
                   dxyz=[[Lsol / discret_solXY[0]] * discret_solXY[0], [largsol / discret_solXY[1]] * discret_solXY[1],
                         [dz_sol / 100.] * ncouches_sol], vDA=vDA, vCN=vCN, vMO=vMO, vARGIs=vARGIs, vNO3=vNO3,
                   vNH4=vNH4, vCALCs=vCALCs, Tsol=meteo_j['Tsol'], obstarac=None, pattern8=pattern8)



# simulation d'un sol 1D
## pour sol nu
epsilon = 0.5#10-10
R1 = S.m_1 * epsilon
ls_roots = [R1]
ls_epsi = [0.2]
ls_N = [0.9]
ls_paramP = [ParamP]


# initialisation de variables de sorties
cumET0, cumPP = [], []

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

    res_soil_step = solN.step_bilanWN_solVGL(S, par_SN, meteo_j, mng_j, ls_paramP, ls_epsi, ls_roots, ls_N, opt_residu, opt_Nuptake=1)
    S, stateEV, ls_ftsw, ls_transp, ls_Act_Nuptake_plt, temps_sol = res_soil_step

    # sorties
    cumET0.append(meteo_j['Et0'] * surfsolref)
    cumPP.append(meteo_j['Precip'] * surfsolref)

    # S.bilanW['cumEV'] # soil evaporation
    # S.bilanW['cumTransp'] # transpiration
    # S.bilanW['cumD'] # drainage
    # S.bilanN['cumMinN'] # soil mineralisation
    # S.bilanN['cumLix'] #lixiviation nitrates
    # S.bilanN['azomes'] # total minearal N


##termes du bilan hydrique global
S.CloseWbalance()  # -> equilibre
S.CloseCbalance()  # -> equilibre
S.CloseNbalance()  # -> equilibre


#?? tourne mis pas d'uptake plant??






