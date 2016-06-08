from PySteppables import *
from PySteppablesExamples import MitosisSteppableBase
import CompuCell
import CompuCellSetup
import sys
from random import *
import math
import time
# from math import *
import numpy
from datetime import datetime
import cProfile
# from scipy import integrate

""" "TO GET ALL CELL ATTRIBUTES" (to see what cell attributes can be accessed/changed in Python):
print dir(cell)
output 7-2-2012:
['__class__', '__del__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattr__', '__getattribute__', '__hash__', '__init__', '__module__',
'__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__swig_destroy__', '__swig_getmetho
ds__', '__swig_setmethods__', '__weakref__', 'angle', 'averageConcentration', 'clusterId', 'clusterSurface', 'ecc', 'extraAttribPtr', 'flag', 'fluctAm
pl', 'iXX', 'iXY', 'iXZ', 'iYY', 'iYZ', 'iZZ', 'id', 'lX', 'lY', 'lZ', 'lambdaClusterSurface', 'lambdaSurface', 'lambdaVecX', 'lambdaVecY', 'lambdaVec
Z', 'lambdaVolume', 'pyAttrib', 'subtype', 'surface', 'targetClusterSurface', 'targetSurface', 'targetVolume', 'this', 'type', 'volume', 'xCM', 'xCOM'
, 'xCOMPrev', 'yCM', 'yCOM', 'yCOMPrev', 'zCM', 'zCOM', 'zCOMPrev']
"""


#STRINGS FOR OUTPUT FILES
drug='cisplatin'
#size='lgSphere'
#vess='451v_0.99pctA'
#vess='151v_2.17pctA'
# size='smSphereIV'
# size='smSphereIP'
# size='smSphereIPIV'
# size='smSphereIV'
# size='smSphereIVinfus'
# size='smSphereIPinfus'
# size='R140CubeEdgeIVinfus'
# size='R140CubeEdgeIPinfus'
size='R140CubeCenterIVinfus'
# size='R140CubeCenterIPinfus'
# size='smSphereIPIV'
# vess='18v_2.26pctA'
# vess='17v_2.26pctA'
# vess='16v_7.64pctA'
# vess='11v_10.89pctA'
# vess='11v_6.65pctA'
# vess='11v_4.53pctA'
# vess='10v_2.26pctA'
# vess='8v_4.38pctA'
# vess='6v_3.82pctA'
# vess='3v_10.04pctA'
# vess='2v_1.27pctA'

# vess='1v_3pctvol_2pctA'
# vess='7v_5pctvol_4pctA'
# vess='6v_8pctvol_6pctA'
# vess='15v_10pctvol_8pctA'
# vess='19v_12pctvol_10pctA'
# vess='0v_0.0pctA'
vess='_10pctA' # for cubes



#INITIALIZE COUNTING VARIABLES
numPCancer=0
numIC50Cancer=0


#INITIALIZE CHEMICAL THRESHOLDS

# FINAL CISPLATIN IC50 = MISTRY 1992
cisplatinIC50=52.71 # FINAL USED: muM (equiv to (equitoxic) 2h-IC-50 Mistry, 1992)
# cisplatinIC50=38.3 # (MATCHES 2H DRUG TIME COURSE): muM (SD=12.6, in SKOV-3, 2h exposure; Table 1, Mistry, 1992)
# cisplatinIC50=0.0001 #test
# cisplatinIC50=126.63 # muM (+/- 12.06 micromols/L, in SKOV3ip1, 48h exposure; Fig. 3, Xu, 2008)
# cisplatinIC50  =3.33 # muM (Nakaro 1997, minimum effective concentration for DNA damage)

#IV CISPLATIN
# first five minutes (no data, linear fit from t=0, [iv]=0 to first time point
# cFirst5Mins=5*465.189 # five minutes worth of MCSes in normal tissue (Swabb 1974?)
cFirst5Mins=5*1207.183 # five minutes worth of MCSes in tumor tissue in vivo (window chamber) for sodium fluorescein (~376Da,like Cisplatin,300Da (Nugent, 1984))
cIVFirstPoint=(5.742+15)*1207.183 # 6931.79 MCS; five minutes worth of MCSes in tumor tissue in vivo (window chamber) for sodium fluorescein (~376Da,like Cisplatin,300Da (Nugent, 1984))
cFirst20Mins=20*1207.183 # five minutes worth of MCSes in tumor tissue in vivo (window chamber) for sodium fluorescein (~376Da,like Cisplatin,300Da (Nugent, 1984))
Infusion15Mins=15*1207.183 # 18107.745 MCS; subtract from runtime once begin using concentration time course
# x-intercept of cisplatin IV concentration function
# cIVZero=57751. #mcs hrs, MCS scaled for diffusion in normal tissue
# cIVZero=math.exp(-16.094/-3.338)*1207.183 ##mcs mins, MCS scaled for diffusion in tumor tissue;Excel check, yintercept = 124.1449668
# x-intercept of cisplatin IP-post-IV concentration function
# cIPpostIV_zero=(0.386/0.5431)*465.... #mcs hrs, MCS scaled for diffusion in normal tissue
# cIPpostIVZero= (0.386/0.0091)*1207.183 #mcs mins, MCS scaled for diffusion in tumor tissue; IP will become negative after this time point (42.42mins = 51205.8mcs)
#  cIPpostIVZero= 896.305 #mcs mins, MCS scaled for diffusion in tumor tissue; IP will become negative after t = 896.305 (Wolfram Alpha cubic polynomial solver).  Will not use because becomes zero after IV goes to zero; then soln. is zero
#  cIVpostPZero=~ 7hrs #mcs mins, MCS scaled for diffusion in tumor tissue; IP will become negative after t = 7 hrs

#ip CISPLATIN
# cIPFirstPt=0.2126*60*1207.183 #hrs, first time point
# cIPFirstPt=7.386*1207.183 #7.386 min, first time point
cIPFirstPtPlus15=(7.386+15.0)*1207.183 #= 27023.998638 MCS = 7.386 min, first time point, plus 15 mins infusion time; initialize with floats
#cIVFirstPtPlus15=(5.742+15.0)*1207.183 # = 25039.389786 MCS = 5.742 min, first time point, plus 15 mins infusion time


Cisp1Min = 65.678 # 65.678 mcs = 1 min for cells of diameter T24 bladder cancer cell line

## cisplatin, platinum accumulation per cell per time step, based on IC50 of bladder cancer cell line;
## also = concentration removed from voxel; µM/MCS * siteConcCis(µM)
cispAccumFrac_SCSG_BFTC_905 = 0.3147446166  # (sens cis and gem)	2.575477619	IC50 µM	cisplatin			
cispAccumFrac_SCSG_J82 = 0.3033715255	    # (sens cis and gem)	5.42972235	IC50 µM	cisplatin				
cispAccumFrac_RCRG_RT4 = 0.2154448963	    # (resist cis and gem)	27.49620513	IC50 µM	cisplatin				
cispAccumFrac_RCRG_HT_1197 = 0.1498013738   # (resist cis and gem)	43.97041406	IC50 µM	cisplatin				
cispAccumFrac_SCRG_SW780 = 0.2691142767     # (sens cis resist gem)	14.02708355	IC50 µM	cisplatin					
cispAccumFrac_SCRG_KU_19_19 = 0.2847440491  # (sens cis resist gem)	10.10456195	IC50 µM	cisplatin
cispAccumFrac_RCSG_LB831_BLC = 0.02925373453# (resist cis sens gem)	225.1619678	IC50 µM	cisplatin
cispAccumFrac_RCSG_DSH1 = 0.04556319402     # (resist cis sens gem)	144.5646771	IC50 µM	cisplatin				

## gemcitabine (possibly dFdCtP) accumulation per cell per time step, based on IC50 of bladder cancer cell line;
## also = concentration removed from voxel; µM/MCS * siteConcGem(µM)
gemAccumFrac_SCSG_BFTC_905 = 4.41575E-03  # (sensitive cis and gem)	5.15E-06	IC50 µM	gemcitabine			
gemAccumFrac_SCSG_J82 = 4.41565E-03       # (sensitive cis and gem)	0.007409799	IC50 µM	gemcitabine				
gemAccumFrac_RCRG_RT4 = 4.22858E-03       # (resistant cis and gem)	13.84278281	IC50 µM	gemcitabine				
gemAccumFrac_RCRG_HT_1197 = 4.39190E-03   # (resistant cis and gem)	1.764248816	IC50 µM	gemcitabine				
gemAccumFrac_SCRG_SW780 = 2.68443E-03     # (sens cis resist gem)	128.0487435	IC50 µM	gemcitabine				
gemAccumFrac_SCRG_KU_19_19 = 4.18843E-03  # (sens cis resist gem)	16.81235445	IC50 µM	gemcitabine				
gemAccumFrac_RCSG_LB831_BLC = 4.41518E-03 # (resist cis sens gem)	0.041854289	IC50 µM	gemcitabine				
gemAccumFrac_RCSG_DSH1 = 4.41445E-03      # (resist cis sens gem)	0.096498675	IC50 µM	gemcitabine				






# PRINT SIMULATION START TIME
now = datetime.now()
print "SIMULATION START TIME =",now




class CispIC50VisualizationSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.scalarCLField=self.createScalarFieldCellLevelPy("NetAccumulatedCispl")
    def step(self,mcs):
        self.scalarCLField.clear()
        for cell in self.cellList:
            if cell.type==2 or cell.type==3 or cell.type==8:
                dictionaryAttrib = CompuCell.getPyAttrib(cell)
                self.scalarCLField[cell]=dictionaryAttrib[3]




# SETCELLCLOCKS
# SETS UNIQUE INITIAL AGE FROM UNIFORM DISTRIBUTION
# STARTS INTERNAL CLOCK IN CELL
# SETS DIVISION RATE FOR CELL, SELECTED FROM GAUSSIAN DISTRIBUTION
# SETS INITIAL CISPLATIN CONCENTRATION
# SETS INITIAL PERTUZUMAB CONCENTRATION
class SetCellDictionaries(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

    def start(self):
        print "This function (SetCellDictionaries) is called once before simulation"
        self.cellList=CellList(self.inventory)
        for cell in self.cellList:

            # SAME CLOCK IN ALL CELLS
            # ALSO CHANGE IN MITOSISSTEPPABLE CLASS (~LINE 361)
            x = gauss(24,1)
            y = uniform(0,24) # age of cells initialized into simulation
            dictionaryAttrib[0:2]=[x,y]
            print dictionaryAttrib[0],dictionarAttrib[1],dictionaryAttrib[2]

            # # # DIFFERENT CLOCK FOR SELECTED CELL TYPES
            # # PCancerGFP,RFP clocks
            # if cell.type==2:
            #     #mitosis clock --- 24+/-6 hours (Abbas) -- dictionaryAttrib[0:2]=[age_to_divide,current_age]
            #     dictionaryAttrib = CompuCell.getPyAttrib(cell)
            #     # x = gauss(25.5,1)
            #     # y = uniform(0,25.5) # age of cells initialized into simulation
            #     # dictionaryAttrib[0:2]=[x,y]
            #     # INITIAL CISPLATIN ACCUMULATED IN CELL
            #     dictionaryAttrib[3]=0.0
            #     # INITIAL ANTIBODY ACCUMULATED ON CELL (ADSORBED BY ErbB2)
            #     # dictionaryAttrib[4]=0.0
            # elif cell.type==3:
            #     # INITIAL CISPLATIN ACCUMULATED IN CELL
            #     dictionaryAttrib[3]=0.0
            #     # INITIAL ANTIBODY ACCUMULATED ON CELL (ADSORBED BY ErbB2)
            #     # dictionaryAttrib[4]=0.0



# SECRETIONSTEPPABLE
# SET CELL SECRETION
# OCCURS PRIOR TO DIFFUSION IN MCS (XML OCCURS AFTER)
# CISPLATIN SOURCES AND SINKS IN EACH TIME STEP:
   #"SECRETION" BY CELL TYPES VESSELWALL, PCANCERGFP, PCANCERRFP, MEDIUM;
   # ACCUMULATION BY TUMOR CELLS;
   # EFFLUX OF ACCUMULATED CISPLATIN? ADD IF NEEDED.
### SOURCES:
  # <!-- CISPLATIN -->
  # <!-- (DIFFUSION COEFFICIENT OF CREATININE (USED AS PROXY FOR CISPLATIN in Morrison, 1986) = 1.9 x 10^- 6 cm^2/s = 0.000114 cm^2/min = 11400 um^2/min = ) -->
  # <!-- D_{NORMAL TISSUE}  (Swabb 1974 via Thurber 2011)
  #      = 1.778 x 10^{-4} (300.5)^{-0.75} = 0.000002466246792 cm^2/s = 246.624679177701 um^2/s
  #      = 4652 (SKOV3.ip1 cell diameter)^2 / (10 min)
  #      = 1 celldiam^2/(1/465.18914889411 min) = 227911.3489336466 MCS/hour = 55822.6978672932 MCS/2 hours
  # D(VX2 carcinoma for sodium fluorescein, MW376 (Nugent 1984)
  #      = 1207.18273728686 cell diam^2 / 1 min (= 5.64um^2 (voxel edge) / 1/60hr)
  #      = 1 cell diam^2 / 1/1207.183 min (= 5.64um^2 (voxel edge) / 1/1207.183min)
  #      = 72430.9642372114 MCS/hour = 144861.928474423 2 hours-->
  # <!-- CISPLATIN (PT) ACCUMULATION AT 5 $\MU$M IN SKOV3  (Fig.3, Mistry, 1992)
  #      = 9.4E-15 $\mu$M/cell/min = 1.57E-16 muM/cell/sec = 9.4E-14 $\mu$M/cell/10 min
  #      accumulation per SKOV3 cell, with respect to current concentration (uM)
  #      = (- 1.2e-07*x^{3} + 1.3e-05*x^{2} + 0.00058*x + 0.076)/10/465.189 uM, x = current concentration -->
  # <!-- ACCUMULATION at 100 micromolar for 50min in SKOV3 
  #      = 8 +/- 0.8 nmol/mg protein/50 mins: 8 +/- 0.8 nM * 1e-3uM/nM /mg protein * 0.21mg protein/uL * 1uL/10^9um^3 * 179.1 um^3/ SKOV3 cell /50min = 6.02E-12 uM Pt/SKOV3 cell/min = 1.00E-13 uM Pt/SKOV3 cell/s -->
  # <!-- PEAK AVE INTRAPERITONEAL
  #      (90mg/m^2 IP for 4 hour dwell = peak ave IP 40.5mug/mL; peak ave IV 1.6mug/mL; Howell, 1982): 
  #      40.5mug/mL*1g/10^6mug*1000mL/L/300.5g/mol cisplatin*10^6umol/mol = 134.977503749375 muM = 8098.65022496251 (concentration * 60 sec/MCS)-->
  # <!-- PEAK VALUES IV
  #      (90mg/m^2 IP = peak ave IP 40.5mug/mL; peak ave IV 1.6mug/mL; 1.6mug/mL*1g/10^6mug*1000mL/L/300.5g/mol cisplatin*10^6umol/mol =  5.3324445925679 uM,  
  #      = 319.946675554074 (concentration * 60 sec/MCS) -->
  # <!-- MINIMUM EFFECTIVE CONCENTRATION OF CISPLATIN (BASED ON DNA DAMAGING LEVELS in clear cell ovarian carcinoma, Takatori, 2012)
  #      = 1 mug/mL; 1ug/ml*1000ml/L*1g/10^6ug* 1/300.5g/mol*10^6umol/mol = 3.33 muM-->
  # <!-- MINIMUM EFFECTIVE CONCENTRATION OF CISPLATIN (IN SERUM, clinical trials, lung cancer, Nakano, 1997)
  #      = 0.3micrograms/ml; 0.3ug/ml*1000ml/L*1g/10^6ug* 1/300.05g/mol*10^6umol/mol = 0.999833361106482 muM -->
  # <!-- IC50 for SKOV3.ip1 cells (Xu, W et al. "Antisense Oligodeoxynucleotides... " IntJGynCanc)  -->
  # <!-- SECRETION = micromolar = umol/
  #      (90mg/m^2IP = peak ave IP 40.5mug/mL; peak ave IV 1.6mug/mL 1.6 1.6mug/mL*1g/10^6mug*1000mL/L/Xg/mol cisplatin*10^6umol/mol 5.3324445925679 uM 319.946675554074 (concentration * 60 sec/MCS) -->
  #      EFFLUX  (Fig.4, Mistry, 1992 ) = 3.58E-04 % accumulated Pt effluxed/ 1/465.189min = 3.58E-06 fraction of accumulated Pt effluxed/ 1/465.189min for D in normal tissue OR 1.38E-06 fraction of accumulated Pt effluxed/ 1/1207.183min for D in tumor tissue
  # <!-- DECAY RATE (not used when empirically determined effective diffusion coefficient is used) = ln(2)/half life = 0.23 / 10 min = 0.000049667766485 /(1/465.189min) (Go, Ata, RXMed:Platinol) -->
# Cisplatin intraperitoneal concentration fitted curve (dosage = 270 mg/m^2 ) 451.51e-0.551*t muM, t=   Howell, 1982 (Fig. 2)
# Cisplatin post-IP plasma concentration fitted curve [IP(t)]*-0.000889*t+0.06288 muM, t= Sugarbaker, 1986, 1990, 1997
# Cisplatin intravenous concentration fitted curve (dosage =100mg/m^2) -3.338ln(t) + 16.094 muM, t=min   Himmelstein, 1981 (Fig. )
# Cisplatin post-IV peritoneal concentration fitted curve [IV(tMins)]*( -2E-07*tMins^3 + 0.0002*tMins^2 - 0.0197*tMins + 0.9963) # goes to zero when t = 896.305 (Wolfram Alpha cubic polynomial solver); Linear solution:  (0.5431*t - 0.386) muM, t=h; y = [IV(t)]*0.0091t - 0.386, R^2= 0.9509, t=,min; Sugarbaker, 1996

class SecretionSteppableCisplatin(SecretionBasePy,SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SecretionBasePy.__init__(self,_simulator, _frequency)
    def start(self):
        pass
    def step(self,mcs):
        ## START PROFILER
        # profile = cProfile.Profile()
        # profile.enable()

        for cell in self.cellList:
            #if cell.type==2:  # cancer cell that has reached the therapeutic threshold
                if cell.type==2 or cell.type==8:  #any cancer cell -- limit of cisplatin absorption not known
                    comPt=CompuCell.Point3D()
                    field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                    # WORKS WHEN cell vol = 1 voxel; changed for tiny speed-up; otherwise, use comPt.x=int(round(cell.xCM/float(cell.volume))) (int(averageCOM))(see CC3D manual)
                    comPt.x=int(cell.xCM)
                    comPt.y=int(cell.yCM)
                    comPt.z=int(cell.zCM)
                    cisplatin=field.get(comPt) # get concentration at center of mass
                    attrSecretor=self.getFieldSecretor("Cisplatin")

                    #ACCUMULATE (ADD TO DICTIONARY, REMOVE FROM FIELD)
                    if cisplatin > 0:
                        dictionaryAttrib = CompuCell.getPyAttrib(cell)
                        # ADD EMPIRICALLY-DETERMINED FRACTION OF CURRENT CONCENTRATION AT CELL COM TO ACCUMULATED CONCENTRATION IN CELL
                        # REMOVE ACCUMULATED DRUG FROM EXTERNAL CONCENTRATION
                        # ******fraction of what would have accumulated in 2 hrs at current concentration at COM (Mistry, 1992)
                        MCSFrac2Hrs=1/Cisp1Min/60/2 # MCSFrac2Hrs= 1 mcs / number of mcs in 2 h  =  1 mcs / (Cisp1Minmcs/min * 60min/hr * 2hr)
                        accumC=(0.4755*cisplatin**1.289)*MCSFrac2Hrs #  Accumulation (micromolar) = 0.4755*IncubationConc(micromolar)^1.289, R-square: 0.9999
                        dictionaryAttrib[3]+=accumC
                        attrSecretor.uptakeInsideCellAtCOM(cell,accumC,1.0) # uM secretion from pixels at outer boundary of cell

                 #     # REMOVE EMPIRICALLY-DETERMINED EFFLUXED DRUG FROM ACCUMULATED DRUG IN CELL
                #     # IF IMPLEMENTING EFFLUX AFTER EXTERNAL DRUG = 0, ADD EFFLUX FROM TOTAL DRUG ACCUMULATED PER CELL BACK INTO EXTERNAL DRUG CONCENTRATION. Net accumulation empirical measurements include efflux.
                # #EFFLUX (REMOVE FROM DICTIONARY, ADD TO FIELD)
                # dictionaryAttrib = CompuCell.getPyAttrib(cell)
                # # efflux=3.58E-06*dictionaryAttrib[3] # for normal tissue
                # efflux=1.38E-06*dictionaryAttrib[3] # for tumor tissue
                # if efflux>0:
                #     dictionaryAttrib[3]-=efflux
                #     attrSecretor.secreteInsideCellAtBoundary(cell,efflux) # uM secretion from pixels at outer boundary of cell

#                 MEDIUM is a null pointer, not a true cell type, and cannot secrete
#                 (conversation w/CC3D developer Maciek Swat, July 2014); must add another medium-like cell type to secrete




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
class DiffusionSolverFESteeringCisplatinIV(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.simulator=_simulator
    def start(self):
        pass
    def step(self,mcs):

        ##### DRUG CONCENTRATIONS AFTER IV DELIVERY:
        # INTRAVENOUS DRUG CONCENTRATION
        # tMins=mcs/465.189 # diffusion time for one cell diameter in normal tissue
#         if 0<=mcs<cFirst5Mins:
        if 0<=mcs<Infusion15Mins:
            tMins=mcs/1207.183 # diffusion time for one cell diameter in tumor tissue
            IVtMins = 0.3725*tMins # linear fit for 15 min infusion(Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall']))
            IVxml=IVtMins            # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall'])
            self.updateXML()

# IP-POST-IV = 0
# NO DATA FOR IP-POST-IV PRIOR TO 5 MINUTES AFTER IV INFUSION (THOUGH THEY MAY HAVE MARKED 0 CONC @ 0 MINUTES AFTER INFUSION)
        elif Infusion15Mins<=mcs<cIVFirstPoint:

            IVtMins = 5.59 # highest and first data point after infusion(Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall']))
            IVxml=IVtMins            # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall'])

            # IP-post-IV CONCENTRATION
            tMins=(mcs/1207.183) - 15.0 # diffusion time for one cell diameter in tumor tissue; take away added infusion time so fit is correct; use floats

            IPpostIV =  2.366e-07*tMins**3 - 0.0001287*tMins**2 + 0.01389*tMins + 0.4209 # Casper, 1984
            IPpostIVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium']))
            IPpostIVxml=IPpostIV # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IPpostIVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium'])
            self.updateXML()
            
        elif cIVFirstPoint<=mcs<(170.862*1207.183+Infusion15Mins):        # prior to end of IV data set

#             tMins=(mcs/1207.183) - Infusion15Mins # diffusion time for one cell diameter in tumor tissue; take away added infusion time so fit is correct
            tMins=(mcs/1207.183) - 15.0 # diffusion time for one cell diameter in tumor tissue; take away added infusion time so fit is correct; use floats
            IVtMins = -1.154e-06*tMins**3 + 0.0005737*tMins**2 - 0.09922*tMins + 5.973 # Casper, 1984
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall']))
            IVxml=IVtMins             # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall'])

            # IP-post-IV CONCENTRATION
            IPpostIV =  2.366e-07*tMins**3 - 0.0001287*tMins**2 + 0.01389*tMins + 0.4209 # Casper, 1984
            IPpostIVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium']))
            IPpostIVxml=IPpostIV            # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IPpostIVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium'])
            self.updateXML()

        else:        # IP-only concentration after end of IV data set; IV = 0 in XML
        
            tMins=(mcs/1207.183) - 15.0 # diffusion time for one cell diameter in tumor tissue; take away added infusion time so fit is correct; use floats
            # IP-post-IV CONCENTRATION
            IPpostIV =  2.366e-07*tMins**3 - 0.0001287*tMins**2 + 0.01389*tMins + 0.4209 # Casper, 1984
            IPpostIVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium']))
            IPpostIVxml=IPpostIV            # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IPpostIVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium'])
            self.updateXML()

    def finish(self):
        # Finish Function gets called after the last MCS
        pass

        # # POST-INTRAVENOUS DOSING INTRAPERITONEAL CONCENTRATION -- SUGARBAKER
        # # use "if" statements when post-IV IP does not go negative when IV = 0
        #     # tHrs=mcs/465.189/60 # diffusion time for one cell diameter in normal tissue
        # # tHrs=mcs/60/1207.183 # diffusion time for one cell diameter in tumor tissue
        # if mcs < cIPpostIVZero: # x-intercept (0) of [IP]/[IV]
        #     # IV AND IP GO TO ZERO SIMULTANEOUSLY
        #     if mcs<cIVZero:

        #         # IP CONCENTRATION FROM FIT OF RATIO OF IP CONCENTRATION TO IV CONCENTRATION
        #         # calculate [IP] for when [IV]<=0
        #         # IVtHrs=(-3.338*math.log(tMins) + 16.094)/60
        #         # IVtMins=-3.338*math.log(tMins) + 16.094 # fit for patient [cisplatin IV], t=min
        #         #IPpostIV = (IVtHrs)*(0.5431*tHrs - 0.386) # IP proportion of IV fit from Sugarbaker, 1996;
                # linear fit for ratio: (0.0091*tMins - 0.386)
        # IPpostIV = (IVtMins)*( -2E-07*tMins**3 + 0.0002*tMins**2 - 0.0197*tMins + 0.9963) # goes to zero when t = 896.305 (Wolfram Alpha cubic polynomial solver)
        # # get xml IP fluid ("Medium) constant concentration value
        # IPpostIVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField'],['SecretionData'],['ConstantConcentration','Type','Medium']))
        # # set new xml value for Medium constant concentration
        # # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
        # IPpostIVxml=IPpostIV
        # self.setXMLElementValue(IPpostIVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField'],['SecretionData'],['ConstantConcentration','Type','Medium'])

            # IVtMins=6.164*math.exp( -0.02012*tMins) # fit for patient [cisplatin IV], t=min (Casper 1984)
            # IPpostIV = (IVtMins)*(-2E-07*tMins**3 + 0.0002*tMins**2 - 0.0197*tMins + 0.9963) # goes to zero when t = 896.305 (Wolfram Alpha cubic polynomial solver) #(Sugarbaker)
            # IPpostIV = (IVtMins)*(1.524e-06*tMins**3 - 0.0001449*tMins**2 + 0.009744*tMins + 0.02845)#(Old IV/IP ratio: Casper,1984)
            # IPpostIV = (IVtMins)*( -2E-07*tMins**3 + 0.0002*tMins**2 - 0.0197*tMins + 0.9963) # goes to zero when t = 896.305 (Wolfram Alpha cubic polynomial solver) #(Sugarbaker)
            # IPpostIV = (IVtMins)*(1.524e-06*tMins**3 - 0.0001449*tMins**2 + 0.009744*tMins + 0.02845) # (Old IP/IV ratioCasper 1984) 
            # IVtMins=-3.338*math.log(tMins) + 16.094 # fit for patient [cisplatin IV], t=min (Sugarbaker)

            #IVtMins=-3.338*math.log(tMins) + 16.094 # fit for patient [cisplatin IV], t=min #(Sugarbaker)
            # IVtMins=2.1996*tMins # linear fit for first 5 min (from [C]=0.0 to [C]=~11muM, t=min)#(Sugarbaker)
#             IVtMins = 0.9731*tMins # linear fit for first 5 min (Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)

#             # IP CONCENTRATION FROM FIT OF RATIO OF IP CONCENTRATION TO IV CONCENTRATION
#             # IPpostIV = (IVtMins)*( -2E-07*tMins**3 + 0.0002*tMins**2 - 0.0197*tMins + 0.9963) # goes to zero when t = 896.305 (Wolfram Alpha cubic polynomial solver) #(Sugarbaker)
#             IPpostIV = (IVtMins)*(1.524e-06*tMins**3 - 0.0001449*tMins**2 + 0.009744*tMins + 0.02845) # (Casper 1984)
#             IPpostIVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium']))
#             IPpostIVxml=IPpostIV # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
#             self.setXMLElementValue(IPpostIVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium'])
#             self.updateXML()
            #IVtMins=-3.338*math.log(tMins) + 16.094 # fit for patient [cisplatin IV], t=min #(Sugarbaker)
            # IVtMins=2.1996*tMins # linear fit for first 5 min (from [C]=0.0 to [C]=~11muM, t=min)#(Sugarbaker)
#             IVtMins = 0.9731*tMins # linear fit for first 5 min (Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)#(Casper 1984)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class DiffusionSolverFESteeringCisplatinIP(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.simulator=_simulator
    def start(self):
        pass
    def step(self,mcs):

        ##### DRUG CONCENTRATIONS AFTER IP DELIVERY:
        if 0<=mcs<Infusion15Mins:
            
            # INTRAPERITONEAL DRUG CONCENTRATION
            IPtMins=64.6338943509 #@7.386 min, first time point
            IPxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium']))
            IPxml=IPtMins # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IPxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium'])

            self.updateXML()

        elif Infusion15Mins<=mcs<cIPFirstPtPlus15:
            
            # INTRAPERITONEAL DRUG CONCENTRATION
            IPtMins=64.6338943509 #@7.386 min, first time point
            IPxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium']))
            IPxml=IPtMins # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IPxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium'])

            # IV-post-IP CONCENTRATION
            tMins=(mcs/1207.183) - 15.0 # diffusion time for one cell diameter in tumor tissue; take away added infusion time so fit is correct; use floats
            IVpostIP= -7.245e-09*tMins**4 + 2.772e-06*tMins**3 - 0.0003881*tMins**2 + 0.02126*tMins + 0.1941 # Casper, 1984
            IVpostIPxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall']))
            IVpostIPxml=IVpostIP            # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVpostIPxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall'])

            self.updateXML()

        elif cIPFirstPtPlus15<=mcs<(1207.183*180.378 + Infusion15Mins):

            tMins=(mcs/1207.183) - 15.0 # diffusion time for one cell diameter in tumor tissue; take away added infusion time so fit is correct; calculate using floats
            # INTRAPERITONEAL DRUG CONCENTRATION
            IPtMins=163.3*tMins**(-0.2859) - 28.31 # fit for patient [cisplatin IP] after first time point, t=min (Casper 1984)
            IPxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium']))
            IPxml=IPtMins
            self.setXMLElementValue(IPxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium'])

            # IV-post-IP CONCENTRATION
            IVpostIP= -7.245e-09*tMins**4 + 2.772e-06*tMins**3 - 0.0003881*tMins**2 + 0.02126*tMins + 0.1941 # Casper, 1984
            IVpostIPxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall']))
            IVpostIPxml=IVpostIP            # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVpostIPxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall'])

            self.updateXML()


        else: # IV conc = 0 after 180 mins

            tMins=(mcs/1207.183) - 15.0 # diffusion time for one cell diameter in tumor tissue; take away added infusion time so fit is correct; calculate using floats
            # INTRAPERITONEAL DRUG CONCENTRATION
            IPtMins=163.3*tMins**(-0.2859) - 28.31 # fit for patient [cisplatin IP] after first time point, t=min (Casper 1984)
            IPxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium']))
            IPxml=IPtMins
            self.setXMLElementValue(IPxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium'])

            self.updateXML()

    def finish(self):
        # Finish Function gets called after the last MCS
        pass







### EQNS FOR IP + IV DELIVERY
# class DiffusionSolverFESteeringCisplatinIV(SteppableBasePy):
#         if 0<mcs<cFirst5Mins:
#             IVtMinspIV=2.1996*tMins # linear fit for first 5 min (from [C]=0.0 to [C]=~11muM, t=min
#             IPpostIV = (IVtMins)*( -2E-07*tMins**3 + 0.0002*tMins**2 - 0.0197*tMins + 0.9963) # goes to zero when t = 896.305 (Wolfram Alpha cubic polynomial solver)
#         elif cFirst5Mins<=mcs<cIVZero: # x-intercept (drug conc. = 0, becoming negative) > mcs >= first IV data point for Cisplatin
#             IVtMins=-3.338*math.log(tMins) + 16.094 # fit for patient [cisplatin IV], t=min
#             IPpostIV = (IVtMins)*( -2E-07*tMins**3 + 0.0002*tMins**2 - 0.0197*tMins + 0.9963) # goes to zero when t = 896.305 (Wolfram Alpha cubic polynomial solver)
# class DiffusionSolverFESteeringCisplatinIP(SteppableBasePy):
#         if 0<mcs<cIPFirstPt2Hrs:
#             IPtHrs=1856.1*tHrs # linear fit for first 0.213 hrs (from [C]=0.0 to [C]=~400muM, t=hrs
#             IVpostIP = (IPtHrs)*(-0.000889*tHrs+0.06288)
#         else:
#             IPtHrs=451.51*math.exp(-0.551*tHrs) # fit for patient [cisplatin IP], t=hrs
#             IVpostIP = (IPtHrs)*(-0.000889*tHrs+0.06288)
#
#            # IPtHrs=1856.1*tHrs # linear fit for first 0.213 hrs (from [C]=0.0 to [C]=~400muM, t=hrs
#            # IPtHrs=451.51*math.exp(-0.551*tHrs) # fit for patient [cisplatin IP], t=hrs
#        # tHrs=mcs/1207.183/60  # diffusion time for one cell diameter in tumor tissue
#            # IPtHrs=1856.1*tHrs # linear fit for first 0.213 hrs (from [C]=0.0 to [C]=~400muM, t=hrs
#            IVpostIP = (IPtMins)*(1.698e-09*tMins**3 + -1.955e-06*tMins**2 +  0.0004822*tMins + 0.001835)#IP conc times the IV/IP ratio at current MCS
#








#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
class DiffusionSolverFESteeringCisplatinIPplusIV(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.simulator=_simulator
    def start(self):
        pass
    def step(self,mcs):

        ##### DRUG CONCENTRATIONS AFTER IV AND IP DELIVERY:
        tMins=mcs/1207.183 # diffusion time for one cell diameter in tumor tissue
        # tHrs=mcs/1207.183/60  # diffusion time for one cell diameter in tumor tissue


        if 0<=mcs<cFirst5Mins:
            # IVtMins=2.1996*tMins # linear fit for first 5 min (from [C]=0.0 to [C]=~11muM, t=min
            # IVtMins = 0.9731*tMins # linear fit for first 5 min (Casper 1984; from [C]=0.0 to [C]=~5.6muM, t=min)
            IPtMins=64.6338943509 #7.386 min, first time point
            # IPtHrs=1856.1*tHrs # linear fit for first 0.213 hrs (from [C]=0.0 to [C]=~400muM, t=hrs
            # IPpostIV = (IVtMins)*( -2E-07*tMins**3 + 0.0002*tMins**2 - 0.0197*tMins + 0.9963) # goes to zero when t = 896.305 (Wolfram Alpha cubic polynomial solver)
            IPpostIV = (IVtMins)*(1.524e-06*tMins**3 - 0.0001449*tMins**2 + 0.009744*tMins + 0.02845) # 
            IVpostIP = (IPtMins)*(1.698e-09*tMins**3 + -1.955e-06*tMins**2 +  0.0004822*tMins + 0.001835)

            IVAll=IVtMins+IVpostIP
            IPAll=IPtMins+IPpostIV

            #SET IV CONCENTRATION = IV + IVpostIP
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall']))
            IVxml=IVAll            # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall'])

            #SET IP CONCENTRATION = IP + IPpostIV
            IPpostIVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium']))
            IPpostIVxml=IPAll
            self.setXMLElementValue(IPpostIVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium'])


        # post-five minutes to IP fit > 0 (0.2 hrs)
        elif cFirst5Mins<=mcs<cIPFirstPt:
            # IVtMins=-3.338*math.log(tMins) + 16.094 # fit for patient [cisplatin IV], t=min
            IVtMins=6.164*math.exp( -0.02012*tMins) # fit for patient [cisplatin IV], t=min# fit for patient [cisplatin IV] (Casper, 1984), t=min
            # IPpostIV = (IVtMins)*( -2E-07*tMins**3 + 0.0002*tMins**2 - 0.0197*tMins + 0.9963) # goes to zero when t = 896.305 (Wolfram Alpha cubic polynomial solver)
            IPpostIV = (IVtMins)*(1.524e-06*tMins**3 - 0.0001449*tMins**2 + 0.009744*tMins + 0.02845)
            # IPtHrs=1856.1*tHrs # linear fit for first 0.213 hrs (from [C]=0.0 to [C]=~400muM, t=hrs
            IPtMins=64.6338943509 #7.386 min, first time point
            # IVpostIP = (IPtHrs)*(-0.000889*tHrs+0.06288)
            IVpostIP = (IPtMins)*(1.698e-09*tMins**3 + -1.955e-06*tMins**2 +  0.0004822*tMins + 0.001835)

            IVAll=IVtMins+IVpostIP
            IPAll=IPtMins+IPpostIV

            #SET IV CONCENTRATION = IV + IVpostIP
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall']))
            IVxml=IVAll            # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall'])

            #SET IP CONCENTRATION = IP + IPpostIV
            IPpostIVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium']))
            IPpostIVxml=IPAll
            self.setXMLElementValue(IPpostIVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium'])


        # ALTERNATE: elif cIPFirstPt2Hrs<=mcs<cIVZero: # x-intercept (drug conc. = 0, becoming negative) > mcs >= first IV data point for Cisplatin
        # SIM GOES TO 120 MINS, IV GOES TO 0 AT 896.305 MINS, SO NO LIMIT AT END, JUST AVOID ZERO MCS DUE TO LOG OPERATION
        else:
            # IVtMins=-3.338*math.log(tMins) + 16.094 # fit for patient [cisplatin IV], t=min
            IVtMins=6.164*math.exp( -0.02012*tMins) # fit for patient [cisplatin IV], t=min# fit for patient [cisplatin IV] (Casper, 1984), t=min
            # IPpostIV = (IVtMins)*( -2E-07*tMins**3 + 0.0002*tMins**2 - 0.0197*tMins + 0.9963) # goes to zero when t = 896.305 (Wolfram Alpha cubic polynomial solver)
            IPpostIV = (IVtMins)*(1.524e-06*tMins**3 - 0.0001449*tMins**2 + 0.009744*tMins + 0.02845)
            IPtMins=163.3*tMins**(-0.2859) - 28.31 # fit for patient [cisplatin IP] after first time point, t=min (Casper 1984)
            # IVpostIP = (IPtHrs)*(-0.000889*tHrs+0.06288)
            IVpostIP = (IPtMins)*(1.698e-09*tMins**3 + -1.955e-06*tMins**2 +  0.0004822*tMins + 0.001835)

            IVAll=IVtMins+IVpostIP
            IPAll=IPtMins+IPpostIV

            #SET IV CONCENTRATION = IV + IVpostIP
            IVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall']))
            IVxml=IVAll            # SET VARIABLE NEEDS TO BE SAME NAME (CAN BE + OR - ALSO) AS GOTTEN VARIABLE, FOR STEERING
            self.setXMLElementValue(IVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','VesselWall'])

            #SET IP CONCENTRATION = IP + IPpostIV
            IPpostIVxml=float(self.getXMLElementValue(['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium']))
            IPpostIVxml=IPAll
            self.setXMLElementValue(IPpostIVxml,['Steppable','Type','DiffusionSolverFE'],['DiffusionField','Name','Cisplatin'],['SecretionData'],['ConstantConcentration','Type','Medium'])


        self.updateXML()
    def finish(self):
        # Finish Function gets called after the last MCS
        pass



#CHANGEWITHCISPLATINSTEPPABLE
#CELL TYPES CHANGE AT DRUG IC50 THRESHOLD
class ChangeWithCisplatinSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator, _frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)
          # 2 lines following are mentioned in the manual, but are not in Listing 17,
          #    or the demo (cellsort_2D_field_modules.py)
          # self.dim=self.simulator.getPotts().getCellFieldG().getDim()
          # self.fieldName="FGF"
    def start(self):
        pass
    def step(self,mcs):
        for cell in self.cellList:
            dictionaryAttrib = CompuCell.getPyAttrib(cell)
            if cell.type==2: # green proliferating cancer cells
                if (dictionaryAttrib[3]>= cisplatinIC50):
                    cell.type=8 #QCancerGFP
            # if cell.type==3:  #red proliferating cancer cell
            #     if (dictionaryAttrib[3]>= cisplatinIC50):
            #         cell.type=8 #QCancerGFP




# ***************************** CODE TO ADD MITOSIS IN UNFROZEN CELLS (VOL PARAMS ARE OLD); 
# ***************************** ALSO UNCOMMENT CELL DICTIONARY KEYS 1 AND 2 FOR CELL AGE ABOVE
# VOLUMEPARAMSTEPPABLE
# SETS CELL TARGET VOLUMES OR COMPRESSIBILITIES
class VolumeParamSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

    def start(self):
        print "This function (VolumeParamSteppable) is called at every MCS"
        self.cellList=CellList(self.inventory)
        for cell in self.cellList:

        # VOLUME CONSTRAINTS:
            # MESOTHELIAL CELLS:
            # Radius Ave=10um(Human), Height Ave=0.7um (rabbit);
            # polygonal, will estimate as cylindrical;
            # V= Pi*h*r^2 = Pi*0.7*10^2= 219um^3
            # if not restrained to be flat (not using elongation or freezing),
            #    just use small cuboidal cells = 1x1x1 or = 2x2x2
            # visceral mesothelial cells
            if cell.type==1:
                cell.targetVolume=5
                cell.lambdaVolume=2.0
            # CANCER CELLS
            # proliferative skov3ip cancer cells in wilson lab, nude mouse mesentery:
            # R = 3.5 --> V = 179.59 --> cube edge = 5.64
            # (also: SKOV3 in vitro (mouse) = 11um in diameter;
            # vol=~4/3Pi(11/2)^3 = 696.92))
            #green PCancer
            elif cell.type==2:
                cell.targetVolume=cell.volume #initial target vol from current vol
                # cell.targetVolume=180
                cell.lambdaVolume=2000.0 # smaller than endothelial so won't eat cell
                # print "GFP cell.volume=",cell.volume
            #red PCancer
            elif cell.type==3:
                cell.targetVolume=cell.volume #initial target vol from current vol
                # cell.targetVolume=180
                cell.lambdaVolume=2000.0 # smaller than endothelial so won't eat cell
                # print "RFP cell.volume=",cell.volume
            # necrotic cancer cells
            elif cell.type==10:
                cell.targetVolume=0 # cell dissolves
                cell.lambdaVolume=2.0 # cell dissolves slowly?
            # ECM
            # average length of collagen fiber = 110, x 1um diameter = 110um^3
            elif cell.type==6:
                cell.targetVolume=cell.volume # fibers random sizes; keep init size
                cell.lambdaVolume=500.0
            # SMOOTH MUSCLE OF SMALL INTESTINE
            elif cell.type==7:
                cell.targetVolume=cell.volume # keep init size (why?)
                cell.lambdaVolume=2000000.0



# *****************************
# MITOSISDATA
# DEFINE CELL TYPES FOR MITOSIS
class MitosisData:
   def __init__(self,_MCS,_parentId,_parentType,_offspringId,_offspringType):
      self.MCS=_MCS
      self.parentId=_parentId
      self.parentType=_parentType
      self.offspringId=_offspringId
      self.offspringType=_offspringType
   def __str__(self):
      return "Mitosis time="+str(self.MCS)+"parentId="+str(self.parentId)+"offspringId="+str(self.offspringId)



# *****************************
# MITOSISSTEPPABLE
# CONDITION-DEPENDENT MITOSIS:
#   CELL AGE = DIVISION AGE,
#   CELL VOLUME = TARGET VOLUME = DOUBLE ORIGINAL VOLUME
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
        # self.cellList=CellList(self.inventory)
    def step(self,mcs):
        print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        # print "cells_to_divide =",cells_to_divide
        for cell in self.cellList:
#            dictionaryAttrib = CompuCell.getPyAttrib(cell)
            if cell.type==2:
                if dictionaryAttrib[1]>dictionaryAttrib[0] and cell.volume>=cell.targetVolume:
                #     print "high enough cisplatin to mark for division apoptosis?",dictionaryAttrib
                #     if dictionaryAttrib[2]>0 and random()<(dictionaryAttrib[2]/(2*cisplatinEC50)):
                #         print "cisplatin high enough for apoptosis"
                #         # turn cell necrotic at division time w/prob =%>EC50 if
                #         # cisplatin conc has been >=EC50 at any time
                #         cell.type==10
                #     else:
                        print "CANCGFP is old enough #############"
                        cells_to_divide.append(cell)
                        # print cells_to_divide
            if cell.type==3:
                if dictionaryAttrib[1]>dictionaryAttrib[0] and cell.volume>=cell.targetVolume:
                    # # print "DIVIDING CANCER CELL ATTRIBUTES (for programming):"
                    # # print dir(cell)
                    # if (dictionaryAttrib[2]>0) and random()<(dictionaryAttrib[2]/(2*cisplatinEC50)):
                    #     # turn cell necrotic at division time w/prob =%>EC50 if
                    #     # cisplatin conc has been >=EC50 at any time
                    #     cell.type==10
                    # else:
                        print "CANCRFP is old enough #############"
                        cells_to_divide.append(cell)
                        # print cells_to_divide
        print "cells_to_divide =",cells_to_divide
        for cell in cells_to_divide:
            # TO CHANGE MITOSIS MODE LEAVE ONE OF THE BELOW LINES UNCOMMENTED
            # self.divideCellRandomOrientation(cell)
            # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
            # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
            # print "CELL TO DIVIDE: ID",cell.id,"CELL TYPE",cell.type,"CELL TARGET VOLUME",cell.targetVolume,"CELL VOLUME",cell.volume
            self.divideCellAlongMinorAxis(cell)                               # this is a valid option
            print "CELL HAS DIVIDED: ID",cell.id,"CELL TYPE",cell.type,"TARGET VOLUME",cell.targetVolume,"VOLUME",cell.volume

    def updateAttributes(self):
        parentCell=self.mitosisSteppable.parentCell
        childCell=self.mitosisSteppable.childCell

        # # UPDATES FOR ALL CELL TYPES THAT HAVE DIVIDED
        # dictionaryAttribParentCell = CompuCell.getPyAttrib(parentCell)
        # # xc = gauss(0.22,0.03)
        # # xc = gauss(2.2,0.3)
        # xc = gauss(22,3)
        # dictionaryAttribParentCell[0:2]=[xc,0]
        # dictionaryAttribChildCell = CompuCell.getPyAttrib(childCell)
        # # xp = gauss(0.22,0.03)
        # # xp = gauss(2.2,0.3)
        # xp = gauss(22,3)
        # dictionaryAttribChildCell[0:2]=[xp,0]
        # print "dictionaryAttrib[0:2]", dictionaryAttribChildCell[0:2]
        # print "PCancerGFP AGE *****",dictionaryAttribChildCell[1]
        # # reset target volumes of child and parent to original vol. (cell growth increments set to make cells approx. 2xtargetvol @ division time 24 hrs.)
        # print "PARENT TARGET VOLUME",parentCell.targetVolume,"CHILD TARGET VOLUME",childCell.targetVolume
        # childCell.targetVolume=parentCell.targetVolume/2
        # # childCell.targetVolume=parentCell.targetVolume
        # print "PARENT TARGET VOLUME",parentCell.targetVolume,"CHILD TARGET VOLUME",childCell.targetVolume
        # parentCell.targetVolume=childCell.targetVolume
        # print "PARENT TARGET VOLUME",parentCell.targetVolume,"CHILD TARGET VOLUME",childCell.targetVolume
        # childCell.lambdaVolume=parentCell.lambdaVolume
        # if parentcell.type==2:
        #     childCell.type=2
        # elif parentcell.type==3:
        #     childCell.type=3
        # elif parentcell.type==11:
        #     childCell.type=11

        # UPDATES FOR SPECIFIC CELL TYPES
        # for cell in self.cellList:
            # parentCell=self.mitosisSteppable.parentCell
            # childCell=self.mitosisSteppable.childCell
        if parentCell.type==2:
            # print "PARENT CELL ID",parentCell.id,"PARENT CELL TYPE",parentCell.type,"PARENT CELL TARGET VOLUME",parentCell.targetVolume,"PARENT CELL VOLUME",parentCell.volume
            # parentCell=self.mitosisSteppable.parentCell
            # childCell=self.mitosisSteppable.childCell
            dictionaryAttribParentCell = CompuCell.getPyAttrib(parentCell)
            # xc = gauss(0.22,0.03)
            xc = gauss(25.5,1)
            dictionaryAttribParentCell[0:2]=[xc,0]
            print "CANCER PARENT MITOSIS TIME IS",dictionaryAttribParentCell[0]
            dictionaryAttribChildCell = CompuCell.getPyAttrib(childCell)
            # xp = gauss(0.22,0.03)
            xp = gauss(25.5,1)
            dictionaryAttribChildCell[0:2]=[xp,0]
            print "CANCER CHILD MITOSIS TIME IS",dictionaryAttribChildCell[0]
            dictionaryAttribChildCell.append(0.0) #initialize 0 concentration cisplatin in element 2 of dictionary
            # reset target volumes of child and parent to original vol. (cell growth increments set to make cells approx. 2xtargetvol @ division time 24 hrs.)
            # print "OLD PARENT TARGET VOLUME",parentCell.targetVolume,"INIT CHILD TARGET VOLUME",childCell.targetVolume
            # childCell.targetVolume=parentCell.targetVolume/2
            childCell.targetVolume=parentCell.volume
            # print "PARENT TARGET VOLUME",parentCell.targetVolume,"CHILD TARGET VOLUME",childCell.targetVolume
            parentCell.targetVolume=childCell.targetVolume
            print "NEW PARENT TARGET VOLUME",parentCell.targetVolume,"NEW CHILD TARGET VOLUME",childCell.targetVolume
            childCell.lambdaVolume=parentCell.lambdaVolume
        elif parentCell.type==3:
            # print "PARENT CELL ID",parentCell.id,"PARENT CELL TYPE",parentCell.type,"PARENT CELL TARGET VOLUME",parentCell.targetVolume,"PARENT CELL VOLUME",parentCell.volume
            # parentCell=self.mitosisSteppable.parentCell
            # childCell=self.mitosisSteppable.childCell
            dictionaryAttribParentCell = CompuCell.getPyAttrib(parentCell)
            # xc = gauss(0.22,0.03)
            xc = gauss(25.5,1)
            dictionaryAttribParentCell[0:2]=[xc,0]
            print "CANCER PARENT MITOSIS TIME IS",dictionaryAttribParentCell[0]
            dictionaryAttribChildCell = CompuCell.getPyAttrib(childCell)
            # xp = gauss(0.22,0.03)
            xp = gauss(25.5,1)
            dictionaryAttribChildCell[0:2]=[xp,0]
            print "CANCER CHILD MITOSIS TIME IS",dictionaryAttribChildCell[0]
            dictionaryAttribChildCell.append(0.0) #initialize 0 concentration cisplatin in element 2 of dictionary
            # reset target volumes of child and parent to original vol. (cell growth increments set to make cells approx. 2xtargetvol @ division time 24 hrs.)
            # print "OLD PARENT TARGET VOLUME",parentCell.targetVolume,"INIT CHILD TARGET VOLUME",childCell.targetVolume
            # childCell.targetVolume=parentCell.targetVolume/2
            childCell.targetVolume=parentCell.volume
            # print "PARENT TARGET VOLUME",parentCell.targetVolume,"CHILD TARGET VOLUME",childCell.targetVolume
            parentCell.targetVolume=childCell.targetVolume
            print "NEW PARENT TARGET VOLUME",parentCell.targetVolume,"NEW CHILD TARGET VOLUME",childCell.targetVolume
            childCell.lambdaVolume=parentCell.lambdaVolume
            # childCell.type=11
            # childCell.type=3
            # elif cell.type==2:
            #     print "CELL ID",cell.id,"CELL TYPE",cell.type,"CELL TARGET VOLUME",cell.targetVolume,"CELL VOLUME",cell.volume
            #     # parentCell=self.mitosisSteppable.parentCell
            #     # childCell=self.mitosisSteppable.childCell
            #     dictionaryAttribParentCell = CompuCell.getPyAttrib(parentCell)
            #     xc = gauss(0.22,0.03)
            #     # xc = gauss(22,3)
            #     dictionaryAttribParentCell[0:2]=[xc,0]
            #     dictionaryAttribChildCell = CompuCell.getPyAttrib(childCell)
            #     xp = gauss(0.22,0.03)
            #     # xp = gauss(22,3)
            #     dictionaryAttribChildCell[0:2]=[xp,0]
            #     # reset target volumes of child and parent to original vol. (cell growth increments set to make cells approx. 2xtargetvol @ division time 24 hrs.)
            #     print "childCell.targetVolume=",childCell.targetVolume
            #     print "parentCell.targetVolume=",parentCell.targetVolume
            #     childCell.targetVolume=parentCell.targetVolume/2
            #     # childCell.targetVolume=parentCell.targetVolume
            #     print "childCell.targetVolume=",childCell.targetVolume
            #     parentCell.targetVolume=childCell.targetVolume
            #     print "parentCell.targetVolume=",parentCell.targetVolume
            #     childCell.lambdaVolume=parentCell.lambdaVolume        # dictionaryAttrib[1]=0
                 # childCell.type=3

        if parentCell.type==2:
            childCell.type=2
        elif parentCell.type==3:
            childCell.type=3




# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& DATA RECORDING FUNCTIONS
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class PrintAllCells(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # any code in the start function runs before MCS=0
        pass
    def step(self,mcs):
        #type here the code that will run every _frequency MCS
        for cell in self.cellList:
            print "cell.type=", cell.type, ";","cell.id=",cell.id
    def finish(self):
        # Finish Function gets called after the last MCS
        pass



# INFORPRINTERSTEPPABLE
class InfoPrinterSteppable(SteppablePy):
    def __init__(self,_simulator,_frequency=60):
        SteppablePy.__init__(self,_frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

    def start(self):
        print "This function (InfoPrinterSteppable) is called once before simulation"
        # numPCancer=0
        # numIC50Cancer=0
        # totalVolumeECM=0
        # totalVolumeMesothelial=0
        # totalVolumeEndothelial=0
    def step(self,mcs):
#        print "This function is called every 60 MCS"
        for cell in self.cellList:
#            dictionaryAttrib = CompuCell.getPyAttrib(cell)
            if cell.type==2:
                numPCancer+=1
#                dictionaryAttrib = CompuCell.getPyAttrib(cell)
#                print "CELL ID=",cell.id,"TYPE=",cell.type,"TARGET VOLUME",cell.targetVolume,"VOLUME=",cell.volume,"AGE=",dictionaryAttrib[1],"MITOTIC RATE *****",dictionaryAttrib[0]
            if cell.type==8:
                numIC50Cancer+=1
#                dictionaryAttrib = CompuCell.getPyAttrib(cell)
#                print "CELL ID=",cell.id,"TYPE=",cell.type,"TARGET VOLUME",cell.targetVolume,"VOLUME=",cell.volume,"AGE=",dictionaryAttrib[1],"MITOTIC RATE *****",dictionaryAttrib[0]
            # if cell.type==1:
            #     totalVolumeMesothelial+=cell.volume
            # if cell.type==11:
                # totalVolumeEndothelial+=cell.volume
 #               dictionaryAttrib = CompuCell.getPyAttrib(cell)
 #               print "CELL ID=",cell.id,"TYPE=",cell.type,"TARGET VOLUME",cell.targetVolume,"VOLUME=",cell.volume,"AGE=",dictionaryAttrib[1],"MITOTIC RATE *****",dictionaryAttrib[0]
            # if cell.type==9:
            #     totalVolumeECM+=cell.volume
        # print "TOTAL # PCancer cells=",numPCancer
        # print "total ECM vol = ",totalVolumeECM
        # print "total Mesothelium vol = ",totalVolumeMesothelial
        # print "total new capillary (Endothelial) vol = ",totalVolumeEndothelial

            # file.write("CELL ID=%d CELL TYPE=%d volume=%d\n" %(cell.id,cell.type,cell.volume))




class CisplatinToFileSteppable(SteppableBasePy):
    # def __init__(self,_simulator,_frequency=4650): # ten minutes of diffusion, 1 cell diameter/MCS in normal tissue
    def __init__(self,_simulator,_frequency=1207): # one minutes of diffusion, 1 cell diameter/MCS in tumor tissue
    # def __init__(self,_simulator,_frequency=24240): # twenty minutes of diffusion, 1 cell diameter/MCS in tumor tissue
        SteppableBasePy.__init__(self,_simulator, _frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)
        datestring=now.strftime("%Y-%m-%d_%H_%M_%S")
        fileName=(datestring+'_'+ drug+'_' +size+'_' +vess+'_' +'Concentrations_')
#        self.openFileInSimulationOutputDirectory('%s.txt' %fileName,'w')
        # TO OUTPUT FILE TO SAME DIRECTORY AS SCREENSHOTS
        # self.file=open(CompuCellSetup.getScreenshotDirectoryName() + '%s.txt' %fileName,'w')
        self.file=open(CompuCellSetup.getScreenshotDirectoryName() + '%s.txt' %fileName,'w')
        print 'output file =',fileName
    def start(self):
        # self.file.write('CellID CellType CellVolume\n')
        self.file.write('MCS CellType CellID COMx COMy COMz Cisplatin(uM)\n')
    def step(self,mcs):
        # print "cisplatin-concentration-to-file function is called every 4650 MCS (10 real-time minutes for Cisplatin in normal tissue)"
        print "cisplatin-concentration-to-file function is called every 12070 MCS (10 real-time minutes for Cisplatin in tumor tissue)"
        for cell in self.cellList:
            comPt=CompuCell.Point3D()
            field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
            # WORKS WHEN cell vol = 1 voxel; otherwise, use comPt.x=int(round(cell.xCM/float(cell.volume))) (int(averageCOM))(see CC3D manual)
            # comPt.x=int(cell.xCM)
            # comPt.y=int(cell.yCM)
            # comPt.z=int(cell.zCM)
            comPt.x=int(round(cell.xCM/float(cell.volume))) # "divide by float = 0 error, 4-19-2013"
            comPt.y=int(round(cell.yCM/float(cell.volume)))
            comPt.z=int(round(cell.zCM/float(cell.volume)))
            cisplatin=field.get(comPt) # get concentration at center of mass
            # print 'cisplatin =',cisplatin
            # self.file.write('CELL ID=%d CELL TYPE=%d volume=%d\n' %(cell.id,cell.type,cell.volume))
            self.file.write('%d %d %d %d %d %d %f \n' %(mcs,cell.type,cell.id,comPt.x,comPt.y,comPt.z,cisplatin))
    def finish(self):
        # pass
        self.file.close() # close the file




class PlotTreatedCellsSteppable(SteppableBasePy):
# supported SVG color names at https://pythonhosted.org/ete2/reference/reference_svgcolors.html

    # def __init__(self,_simulator,_frequency=4650): # once per 10 minute approx. (1/465.189 min for Cisplatin diffusion 1 cell diam in normal tissue)
    def __init__(self,_simulator,_frequency=1207): # once per minute approx. (1/465.189 min for Cisplatin diffusion 1 cell diam in tumor tissue)
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
            # self.pW=self.addNewPlotWindow(_title='Number of cells exposed to cumulative therapeutic drug quantity',_xAxisTitle='MonteCarlo Step (MCS) (1 MCS = 10 mins)',_yAxisTitle='Tumor cells')
            # self.pW.addPlot('TreatedCells',_style='Dots',_color='Red',_size=5)
            # self.pW=self.addNewPlotWindow(_title='Untreated Cells',_xAxisTitle='MonteCarlo Step (MCS) (1 MCS = 10 mins)',_yAxisTitle='Tumor cells')
            # self.pW.addPlot('UntreatedCells',_style='Dots',_color='Blue',_size=5)
            # self.pW=self.addNewPlotWindow(_title='Intravenous (IV) Cisplatin',_xAxisTitle='MonteCarlo Step (MCS) (1 MCS = 1/465.189 mins)',_yAxisTitle='Plasma Cisplatin (micromolar)')
            self.pW=self.addNewPlotWindow(_title='Intravenous (IV) Cisplatin',_xAxisTitle='MonteCarlo Step (MCS) (1 MCS = 1/1207.183 mins)',_yAxisTitle='Plasma Cisplatin (micromolar)')
            self.pW.addPlot('Cisplatin',_style='Dots',_color='MediumSpringGreen',_size=4)
            # self.pW=self.addNewPlotWindow(_title='SKOV3 Cisplatin Accumulation',_xAxisTitle='MonteCarlo Step (MCS) (1 MCS = 10 mins)',_yAxisTitle='Plasma Cisplatin (micromolar)')
            # self.pW.addPlot('SKOV3accumPerConc',_style='Dots',_color='MediumSpringGreen',_size=5)
            # self.pW.addPlot('TreatedCells',_style='Dots',_color='red',_size=5)
            # self.pW.addPlot('Cell1Vol',_style='Steps',_color='black',_size=5)
            # self.pW.setYAxisLogScale()
            # self.pW.set_xlim([5, 15])

    def step(self,mcs):
        # diffusion time for one cell diameter
        # tMins=mcs/465.189 # in normal tissue
        # tMins=mcs/1207.183 # in tumor tissue
        tHrs=mcs/1207.183/60 # in tumor tissue
        if 0<=mcs<cIPFirstPt2Hrs: # x-intercept (drug conc. = 0, becoming negative) > mcs > first IV data point for Cisplatin
            cisplatinIPtHrs=1856.1*tHrs
            self.pW.addDataPoint("Cisplatin",mcs,cisplatinIPtHrs) # name of the data series, x, y
            self.pW.showAllPlots()
            fileName="CisplatinIPplot"+str(mcs)+".txt"
            self.pW.savePlotAsPNG(fileName,1000,1000) # here we specify size of the image
        #if cIVZero>mcs>cFirst5Mins: # x-intercept (drug conc. = 0, becoming negative) > mcs > first IV data point for Cisplatin
        else:
            cisplatinIPtHrs=451.51*math.exp(-0.551*tHrs)
            # field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
            # cisplatinIV=max(field[])
            #Skov3accumperconc = -1.2e-07*cisplatin**3 + 1.3e-05*cisplatin**2 + 0.00058*cisplatin + 0.076
            # self.pW.addDataPoint("Treated Cells",mcs,treatedCells) # name of the data series, x, y
            # self.pW.addDataPoint("UntreatedCells",mcs,untreatedCells) # name of the data series, x, y
            self.pW.addDataPoint("Cisplatin",mcs,cisplatinIPtHrs) # name of the data series, x, y
            # self.pW.addDataPoint("Cisplatin",mcs,SKOV3accumPerConc) # name of the data series, x, y
            # self.pW.addDataPoint("Cell1Vol",mcs,cell1.volume) #name of the data series, x, y
            self.pW.showAllPlots()
            # fileName="NumCellsIC50Cisplatin_Plot"+str(mcs)+".png"
            fileName="CisplatinIPplot"+str(mcs)+".png"
            self.pW.savePlotAsPNG(fileName,1000,1000) # here we specify size of the image
            # fileName="CisplatinMax_txt"+str(mcs)+".txt"
            # self.pW.savePlotAsData(fileName)
            # fileName="CisplatinIVConc_txt"+str(mcs)+".txt"
            # self.pW.savePlotAsData(fileName)



class CellAccumToFileSteppable(SteppableBasePy):
    # def __init__(self,_simulator,_frequency=4650): # ten minutes of Cisplatin diffusion, 1 cell diameter/MCS in normal tissue
    def __init__(self,_simulator,_frequency=121): # 0.1 minute of Cisplatin diffusion, 1 cell diameter/MCS in tumor tissue
    # def __init__(self,_simulator,_frequency=113): # 0.1 hours of IgG diffusion, 1 cell diameter/MCS in tumor tissue
    # def __init__(self,_simulator,_frequency=1): # for testing
        SteppableBasePy.__init__(self,_simulator, _frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

        datestring=now.strftime("%Y-%m-%d_%H_%M_%S")
        fileName=(datestring + '_' + 'Accumulated' + drug + '_' + size + '_' + vess)
        # TO OUTPUT FILE TO SAME DIRECTORY AS SCREENSHOTS
        self.file=open(CompuCellSetup.getScreenshotDirectoryName() + '%s.txt' %fileName,'w')
        print 'output file =',fileName
    def start(self):
        self.file.write('MCS cell.type CisplatinConcentrationMcrmolar CisplatinAccumulatedMcrmolar\n')
    def step(self,mcs):
        print "cell-list-to-file function is called every 113 MCS (1 real-time minutes for IgG in tumor)"
        for cell in self.cellList:
            if cell.type==2:
                dictionaryAttrib = CompuCell.getPyAttrib(cell)
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                    # WORKS WHEN cell vol = 1 voxel; changed for tiny speed-up; otherwise, use comPt.x=int(round(cell.xCM/float(cell.volume))) (int(averageCOM))(see CC3D manual
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                cisplatin=field.get(comPt) # get concentration at center of mass
                self.file.write('%d %d %f %f \n' %(mcs,cell.type,cisplatin,dictionaryAttrib[3])) # WRITE A LINE FOR EACH CELL 
            # if cell.type==3:
            #     dictionaryAttrib = CompuCell.getPyAttrib(cell)
            #     self.file.write('%d %d %d \n' %(mcs,cell.type,dictionaryAttrib[4]))
            if cell.type==8:
                dictionaryAttrib = CompuCell.getPyAttrib(cell)
                dictionaryAttrib = CompuCell.getPyAttrib(cell)
                comPt=CompuCell.Point3D()
                field=CompuCell.getConcentrationField(self.simulator,"Cisplatin")
                    # WORKS WHEN cell vol = 1 voxel; changed for tiny speed-up; otherwise, use comPt.x=int(round(cell.xCM/float(cell.volume))) (int(averageCOM))(see CC3D manual
                comPt.x=int(cell.xCM)
                comPt.y=int(cell.yCM)
                comPt.z=int(cell.zCM)
                cisplatin=field.get(comPt) # get concentration at center of mass
                self.file.write('%d %d %f %f \n' %(mcs,cell.type,cisplatin,dictionaryAttrib[3])) # WRITE A LINE FOR EACH CELL 
    def finish(self):
        # pass
        self.file.close() # close the file




# PYTHONOUTPUTTOFILESTEPPABLE
class CellListToFileSteppable(SteppableBasePy):
    # def __init__(self,_simulator,_frequency=4650): # ten minutes of diffusion, 1 cell diameter/MCS in normal tissue
    def __init__(self,_simulator,_frequency=121): # 0.1 minute of diffusion, 1 cell diameter/MCS in tumor tissue
        SteppableBasePy.__init__(self,_simulator, _frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)
        # 2 lines following are mentioned in the manual
        # file=open('CellTotals.txt','a')
        ## FILE WILL BE IN ROOT OR C DRIVE:
        ## FILENAME "CellTotals_year-month-day_hour_minute_second.txt"

        datestring=now.strftime("%Y-%m-%d_%H_%M_%S")
        fileName=(datestring+'_'+ drug+'_' +size+'_' +vess)
        # self.file=open('%s.txt' %fileName,'w')
        # self.openFileInSimulationOutputDirectory('%s.txt' %fileName,'w')
        # TO OUTPUT FILE TO SAME DIRECTORY AS SCREENSHOTS
        self.file=open(CompuCellSetup.getScreenshotDirectoryName() + '%s.txt' %fileName,'w')
        print 'output file =',fileName
    def start(self):
        #pass
        self.file.write('MCS #UntreatedCancerCells #TreatedCancerCells\n')
    def step(self,mcs):
        print "cell-list-to-file function is called every 4650 MCS (10 real-time minutes for Cisplatin)"
        numPCancer=0
        numIC50Cancer=0
        for cell in self.cellList:
            if cell.type==2:
                numPCancer+=1
            if cell.type==8:
                numIC50Cancer+=1
            # self.file.write('CELL ID=%d CELL TYPE=%d volume=%d\n' %(cell.id,cell.type,cell.volume))
        self.file.write('%d %d %d \n' %(mcs,numPCancer,numIC50Cancer)) # WRITE AFTER ALL CELLS COUNTED
    def finish(self):
        # pass
        self.file.close() # close the file



